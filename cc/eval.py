# %%
from collections import defaultdict
from collections.abc import Callable, Generator
from functools import partial
from multiprocessing import get_context
from pathlib import Path
from sys import argv
from typing import Literal
from warnings import catch_warnings

import anndata as ad
import numpy as np
import pandas as pd
from scipy.cluster.vq import kmeans2
from scipy.sparse import csc_array
from sklearn.metrics import roc_auc_score


def ls_resid_scaled(A: np.ndarray, b: np.ndarray) -> float | np.floating:
    dn = np.sum(np.square(b - np.mean(b)))

    # return positive infinity if denominator is zero so invalid entries get sorted as large values
    return (np.linalg.norm(np.linalg.lstsq(A, b)[1]) / dn) if (dn > 0) else (np.inf)


def pearsonr(x: np.ndarray, y: np.ndarray) -> float | np.floating:
    x_min_mu = x - np.mean(x)
    y_min_mu = y - np.mean(y)
    x_var = np.mean(np.square(x_min_mu))
    y_var = np.mean(np.square(y_min_mu))

    dn = np.sqrt(x_var * y_var)

    # return negative infinity if denominator is zero so invalid entries get sorted as small values
    # return abs of pearson correlation since embeddings generally have no inherent "direction"
    return (np.abs(np.mean(x_min_mu * y_min_mu)) / dn) if (dn > 0) else (-np.inf)


def get_col(X: csc_array, j: int) -> np.ndarray:
    return X[:, j].toarray()


def mk_pertconf(fg_ad: ad.AnnData, fg_pert_conf_key: str = "pert_conf") -> dict[str, object]:
    # drop delta matrix out of conf
    pert_conf = {k: v for k, v in fg_ad.uns[fg_pert_conf_key].items() if (k not in {"delta", "orig_cells", "orig_genes"})}

    pert_conf["n_cells"] = fg_ad.obs["perturbed"].sum()
    pert_conf["n_genes"] = fg_ad.var["perturbed"].sum()

    return pert_conf


def yield_per_conf(
    data_dir: Path = Path("data") / "sP_bench",
    z_ord: Literal["C", "F"] = "F",
    z_per_col: bool = False,
) -> Generator[tuple[ad.AnnData, ad.AnnData, np.ndarray, dict[str, object]]]:
    for bg_p in (data_dir / "inp" / "bg").glob("*.h5ad"):
        bg = bg_p.stem
        bg_ad = ad.read_h5ad(bg_p)
        for fg_p in (data_dir / "inp" / "fg").glob("*.h5ad"):
            fg = fg_p.stem
            fg_ad = ad.read_h5ad(fg_p)

            pert_conf = mk_pertconf(fg_ad)

            for meth_dir in (data_dir / "out").glob("*/"):
                npz_p = meth_dir / f"{bg}-{fg}.npz"

                z = np.real(np.asarray(np.load(npz_p), order=z_ord))

                if z_per_col:
                    for j in range(z.shape[1]):
                        yield (bg_ad, fg_ad, z[:, j], ({"index": j, "method": meth_dir.stem} | pert_conf))
                else:
                    yield bg_ad, fg_ad, z, ({"method": meth_dir.stem} | pert_conf)


def pearsonr_jacc(bg_ad: ad.AnnData, fg_ad: ad.AnnData, z: np.ndarray) -> float:
    X, pert_mask, z = csc_array(fg_ad.X[:]), fg_ad.var["perturbed"].to_numpy(), np.array(z[bg_ad.n_obs :])  # ty: ignore[not-subscriptable]
    which_true = np.ravel(np.argwhere(pert_mask))

    corrs = np.fromiter(
        map(partial(pearsonr, z), map(partial(get_col, X), range(X.shape[1]))),
        X.data.dtype,
        count=X.shape[1],
    )

    n = which_true.size
    nth_corr = np.partition(corrs, -n)[-n]

    n_inter = np.sum(corrs[which_true] >= nth_corr)

    return n_inter / (which_true.size + which_true.size - n_inter)


def solve_jacc(bg_ad: ad.AnnData, fg_ad: ad.AnnData, z: np.ndarray) -> float:
    X, pert_mask, z = csc_array(fg_ad.X[:]), fg_ad.var["perturbed"].to_numpy(), np.array(z[bg_ad.n_obs :, :])  # ty: ignore[not-subscriptable]
    which_true = np.ravel(np.argwhere(pert_mask))

    resids = np.fromiter(
        map(partial(ls_resid_scaled, z), map(partial(get_col, X), range(X.shape[1]))),
        X.data.dtype,
        count=X.shape[1],
    )

    n = which_true.size
    nth = np.partition(resids, n)[n]

    n_inter = np.sum(resids[which_true] <= nth)

    return n_inter / (which_true.size + which_true.size - n_inter)


def labs_auroc(bg_ad: ad.AnnData, fg_ad: ad.AnnData, z: np.ndarray) -> float:
    return roc_auc_score(fg_ad.obs["perturbed"].to_numpy(), np.array(z[bg_ad.n_obs :]))


def labs_jacc(
    bg_ad: ad.AnnData,
    fg_ad: ad.AnnData,
    z: np.ndarray,
    reps: int = 8,
    agg_fun: Callable[[np.ndarray], float] = np.mean,
) -> float:
    labs, z = fg_ad.obs["perturbed"].to_numpy(), np.array(z[bg_ad.n_obs :])

    rng = np.random.default_rng(reps)

    return agg_fun(
        np.fromiter(
            (  # compute jaccard index on both possible interpretations of kmeans labels, return max
                max(
                    np.sum(np.logical_and(labs, c)) / np.sum(np.logical_or(labs, c)),
                    np.sum(np.logical_and(labs, c_i)) / np.sum(np.logical_or(labs, c_i)),
                )
                for c, c_i in ((c, ~c) for c in (kmeans2(z, 2, minit="++", rng=rng)[1].astype(bool) for _ in range(reps)))
            ),
            float,
            reps,
        )
    )


PER_J_FNS = {
    "program recovery (pearson-r) jacc idx [slice]": pearsonr_jacc,
    "label recovery (k-means) jacc idx [slice]": labs_jacc,
    "auROC": labs_auroc,
}
FULL_FNS = {
    "label recovery [full] (k-means, jacc idx)": labs_jacc,
    "program recovery [full] (solve-resid, jacc idx)": solve_jacc,
}


def run_stats(bg_ad: ad.AnnData, fg_ad: ad.AnnData, z: np.ndarray, conf: dict[str, object]) -> list:

    out = []

    for stat_nm, fn in FULL_FNS.items():
        out.append((fn(bg_ad, fg_ad, z), conf | {"stat": stat_nm}))

    for j in range(z.shape[1]):
        z_slice = z[:, j]
        for stat_nm, fn in PER_J_FNS.items():
            out.append((fn(bg_ad, fg_ad, z_slice), {"index": j} | conf | {"stat": stat_nm}))

    return out


def proc_pmap(x: list) -> tuple[pd.DataFrame, pd.DataFrame]:
    full, per_j = defaultdict(list), defaultdict(list)

    for e in x:
        for stat, conf in e:
            o = per_j if "index" in conf else full

            o["value"].append(stat)
            for k, v in conf.items():
                o[k].append(v)

    return pd.DataFrame(full), pd.DataFrame(per_j)


# %%
if __name__ == "__main__":
    args = argv[1:]
    if len(args) > 0:
        bg_ad, fg_ad = map(ad.read_h5ad, args[:2])
        meth_name = args[2]
        out_dir = Path(args[3])
        z = np.real(np.load(out_dir / "emb.npz"))

        conf = {"method": meth_name} | mk_pertconf(fg_ad)

        with catch_warnings(action="ignore", category=RuntimeWarning, lineno=30):
            outs = run_stats(bg_ad, fg_ad, z, conf)

        full_scores, per_j_scores = proc_pmap([outs])

        full_scores.to_csv(out_dir / "full_scores.tsv", sep="\t", index=False)
        per_j_scores.to_csv(out_dir / "per_j_scores.tsv", sep="\t", index=False)

    else:
        with get_context("spawn").Pool() as p:
            outs = p.starmap(run_stats, yield_per_conf(), chunksize=32)

        full_scores, per_j_scores = proc_pmap(outs)

        # %%
        full_scores.to_csv("data/sP_bench/full_scores.tsv", sep="\t", index=False)
        per_j_scores.to_csv("data/sP_bench/per_j_scores.tsv", sep="\t", index=False)
