from pathlib import Path
from sys import argv

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import spatialpeeler as sP
from sklearn.decomposition import NMF


def run(bg: ad.AnnData, fg: ad.AnnData, k: int) -> tuple[np.ndarray, pd.DataFrame]:
    for ds in (bg, fg):
        sc.pp.scale(ds, zero_center=False)
    adata = ad.concat({"bg": bg, "fg": fg}, label="cond")
    e_key = "X_nmf"
    sP_key = "sP_p_hat_k={k}"
    nmf_mod = NMF(n_components=k, init="nndsvda", max_iter=1000)
    adata.obsm[e_key] = nmf_mod.fit_transform(adata.X)
    adata.varm[e_key] = nmf_mod.components_.T
    for j in range(k):
        sP.run_fact_reg(adata, "cond", "bg", j, sP_key.format(k=j), X_key=e_key, fit_args={"maxiter": 1000})

    return (
        np.array([np.asarray(adata.obsm[sP_key.format(k=j)][sP.PAT_KEY]) for j in range(k)]).T,
        sP.diag.sweep_uns(adata, sP_key),
    )


if __name__ == "__main__":
    bg = ad.read_h5ad(argv[1])
    fg = ad.read_h5ad(argv[2])
    k = int(argv[3])
    pth = Path(argv[4])

    pats, meta = run(bg, fg, k)

    with open(pth / "emb.npz", "wb") as f:
        np.save(f, pats)

    meta.to_csv(pth / "meta.tsv", sep="\t")
