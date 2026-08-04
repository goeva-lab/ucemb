from collections.abc import Iterator
from functools import partial
from itertools import product, starmap
from pathlib import Path
from sys import argv

import anndata as ad
import numpy as np
from plotnine import *
from scipy.sparse import csr_array
from sklearn.decomposition import NMF


def pert(
    adata: ad.AnnData,
    coords: np.ndarray,
    rng: np.random.Generator,
    gene_idxs: np.ndarray,
    cents: tuple[float, float],
    rad: float,
    frac: float,
    noise_mean: float,
) -> ad.AnnData:
    assert all((0 <= e <= 1) for e in (*cents, frac))
    assert len(gene_idxs.shape) == 1
    assert coords.shape == (adata.n_obs, 2)

    # re-compute centroids
    coords_min = np.min(coords, axis=0)
    cents = coords_min + ((np.max(coords, axis=0) - coords_min) * cents)

    dist = np.sqrt(np.sum(np.square(coords - cents), axis=1))
    pert_mask_cells = dist <= rad

    if frac < 1:
        pert_mask_cells[pert_mask_cells] = rng.choice([True, False], size=np.sum(pert_mask_cells), p=[frac, 1 - frac])

    cell_idxs = np.argwhere(pert_mask_cells).ravel()

    noise_arr = rng.poisson(noise_mean, size=(cell_idxs.size, gene_idxs.size))
    x_idx, y_idx = np.repeat(cell_idxs, gene_idxs.size), np.tile(gene_idxs, cell_idxs.size)

    # remove any cells/genes from indices w/o any perturbation
    nz_mask = noise_arr != 0
    cell_idxs_s = cell_idxs[np.sum(nz_mask, axis=1) > 0]
    gene_idxs_s = gene_idxs[np.sum(nz_mask, axis=0) > 0]

    # unravel + drop zero entries prior to sparse array construction
    noise_arr, nz_mask = np.ravel(noise_arr), np.ravel(nz_mask)
    noise_arr, x_idx, y_idx = noise_arr[nz_mask], x_idx[nz_mask], y_idx[nz_mask]

    delta = csr_array(
        (noise_arr, (x_idx, y_idx)),
        shape=adata.shape,
        dtype=adata.X.dtype,  # ty: ignore[unresolved-attribute]
    )

    # sanity check
    for ax, idx in [(1, cell_idxs_s), (0, gene_idxs_s)]:
        assert all(np.sort(np.ravel(np.argwhere(delta.sum(axis=ax) > 0))) == np.sort(idx))

    pert_adata = adata.copy()

    pert_adata.raw = pert_adata
    pert_adata.X = csr_array(adata.X + delta)
    pert_adata.uns["pert_conf"] = {
        "x_cent": cents[0],
        "y_cent": cents[1],
        "rad": rad,
        "frac": frac,
        "λ": noise_mean,
        "delta": delta,
        "orig_cells": cell_idxs,
        "orig_genes": gene_idxs,
    }
    for df, idx in [(pert_adata.obs, cell_idxs_s), (pert_adata.var, gene_idxs_s)]:
        df["perturbed"] = False
        df.iloc[idx, df.columns.get_loc("perturbed")] = True  # ty: ignore[invalid-assignment]

    return pert_adata


def gen_ds(
    inp: ad.AnnData,
    seed: int,
    top_ns: list[int] | int,
    cents: list[tuple[float, float]] | tuple[float, float],
    rads: list[float] | float,
    fracs: list[float] | float,
    noise_means: list[float] | float,
) -> tuple[ad.AnnData, Iterator[ad.AnnData]]:
    # allow for singleton inputs
    if not isinstance(top_ns, list):
        top_ns = [top_ns]
    if not isinstance(cents, list):
        cents = [cents]
    if not isinstance(rads, list):
        rads = [rads]
    if not isinstance(fracs, list):
        fracs = [fracs]
    if not isinstance(noise_means, list):
        noise_means = [noise_means]

    rng = np.random.default_rng(seed)

    y = inp.obsm["spatial"][:, 1]
    y_split = (np.min(y) + np.max(y)) / 2
    spl_mask = y > y_split

    bot = inp[spl_mask, :].copy()
    # flip y around to mirror bottom
    bot_y = bot.obsm["spatial"][:, 1]
    bot.obsm["spatial"][:, 1] = np.max(bot_y) - (bot_y - np.min(bot_y))

    gene_ordr = np.argsort(
        NMF(n_components=15, init="nndsvda", random_state=rng.integers(np.pow(2, 32) - 1), max_iter=1000)
        .fit(inp.X)
        .components_.T,
        axis=0,
    )

    return inp[~spl_mask, :].copy(), starmap(
        partial(pert, bot, np.asarray(bot.obsm["spatial"]), rng),
        product((np.unique(gene_ordr[-top_n:, :]) for top_n in top_ns), cents, rads, fracs, noise_means),
    )


if __name__ == "__main__":
    if len(argv[1:]) > 0:
        adata_p = Path(argv[1])
        adata = ad.read_h5ad(adata_p)
        # we don't need any of the stored metadata for the benchmark
        del adata.obs

        base_p = Path(argv[2])
        seed = int(argv[3])

        top_ns = [int(x) for x in argv[4].split(",")]
        cents = [(float(x), float(y)) for x, y in (e.split(":", maxsplit=1) for e in argv[5].split(","))]

        rads, fracs, noise_means = ([float(x) for x in e.split(",")] for e in argv[6:9])

        start_idx = int(argv[9]) if (len(argv) > 9) else 1

        bg, fgs = gen_ds(adata, seed, top_ns, cents, rads, fracs, noise_means)
        bg.write_h5ad((base_p / "bg" / adata_p.stem).with_suffix(".h5ad"))
        for i, fg in enumerate(fgs, start_idx):
            fg.write_h5ad((base_p / "fg" / f"{adata_p.stem}_{i}").with_suffix(".h5ad"))
