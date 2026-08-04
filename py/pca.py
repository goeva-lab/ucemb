from pathlib import Path
from sys import argv

import anndata as ad
import numpy as np
import scanpy as sc


def run(bg: ad.AnnData, fg: ad.AnnData, k: int) -> np.ndarray:
    for e in [fg, bg]:
        sc.pp.normalize_total(e)
        sc.pp.log1p(e)
        sc.pp.scale(e, zero_center=False)

    adata = ad.concat([fg, bg])
    sc.pp.pca(adata, n_comps=k)
    return adata.obsm["X_pca"][:, :k]


if __name__ == "__main__":
    bg = ad.read_h5ad(argv[1])
    fg = ad.read_h5ad(argv[2])
    k = int(argv[3])

    with open(Path(argv[4]) / "emb.npz", "wb") as f:
        np.save(f, run(bg, fg, k))
