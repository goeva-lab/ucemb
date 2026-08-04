from pathlib import Path
from sys import argv

import anndata as ad
import numpy as np
import scanpy as sc
from contrastive import CPCA


def run(bg: ad.AnnData, fg: ad.AnnData, k: int, c_param: float) -> np.ndarray:
    for e in [fg, bg]:
        sc.pp.normalize_total(e)
        sc.pp.log1p(e)
        sc.pp.scale(e, zero_center=True)
    bg, fg = bg.X, fg.X  # ty:ignore[invalid-assignment]

    mod = CPCA(n_components=k)
    mod.fit(fg, bg, preprocess_with_pca_dim=fg.shape[1])
    # cpca returns complex-valued arrays since it uses numpy.linalg.eig, which we can discard since our inputs are all real-valued
    return np.real(mod.cpca_alpha(np.vstack([bg, fg]), alpha=c_param)[:, :k])  # ty:ignore[no-matching-overload]


if __name__ == "__main__":
    fg = ad.read_h5ad(argv[2])
    bg = ad.read_h5ad(argv[1])
    k = int(argv[3])
    c_param = float(argv[5])

    with open(Path(argv[4]) / "emb.npz", "wb") as f:
        np.save(f, run(bg, fg, k, c_param))
