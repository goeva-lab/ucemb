from pathlib import Path
from sys import argv

import anndata as ad
import numpy as np
import scanpy as sc
from ccpca import CCPCA


def run(bg: ad.AnnData, fg: ad.AnnData, k: int, c_param: float) -> np.ndarray:
    for e in [fg, bg]:
        sc.pp.normalize_total(e)
        sc.pp.log1p(e)
        sc.pp.scale(e, zero_center=True)
    bg, fg = bg.X, fg.X  # ty:ignore[invalid-assignment]

    return CCPCA(n_components=k).fit(fg, bg, alpha=c_param).transform(np.vstack([bg, fg]))[:, :k]  # ty:ignore[no-matching-overload]


if __name__ == "__main__":
    bg = ad.read_h5ad(argv[1])
    fg = ad.read_h5ad(argv[2])
    k = int(argv[3])
    c_param = float(argv[5])

    with open(Path(argv[4]) / "emb.npz", "wb") as f:
        np.save(f, run(bg, fg, k, c_param))
