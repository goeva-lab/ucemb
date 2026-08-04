from pathlib import Path
from sys import argv

import anndata as ad
import numpy as np
import scanpy as sc
from pcpca import PCPCA


def run(bg: ad.AnnData, fg: ad.AnnData, k: int, c_param: float) -> np.ndarray:
    for e in [fg, bg]:
        sc.pp.normalize_total(e)
        sc.pp.log1p(e)
        sc.pp.scale(e, zero_center=True)
    bg, fg = bg.X, fg.X  # ty:ignore[invalid-assignment]

    # need to invert order of returned tuple since PCPCA returns a tuple of (fg_emb, bg_emb)
    # and we want our matrix to be a hconcat of order [..bg_emb, ..fg_emb]
    return np.hstack(PCPCA(n_components=k, gamma=c_param).fit_transform(fg.T, bg.T)[::-1]).T[:, :k]


if __name__ == "__main__":
    bg = ad.read_h5ad(argv[1])
    fg = ad.read_h5ad(argv[2])
    k = int(argv[3])
    c_param = float(argv[5])

    with open(Path(argv[4]) / "emb.npz", "wb") as f:
        np.save(f, run(bg, fg, k, c_param))
