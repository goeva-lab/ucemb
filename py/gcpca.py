from pathlib import Path
from sys import argv

import anndata as ad
import numpy as np
from generalized_contrastive_PCA import gcPCA


def run(bg: ad.AnnData, fg: ad.AnnData, k: int, c_param: float) -> np.ndarray:
    mod = gcPCA(Ncalc=k, alpha=c_param, method="v4.1")
    mod.fit_transform(bg.X.toarray(), fg.X.toarray())  # ty:ignore[unresolved-attribute]
    return np.vstack([mod.Ra_transformed_[:, :k], mod.Rb_transformed_[:, :k]])


if __name__ == "__main__":
    fg = ad.read_h5ad(argv[2])
    bg = ad.read_h5ad(argv[1])
    k = int(argv[3])
    c_param = float(argv[5])

    with open(Path(argv[4]) / "emb.npz", "wb") as f:
        np.save(f, run(bg, fg, k, c_param))
