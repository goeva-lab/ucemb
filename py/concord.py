from pathlib import Path
from sys import argv

import anndata as ad
import concord as ccd
import numpy as np
import scanpy as sc


def run(bg: ad.AnnData, fg: ad.AnnData, k: int) -> np.ndarray:
    adata = ad.concat({"bg": bg, "fg": fg}, label="domain")
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    ccd.Concord(adata=adata, latent_dim=k, domain_key="domain", save_dir="/tmp/cache/concord").fit_transform()
    return adata.obsm["Concord"]


if __name__ == "__main__":
    bg = ad.read_h5ad(argv[1])
    fg = ad.read_h5ad(argv[2])
    k = int(argv[3])
    c_param = float(argv[5])

    with open(Path(argv[4]) / "emb.npz", "wb") as f:
        np.save(f, run(bg, fg, k))
