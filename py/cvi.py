from pathlib import Path
from sys import argv

import anndata as ad
import numpy as np
from scvi.external import ContrastiveVI


def run(bg: ad.AnnData, fg: ad.AnnData, k: int) -> np.ndarray:
    adata, bg_idx, fg_idx = ad.concat([bg, fg]), range(bg.n_obs), range(bg.n_obs, bg.n_obs + fg.n_obs)
    ContrastiveVI.setup_anndata(adata)
    mod = ContrastiveVI(adata, n_salient_latent=k, n_background_latent=10)
    mod.train(background_indices=bg_idx, target_indices=fg_idx)
    return mod.get_latent_representation(adata, representation_kind="salient")[:, :k]


if __name__ == "__main__":
    bg = ad.read_h5ad(argv[1])
    fg = ad.read_h5ad(argv[2])
    k = int(argv[3])

    with open(Path(argv[4]) / "emb.npz", "wb") as f:
        np.save(f, run(bg, fg, k))
