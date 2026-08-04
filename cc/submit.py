#!/bin/env python
from hashlib import blake2b
from pathlib import Path
from subprocess import run
from sys import argv
from tempfile import NamedTemporaryFile
from time import time_ns

CC_CONFS = {
    "rorqual": {"host": "rorqual.alliancecan.ca", "workdir": "links/scratch/ucemb"},
    "fir": {"host": "fir.alliancecan.ca", "workdir": "scratch/ucemb"},
    "nibi": {"host": "nibi.alliancecan.ca", "workdir": "scratch/ucemb"},
}
METH_CONFS = {
    "pca": ("--time=2:0:0", "--mem=8GB"),
    "sp": ("--time=2:0:0", "--mem=8GB"),
    "concord": ("--time=2:0:0", "--mem=8GB", "--gpus-per-node=h100:1"),
    "cvi": ("--time=2:0:0", "--mem=12GB", "--gpus-per-node=h100:1"),
    "pcpca": ("--time=8:0:0", "--mem=70GB"),
    "ccpca": ("--time=30:0:0", "--mem=100GB"),
    "gcpca": ("--time=18:0:0", "--mem=128GB"),
    "cpca": ("--time=18:0:0", "--mem=128GB"),
    "cocost": ("--time=12:0:0", "--mem=300GB"),
}


def main(
    host: str,
    workdir: str,
    sponsor_acc: str,
    batch_file: Path,
    eval_file: Path,
    fgs: list[Path],
    bgs: list[Path],
    ks: list[int],
    c_params: list[float],
    run_arr: list[tuple[Path, tuple[str, ...]]],
):
    run(["rsync", "-vut", *map(str, bgs), f"{host}:{workdir}/data/bg"], check=True)
    run(["rsync", "-vut", *map(str, fgs), f"{host}:{workdir}/data/fg"], check=True)
    run(["rsync", "-vut", *(str(e[0]) for e in run_arr), f"{host}:{workdir}/oci"], check=True)
    run(["rsync", "-vut", str(batch_file), f"{host}:{workdir}/{batch_file.name}"], check=True)
    run(["rsync", "-vut", str(eval_file), f"{host}:{workdir}/{eval_file.name}"], check=True)

    jobs_specs = dict()  # noqa: C408
    for ofile, args in run_arr:
        if args not in jobs_specs:
            jobs_specs[args] = []
        jobs_specs[args].extend(
            sorted(
                ";".join(
                    (
                        f"OCI_ARCHIVE_PATH='oci/{ofile.name}'",
                        f"BG_H5AD_PATH='data/bg/{bg.name}'",
                        f"FG_H5AD_PATH='data/fg/{fg.name}'",
                        f"OUT_DIR='out/{ofile.stem}/{bg.stem}-{fg.stem}'",
                        f"EMB_DIM='{k}'",
                        f"C_PARAM='{c_param}'",
                    )
                )
                for bg in bgs
                for fg in fgs
                for k in ks
                for c_param in c_params
            )
        )

    for args, specs in jobs_specs.items():
        spec_str = "\n".join(specs)
        o = blake2b(digest_size=8)
        # mix both specs and args into hash
        for arg in sorted(args):
            o.update(arg.encode())
        o.update(spec_str.encode())
        out_fn = f"{(time_ns() / (10**6)):.0f}-{o.hexdigest()}.cmd"

        with NamedTemporaryFile("w", delete_on_close=False) as f:
            f.write(spec_str)
            f.close()
            run(["rsync", "-vut", f.name, f"{host}:{workdir}/jobs/{out_fn}"], check=True)

        run(
            [
                "ssh",
                host,
                "sbatch",
                "--parsable",
                f"--account={sponsor_acc}",
                f"--array=1-{len(specs)}",
                f"--chdir={workdir}",
                f"--export=JOB_FILE=jobs/{out_fn}",
                *args,
                f"{workdir}/{batch_file.name}",
            ],
            check=True,
        )


if __name__ == "__main__":
    cc_clust, sponsor_acc, data_dir, tar_dir, ks, c_params, rem = (
        argv[1].lower(),
        argv[2],
        *map(Path, argv[3:5]),
        [int(k) for k in argv[5].split(",")],
        [float(c) for c in argv[6].split(",")],
        argv[7:],
    )

    main(
        **(CC_CONFS[cc_clust]),
        sponsor_acc=sponsor_acc,
        batch_file=Path("cc") / "batch.sh",
        eval_file=Path("cc") / "eval.py",
        bgs=list((data_dir / "bg").glob("*.h5ad")),
        fgs=list((data_dir / "fg").glob("*.h5ad")),
        ks=ks,
        c_params=c_params,
        run_arr=[(tar_dir / f"{m}.tar", METH_CONFS[m]) for m in (rem if len(rem) > 0 else METH_CONFS.keys())],
    )
