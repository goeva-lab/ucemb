#!/bin/env python

from pathlib import Path
from subprocess import run
from sys import argv

NEEDS_KEY = {"sp"}

if __name__ == "__main__":
    oci_exe, ssh_key = argv[1:3]
    df_dir, out_dir = map(Path, argv[3:5])

    for docker_file in df_dir.glob("*.dockerfile"):
        run(
            [
                oci_exe,
                "build",
                "--progress",
                "plain",
                "--platform",
                "linux/amd64",
                *(
                    ["--secret", f"id=ssh_key,src={ssh_key}", "--secret", f"id=ssh_key_pub,src={ssh_key}.pub"]
                    if docker_file.stem in NEEDS_KEY
                    else []
                ),
                "--tag",
                f"oci-archive:{out_dir / f'{docker_file.stem}.tar'}",
                "-f",
                str(docker_file),
                ".",
            ],
            check=True,
        )
