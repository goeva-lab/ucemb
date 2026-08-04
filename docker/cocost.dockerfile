from ghcr.io/r-hub/r-minimal/r-minimal:next
copy --from=ghcr.io/astral-sh/uv:latest /uv /bin/

run ["mkdir", "workdir"]
workdir /workdir
run ["installr", "-d", "-t", "gfortran bzip2-dev", "kernlab", "bioc::rhdf5", "bioc::anndataR", "github::WuLabMDA/CoCo-ST", "Seurat"]

run ["uv", "venv", "--managed-python", "-p", "3.13"]
run ["apk", "add", "build-base"]
run ["uv", "pip", "install", "anndata", "numpy"]
env RETICULATE_PYTHON=/workdir/.venv/bin/python

copy ./R/cocost.R /workdir/main.R
entrypoint ["R", "-q", "--no-restore", "--no-save", "-f", "/workdir/main.R", "--args"]