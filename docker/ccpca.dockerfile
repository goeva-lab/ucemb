from ghcr.io/astral-sh/uv:trixie-slim
run ["apt", "update"]
run ["apt", "install", "-y", "build-essential", "cmake"]
run ["apt", "install", "-y", "curl", "lsb-release", "gnupg"]
run ["curl", "https://apt.llvm.org/llvm.sh", "-o", "/llvm.sh"]
run ["bash", "/llvm.sh", "20"]
run ["mkdir", "workdir"]
workdir /workdir
run ["uv", "venv", "--managed-python", "-p", "3.13"]
run ["uv", "pip", "install", "anndata"]

run ["uv", "pip", "install", "scanpy", "ccpca"]
copy ./py/ccpca.py /workdir/main.py
entrypoint ["/workdir/.venv/bin/python", "/workdir/main.py"]