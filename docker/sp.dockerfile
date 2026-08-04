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

run ["apt", "install", "-y", "git", "ssh"]
run --mount=type=secret,id=ssh_key,target=/root/.ssh/gh_key --mount=type=secret,id=ssh_key_pub,target=/root/.ssh/gh_key.pub \
    GIT_SSH_COMMAND='ssh -o StrictHostKeyChecking=no -i /root/.ssh/gh_key' uv pip install 'scanpy' 'git+ssh://git@github.com/goeva-lab/spatialpeeler.git'
copy ./py/sp.py /workdir/main.py
entrypoint ["/workdir/.venv/bin/python", "/workdir/main.py"]