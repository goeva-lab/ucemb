#!/bin/bash

#SBATCH -o slurm/%A_%a.log
#SBATCH -e slurm/%A_%a.log

# required env vars:
# JOB_FILE: path of job file containing env vars declarations per line

module load StdEnv/2023 apptainer
set -euxo pipefail

mkdir "${SLURM_TMPDIR}/out"
mkdir "${SLURM_TMPDIR}/cache"

eval "$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${JOB_FILE}" | tr -d '\n')"
METHOD_NAME="$(basename "${OCI_ARCHIVE_PATH}" '.tar')"

apptainer build --disable-cache "${SLURM_TMPDIR}/img.sif" "oci-archive://${OCI_ARCHIVE_PATH}"

# use zsh for more advanced time functionality
METHOD_NAME="${METHOD_NAME}" BG_H5AD_PATH="${BG_H5AD_PATH}" FG_H5AD_PATH="${FG_H5AD_PATH}" EMB_DIM="${EMB_DIM}" C_PARAM="${C_PARAM}" zsh - 2> "${SLURM_TMPDIR}/out/run.log" <<'EOF' 
set -euxo pipefail

TIMEFMT="::::[[ ZSH TIME ]]:::: METHOD_NAME=${METHOD_NAME} | BG_H5AD_PATH=${BG_H5AD_PATH} | FG_H5AD_PATH=${FG_H5AD_PATH} | EMB_DIM=${EMB_DIM} | C_PARAM=${C_PARAM} | time=%E | pct=%P | mem=%M ::::[[ ZSH TIME ]]::::"
time apptainer run --nv --writable-tmpfs -C \
    -B "${BG_H5AD_PATH}:/workdir/bg.h5ad" \
    -B "${FG_H5AD_PATH}:/workdir/fg.h5ad" \
    -B "${SLURM_TMPDIR}/out:/workdir/out" \
    -B "${SLURM_TMPDIR}/cache:/tmp/cache" \
    --cwd '/workdir/' \
    --env 'NUMBA_CACHE_DIR=/tmp/cache/numba,TORCH_HOME=/tmp/cache/torch' \
    "${SLURM_TMPDIR}/img.sif" \
    '/workdir/bg.h5ad' '/workdir/fg.h5ad' "${EMB_DIM}" "/workdir/out/" "${C_PARAM}"
EOF

.venv/bin/python ./eval.py "${BG_H5AD_PATH}" "${FG_H5AD_PATH}" "${METHOD_NAME}" "${SLURM_TMPDIR}/out"

mkdir -p "${OUT_DIR}"
cp "${SLURM_TMPDIR}/out/"* "${OUT_DIR}"
