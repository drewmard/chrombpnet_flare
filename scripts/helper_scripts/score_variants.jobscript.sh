#!/bin/bash

source ~/.bashrc

set -e
set -u
set -o pipefail

# module load cuda/11.2
# module load cudnn/8.1

# https://stackoverflow.com/questions/34534513/calling-conda-source-activate-from-bash-script
# eval "$(conda shell.bash hook)"
micromamba activate chrombpnet2

# this is from scg; i dont think i need
# export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

echo "Live"
python "$@"

