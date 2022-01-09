#!/bin/bash

echo "startung conda env"
eval "$(conda shell.bash hook)"
conda activate calour
echo "conda env started. submitting batch"
echo $@
sbatch $@ ~/scripts/get_exp_batch-noget.sh
echo "submitted"
