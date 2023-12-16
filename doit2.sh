#!/bin/bash

# parameters:
# $1 - dir name
# rest of parameters - to pass to get_exp_batch.sh

DIR_NAME=$1
shift

echo "startung conda env"
eval "$(conda shell.bash hook)"
conda activate calour
echo "conda env started. submitting batch"

echo "doit2.sh Processing $DIR_NAME"
echo "additional parameters: $@"

sbatch --mail-type=end --mail-user=amnonim@gmail.com --job-name=$DIR_NAME ~/scripts/get_exp_batch.sh $@
echo "submitted"
