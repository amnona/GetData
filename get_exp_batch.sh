#!/bin/bash
#SBATCH --job-name=process_exp

# parameters - all parameters are passed to process_experiment.py

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "starting process_experiment"
echo "parameters: $@"
python3 -m GetData.process_experiment -i SraRunTable.txt --deblur-path /RG/compbio/groupData/databases/deblur -p ~/bin/sratoolkit.3.0.0-centos_linux64/bin/ "$@"
echo "finished"
