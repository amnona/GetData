#!/bin/bash
#SBATCH --job-name=process_exp

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "starting process_experiment without get from SRA"
python3 -m GetData.process_experiment -i SraRunTable.txt --deblur-path /RG/compbio/groupData/databases/deblur -p ~/bin/sratoolkit.2.9.6-1-centos_linux64/bin/ --skip-get
echo "finished"
