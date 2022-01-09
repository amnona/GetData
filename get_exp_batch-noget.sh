#!/bin/bash
#SBATCH --job-name=process_exp

echo "starting process_experiment without get from SRA"
~/git/GetData/process_experiment.py -i SraRunTable.txt --deblur-path /RG/compbio/groupData/databases/deblur -p ~/bin/sratoolkit.2.9.6-1-centos_linux64/bin/ --skip-get
echo "finished"
