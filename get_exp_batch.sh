#!/bin/bash
#SBATCH --job-name=process_exp

# parameters - all parameters are passed to process_experiment.py

echo "starting process_experiment"
echo "parameters: $@"
~/git/GetData/process_experiment.py -i SraRunTable.txt --deblur-path /RG/compbio/groupData/databases/deblur -p ~/bin/sratoolkit.2.9.6-1-centos_linux64/bin/
echo "finished"
