#!/bin/bash

# $1 - new project name
# $2 - other parameters to pass to process_experiment.py (typically --exp-type its (if needed) or 16s (default)
#
cd ~/Projects

# store the parameter $1 in a new parameter DIR_NAME
DIR_NAME=$1

# get the rest of the parameters to pass
shift

# create the new project directory
mkdir $DIR_NAME

cd $DIR_NAME

# move the SraRunTable.txt file to the new dir
# if the file SraTable.txt is in the Downloads directory copy it to the current directory. otherwise, move the file SraTable.csv to the current directory as SraRunTable.txt
# in any case, it is assumed to be a csv file
if [ -f ~/Downloads/SraRunTable.txt ]; then
    mv ~/Downloads/SraRunTable.txt .
else
    mv ~/Downloads/SraRunTable.csv SraRunTable.txt
fi

# create a tsv mapping file
cat SraRunTable.txt | tr "," "\\t" > map.txt

# copy the analysis notebook to here
cp ../analysis.ipynb .

# make remote directory
ssh amnonam@sheshet.cslab.openu.ac.il "ssh login8.hpc.pub.lan mkdir work/$DIR_NAME"
# ssh amnonam@sheshet.cslab.openu.ac.il "ssh my.hpc.pub.lan mkdir work/$DIR_NAME"

# copy the SraRunTable to the openu server
rsync -avzP -e 'ssh -o "ProxyCommand ssh amnonam@sheshet.cslab.openu.ac.il -W %h:%p"' SraRunTable.txt amnonam@login8.hpc.pub.lan:work/$DIR_NAME/
# copy to sheshet
# rsync -avzP SraRunTable.txt syonat4@sheshet.cslab.openu.ac.il:SraRunTable.txt
# ssh syonat4@sheshet.cslab.openu.ac.il "rsync -avzP SraRunTable.txt my.hpc.pub.lan:work/$1/"


# run the pipeline
# ssh amnonam@shetshet.cslab.openu.ac.il "ssh my.hpc.pub.lan \"cd work/$1 ; /opt/slurm/bin/sbatch --mem-per-cpu=32000 --job-name=$1 --mail-type=end --mail-user=amnonim@gmail.com /home/amnonam/scripts/get_exp_batch.sh\""
# ssh amnonam@shetshet.cslab.openu.ac.il "ssh my.hpc.pub.lan \"cd work/$1 ; /home/amnonam/scripts/doit2.sh --mem-per-cpu=32000 --job-name=$1 --mail-type=end --mail-user=amnonim@gmail.com\""

# ssh amnonam@sheshet.cslab.openu.ac.il "ssh login8.hpc.pub.lan \"cd work/$DIR_NAME ; /home/amnonam/scripts/doit2.sh --job-name=$DIR_NAME --mail-type=end --mail-user=amnonim@gmail.com\"" 
ssh amnonam@sheshet.cslab.openu.ac.il "ssh login8.hpc.pub.lan \"cd work/$DIR_NAME ; /home/amnonam/scripts/doit2.sh $DIR_NAME $@\"" 
