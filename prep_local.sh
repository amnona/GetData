#!/bin/bash

# $1 - new project name

#
cd ~/Projects

# create the new project directory
mkdir $1

cd $1

# move the SraRunTable.txt file to the new dir
mv ~/Downloads/SraRunTable.txt .

# create a tsv mapping file
cat SraRunTable.txt | tr "," "\\t" > map.txt

# copy the analysis notebook to here
cp ../analysis.ipynb .

# make remote directory
ssh amnonam@sheshet.cslab.openu.ac.il "ssh my.hpc.pub.lan mkdir work/$1"

# copy the SraRunTable to the openu server
rsync -avzP -e 'ssh -o "ProxyCommand ssh amnonam@sheshet.cslab.openu.ac.il -W %h:%p"' SraRunTable.txt amnonam@login8.hpc.pub.lan:work/$1/
# copy to sheshet
# rsync -avzP SraRunTable.txt syonat4@sheshet.cslab.openu.ac.il:SraRunTable.txt
# ssh syonat4@sheshet.cslab.openu.ac.il "rsync -avzP SraRunTable.txt my.hpc.pub.lan:work/$1/"


# run the pipeline
# ssh amnonam@shetshet.cslab.openu.ac.il "ssh my.hpc.pub.lan \"cd work/$1 ; /opt/slurm/bin/sbatch --mem-per-cpu=32000 --job-name=$1 --mail-type=end --mail-user=amnonim@gmail.com /home/amnonam/scripts/get_exp_batch.sh\""
# ssh amnonam@shetshet.cslab.openu.ac.il "ssh my.hpc.pub.lan \"cd work/$1 ; /home/amnonam/scripts/doit2.sh --mem-per-cpu=32000 --job-name=$1 --mail-type=end --mail-user=amnonim@gmail.com\""
ssh amnonam@sheshet.cslab.openu.ac.il "ssh login8.hpc.pub.lan \"cd work/$1 ; /home/amnonam/scripts/doit2.sh --job-name=$1 --mail-type=end --mail-user=amnonim@gmail.com\"" 
