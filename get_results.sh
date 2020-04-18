#!/bin/bash

# $1 - project name

cd ~/Projects
cd $1

# copy the delurred biom table and the run summary
rsync -avzP -e 'ssh -o "ProxyCommand ssh amnonam@luna.cslab.openu.ac.il -W %h:%p"' amnonam@my.hpc.pub.lan:work/$1/deblur/all.biom ./
rsync -avzP -e 'ssh -o "ProxyCommand ssh amnonam@luna.cslab.openu.ac.il -W %h:%p"' amnonam@my.hpc.pub.lan:work/$1/process_experiment.log ./

# show the region
grep region ./process_experiment.log

# if we got the biom table, delete the original fasta, trim and revcomp dir
if test -f "all.biom"; then
    echo "cleaning fasta directories from server"
    ssh amnonam@luna.cslab.openu.ac.il "ssh my.hpc.pub.lan rm -r work/$1/fasta"
    ssh amnonam@luna.cslab.openu.ac.il "ssh my.hpc.pub.lan rm -r work/$1/trim"
    ssh amnonam@luna.cslab.openu.ac.il "ssh my.hpc.pub.lan rm -r work/$1/revcomp"
else
	echo "all.biom not found. run failed :("
fi
