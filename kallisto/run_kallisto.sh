#!/bin/bash

set -v
set -e 

echo $PWD


## Can use this one as an argument. 
##fastqfolder=../../fastq/
fastqfolder=../fastq/

kallistoref=./ref/cDNA_introns.idx

#find ${fastqfolder} -name 'SCAIP*fastq.gz' | sed 's/.*\///;s/_S.*.fastq.gz//' | sort | uniq > libList.txt


cat V3list.txt | \
while read sample; 
do 
##    fastqs=`find ${fastqfolder} -name "${sample}*R[12]*fastq.gz"`
##    fastqlist=`echo ${fastqs} | tr ' ' ,`
##time kallisto bus -i ${kallistoref} -o ./bus/${sample} -x 10xv2 -t 20 ${fastqs}
    echo "#################"
    fastqs=`find ${fastqfolder} -name "${sample}*R[12]*fastq.gz" | sort | uniq`
    fastqlist=`echo ${fastqs} | tr ' ' ' '`
    echo "#################"
    echo $sample 
    echo $fastqlist
##    echo $fastqs
    echo "cd $PWD; 
module unload python;
module load bustools anaconda3.python;
##module load kallisto
source activate kallisto 
time kallisto bus -i ${kallistoref} -o ./bus/${sample} -x 10xv3 -t 12 ${fastqlist}" | qsub -q erprq -l nodes=1:ppn=28 -l mem=120g -N $sample
    sleep 0.5;
done



## --jobmode=erprq
## qsub -I -q erprq -l nodes=1:ppn=28 -l mem=120g
