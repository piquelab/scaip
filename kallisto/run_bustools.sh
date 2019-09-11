#!/bin/bash

set -v
set -e 

echo $PWD


## Can use this one as an argument. 
##fastqfolder=../../fastq/

kallistoref=./ref/cDNA_introns.idx


cat libList.txt | \
while read sample; 
do 
##    fastqs=`find ${fastqfolder} -name "${sample}*R[12]*fastq.gz"`
##    fastqlist=`echo ${fastqs} | tr ' ' ,`
##time kallisto bus -i ${kallistoref} -o ./bus/${sample} -x 10xv2 -t 20 ${fastqs}
    echo "#################"
    echo $sample 
#    echo $fastqs
    echo "cd $PWD; 
module unload python;
module load bustools anaconda3.python;
##module load kallisto
source activate kallisto;

cd bus/${sample};

mkdir cDNA_capture/ introns_capture/ spliced/ unspliced/ tmp/ ;
bustools correct -w ../../ref/10xv2_whitelist.txt -p output.bus | bustools sort -o output.correct.sort.bus -t 4 - ;
bustools capture -o cDNA_capture/ -c ../../ref/cDNA_transcripts.to_capture.txt -e matrix.ec -t transcripts.txt output.correct.sort.bus ;
bustools capture -o introns_capture/ -c ../../ref/introns_transcripts.to_capture.txt -e matrix.ec -t transcripts.txt output.correct.sort.bus ;
bustools count -o unspliced/u -g ../../ref/cDNA_introns.t2g.txt -e cDNA_capture/split.ec -t transcripts.txt --genecounts cDNA_capture/split.bus ;
bustools count -o spliced/s -g ../../ref/cDNA_introns.t2g.txt -e introns_capture/split.ec -t transcripts.txt --genecounts introns_capture/split.bus ;" | qsub -q erprq -l nodes=1:ppn=6 -l mem=50g -N $sample
    sleep 0.5;
done



## --jobmode=erprq
## qsub -I -q erprq -l nodes=1:ppn=28 -l mem=120g
