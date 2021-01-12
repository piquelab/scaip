#!/bin/bash

set -v
set -e 

echo $PWD

##generate folder
#find ../fastq/ -name 'SCAIP*fastq.gz' | sed 's/.*\///;s/_S.*.fastq.gz//' | sort | uniq > libList.txt
#cat libList.txt | sed 's/4-/4/;s/-[Ee][tT]OH//;s/ctrl/CTRL/;s/-RNA//;s/-LPS$/-LPS-EtOH/;s/-PHA$/-PHA-EtOH/' > libList.id.txt
#paste libList.txt libList.id.txt > libList.new.txt

##find ${fastqfolder} -name '*fastq.gz' | sed 's/.*\///;s/_S1.*//' | sort | uniq > libList.txt
##find ${fastqfolder} -name '*CC*fastq.gz' | sed 's/.*\///;s/_S.*.fastq.gz//' | sort | uniq > libList.txt
#find ${fastqfolder} -name 'SCAIP6*fastq.gz' | sed 's/.*\///;s/_S.*.fastq.gz//' | sort | uniq > libList.scaip6.v3.txt

## Can use this one as an argument. 
##fastqfolder=../../fastq/
fastqfolder=${PWD}/../fastq/

kref=${PWD}/ref/



mkdir -p bus

cat libList.v3.error | \
while read sample id; 
do 
##    fastqs=`find ${fastqfolder} -name "${sample}*R[12]*fastq.gz"`
##    fastqlist=`echo ${fastqs} | tr ' ' ,`
##time kallisto bus -i ${kallistoref} -o ./bus/${sample} -x 10xv2 -t 20 ${fastqs}
    echo "#################"
    fastqs=`find ${fastqfolder} -name "${sample}_*R[12]*fastq.gz" | sort | uniq`
    fastqlist=`echo ${fastqs} | tr ' ' ' '`
    echo "#################"
    echo $sample 
    echo $fastqlist
##    echo $fastqs
    echo "cd $PWD; 
module unload python;
module load anaconda3.python;
source activate kallisto2;
mkdir -p \${TMPDIR}/${id};
cd \${TMPDIR}/${id};
time kb count --h5ad -i ${kref}/index.idx -o ${id} -g ${kref}/t2g.txt -x 10XV3 -c1 ${kref}/spliced_t2c.txt -c2 ${kref}/unspliced_t2c.txt --lamanno --filter bustools --verbose -t 8 ${fastqlist};
mv ${id} ${PWD}/bus/;" | qsub -q erprq -l nodes=1:ppn=8 -l mem=120g -N $sample
    sleep 5;
done


# kb count --h5ad -i index.idx -g t2g.txt -x 10xv2 -o SRR6470906 \
# -c1 spliced_t2c.txt -c2 unspliced_t2c.txt --lamanno --filter bustools -t 2 \
# SRR6470906_S1_L001_R1_001.fastq.gz \
# SRR6470906_S1_L001_R2_001.fastq.gz \
# SRR6470906_S1_L002_R1_001.fastq.gz \
# SRR6470906_S1_L002_R2_001.fastq.gz


## --jobmode=erprq
## qsub -I -q erprq -l nodes=1:ppn=28 -l mem=120g
