#!/bin/bash
#SBATCH -q primary
#SBATCH --mem=500G
#SBATCH --time=3-10:00:00
#SBATCH -N 1-1
#SBATCH -n 1

cd $PWD

module unload python
module load anaconda3.python

set NUMEXPR_MAX_THREADS=1

../MetaXcan/software/SPrediXcan.py \
  --model_db_path ../GTEx_v8/elastic_net_models/en_Spleen.db \
  --covariance ../GTEx_v8/gtex_v8_expression_elastic_net_snp_smultixcan_covariance.txt.gz \
  --gwas_file HanY_prePMID_asthma_UKBB.txt.gz \
  --snp_column SNP \
  --effect_allele_column EA \
  --non_effect_allele_column NEA \
  --beta_column OR \
  --pvalue_column P \
  --output_file results/Spleen.csv




