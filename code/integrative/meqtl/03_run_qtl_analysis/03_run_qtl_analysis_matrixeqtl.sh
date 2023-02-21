#!/bin/bash

if [ $# -ne 6 ]; then
    echo "$0: usage: Not enough arguments
          First argument: treatment (veh, dex, delta)
          Second argument: source folder path 
	        Third argument: result folder path
	        Fourth argument: partition
	        Fifth argument: node
	        Sixth argument: memory in Gb"
    exit 1
fi

module load R

treatment=$1
src_dir=$2
rslt_dir=$3
partition=$4
node=$5
memory=$6

job_name=meqtl_$1
out_fn=meqtl_$1.out

type="with_na"

sbatch --job-name=$job_name --part=$partition \
  --nodelist=$partition$node --mem=$memory --output=$out_fn \
  --wrap="Rscript --vanilla /home/ahryhorzhevska/mpip/projects/dex-stim-human-array/code/integrative/meqtl/03_run_qtl_analysis/03_qtl_analysis_matrixeqtl.R $treatment $type $src_dir $rslt_dir"

head -n 1 $rslt_dir/me-qtl_cis_result_snps_${type}_${treatment}_beta.csv > $rslt_dir/me-qtl_cis_result_snps_${type}_${treatment}_beta_fdr_005.csv
awk '$6<=0.05' $rslt_dir/me-qtl_cis_result_snps_${type}_${treatment}_beta.csv >> $rslt_dir/me-qtl_cis_result_snps_${type}_${treatment}_beta_fdr_005.csv

# ./03_run_qtl_analysis_matrixeqtl.sh dex /binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/integrative/matrixEQTL/ /binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/with_missingness/ pe 7 800Gb
