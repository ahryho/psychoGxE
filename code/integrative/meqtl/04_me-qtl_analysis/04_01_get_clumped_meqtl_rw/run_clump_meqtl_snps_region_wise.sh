#!/bin/bash

treatment="dex"

for i in 0 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000 110000 120000 130000
do
  sbatch ./clump_meqtl_snps_region_wise.sh $i $treatment
done

# find clumped/${treatment}/. -type f -name "*.txt" -exec cat {} + > clumped/me-qtl_cis_${treatment}_snp_probabilities_beta_clumped_cpg_snp_associations.txt
# find clumped/${treatment}/. -type f -name "*.txt" -exec cat {} + > clumped/me-qtl_cis_${treatment}_snps_with_na_beta_clumped_cpg_snp_associations.txt
