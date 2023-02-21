#!/bin/bash
#
#SBATCH --job-name=clump_mesnps
#SBATCH --output=logs/clump_mesnps_%A_%a.out
#SBATCH --error=logs/clump_mesnps_%A_%a.err
#SBATCH --mem-per-cpu=50M    # Each task uses max 100M of memory
#SBATCH --array=1-10000%100    # Submit 10'000 tasks. Run max 100 concurrently, 5931 - veh_na, 5352 - dex_na
#SBATCH --part=pe
#SBATCH --exclude=pe6,pe5

# snp_dir=/binder/mgp/datasets/2020_DexStim_Array_Human/snps
snp_dir=/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/data/snps/final_imputed_qc_snps/filtered_196_samples
treatment=$2

# Code for meqtls calculated using genotype probabilities
# rslt_dir=/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/snp_probabilities
# cd $rslt_dir
# awk 'NR != 1 {print $2}' me-qtl_cis_result_snps_${treatment}_probabilities_beta_fdr_005.csv | sort | uniq > clumped/list_${treatment}_cpgs_for_clumping.csv

# Code for meqtls calculated using SNPs with NAs
rslt_dir=/binder/mgp/workspace/2020_DexStim_Array_Human/dex-stim-human-array/output/data/integrative/matrixEQTL/meqtls/with_missingness/
# cd $rslt_dir
# awk 'NR != 1 {print $2}' me-qtl_cis_result_snps_with_na_${treatment}_beta_fdr_005.csv | sort | uniq > clumped/list_${treatment}_cpgs_for_clumping.csv

rslt_treatment_dir=${rslt_dir}/clumped/${treatment}
mkdir -p ${rslt_treatment_dir}

SCRATCH_DIRECTORY=${rslt_dir}/scratch/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}

i=$1
n=$(($SLURM_ARRAY_TASK_ID+${i}))
cpg=$(awk -v cpg_id=$n 'NR==cpg_id {print $1}' ${rslt_dir}/clumped/list_${treatment}_cpgs_for_clumping.csv)

echo "SNP P" > ${treatment}_snp_lst_for_ld_clump_${cpg}.txt
awk -v cg=$cpg '$2 ~ cg {print $1 " " $5}' ${rslt_dir}/me-qtl_cis_result_snps_with_na_${treatment}_beta_fdr_005.csv >> ${treatment}_snp_lst_for_ld_clump_${cpg}.txt

plink --bfile $snp_dir/dex_geno_imputed --clump ${treatment}_snp_lst_for_ld_clump_${cpg}.txt --clump-r2 0.2 --clump-kb 200 --clump-p1 0.05 --clump-p2 1 --out ${treatment}_${cpg}

# echo "CpG_ID SNP" >  me-qtl_cis_${treatment}_ind_cpg_snp_associations_${cpg}.txt
awk -v cg=$cpg 'NR != 1 {if($3) print cg " " $3}' ${treatment}_${cpg}.clumped >> ${SCRATCH_DIRECTORY}/me-qtl_cis_${treatment}_ind_cpg_snp_associations_${cpg}.txt

cp ${SCRATCH_DIRECTORY}/me-qtl_cis_${treatment}_ind_cpg_snp_associations_${cpg}.txt ${rslt_treatment_dir}

# Save plink clumped results to the treatment folder
mkdir -p ${rslt_treatment_dir}/plink_results/
cp ${SCRATCH_DIRECTORY}/${treatment}_${cpg}.clumped ${rslt_treatment_dir}/plink_results

cd ${rslt_dir}
rm -rf ${SCRATCH_DIRECTORY}

# awk 'NR !=1 {print $1}' ${treatment}/me-qtl_cis_${treatment}_ind_cpg_snp_associations.txt | uniq |wc -l
# awk 'NR !=1 {print $2}' ${treatment}/me-qtl_cis_${treatment}_ind_cpg_snp_associations.txt | uniq |wc -l
# wc -l ${treatment}/me-qtl_cis_${treatment}_ind_cpg_snp_associations.txt

exit 0