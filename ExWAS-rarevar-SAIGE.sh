#!/bin/bash
#GFAP, same for NEFL
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH -o %j.out
#SBATCH -e %j.err

pheno='GFAP'
module load miniconda3
conda activate RSAIGE
for i in {1..22};do
	Rscript /public/home/Vinceyang/wxr_data/step1_fitNULLGLMM.R   \
		--sparseGRMFile=/public/mig_old_storage/home1/Vinceyang/yl_data/WES/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
		--sparseGRMSampleIDFile=/public/mig_old_storage/home1/Vinceyang/yl_data/WES/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
		--plinkFile=/public/mig_old_storage/home1/Vinceyang/yl_data/WES/unrelated_0_084/ukb_wes_chr${i}_sample_qc_final_unrelated \
		--useSparseGRMtoFitNULL=FALSE   \
		--useSparseGRMforVarRatio=TRUE \
		--phenoFile=/public/home/Vinceyang/lhy_data/GFAP_NFL/ExWAS/rare_var_SAIGE/Step1/${pheno}/${pheno}_Caucasian.csv \
		--phenoCol=${pheno} \
		--covarColList=sex,age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
		--qCovarColList=sex  \
		--sampleIDColinphenoFile=eid \
		--isCovariateOffset=FALSE \
		--traitType=quantitative \
		--nThreads=30   \
		--isCateVarianceRatio=FALSE \
		--cateVarRatioMinMACVecExclude=20 \
    --cateVarRatioMaxMACVecInclude=50000 \
		--outputPrefix=/public/home/Vinceyang/lhy_data/GFAP_NFL/ExWAS/rare_var_SAIGE/Step1/${pheno}/${pheno}_s1_chr${i}_10PC_both	\
		--IsOverwriteVarianceRatioFile=TRUE	
done

for i in {1..22};do
  Rscript /public/home/Vinceyang/wxr_data/step2_SPAtests.R \
    --bedFile=/public/mig_old_storage/home1/Vinceyang/yl_data/WES/unrelated_0_084/ukb_wes_chr${i}_sample_qc_final_unrelated.bed \
    --bimFile=/public/mig_old_storage/home1/Vinceyang/yl_data/WES/unrelated_0_084/ukb_wes_chr${i}_sample_qc_final_unrelated.bim \
    --famFile=/public/mig_old_storage/home1/Vinceyang/yl_data/WES/unrelated_0_084/ukb_wes_chr${i}_sample_qc_final_unrelated.fam \
    --SAIGEOutputFile=/public/home/Vinceyang/lhy_data/GFAP_NFL/ExWAS/rare_var_SAIGE/Step2/${pheno}_new/${pheno}_s2_chr${i}_10PC_both.txt \
    --AlleleOrder=ref-first \
    --minMAF=0 \
    --minMAC=0.5 \
    --GMMATmodelFile=/public/home/Vinceyang/lhy_data/GFAP_NFL/ExWAS/rare_var_SAIGE/Step1/${pheno}/${pheno}_s1_chr${i}_10PC_both.rda \
    --varianceRatioFile=/public/home/Vinceyang/lhy_data/GFAP_NFL/ExWAS/rare_var_SAIGE/Step1/${pheno}/${pheno}_s1_chr${i}_10PC_both.varianceRatio.txt \
    --sparseGRMFile=/public/mig_old_storage/home1/Vinceyang/yl_data/WES/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile=/public/mig_old_storage/home1/Vinceyang/yl_data/WES/GRM/UKB_GRM_relatednessCutoff_0.05_5000_randomMarkersUsed_unrelated_2nd_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --groupFile=/public/mig_old_storage/home1/Vinceyang/yl_data/WES/lof_missense_five/SnpEff_gene_group_chr${i}.txt \
    --annotation_in_groupTest="lof,missense,missense:lof" \
    --maxMAF_in_groupTest=0.00001,0.0001,0.001,0.01 \
    --is_output_markerList_in_groupTest=TRUE \
    --LOCO=FALSE \
    --is_fastTest=TRUE
done
