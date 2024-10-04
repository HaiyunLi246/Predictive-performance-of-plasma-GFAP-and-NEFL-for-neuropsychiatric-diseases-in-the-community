#!/bin/bash

phenolist='GFAP'
for Pheno in $phenolist
do
  mkdir -p /public/home/Vinceyang/lhy_data/GFAP_NFL/GWAS/result/${Pheno}/
  cd /public/home/Vinceyang/lhy_data/GFAP_NFL/GWAS/result/${Pheno}/
  for i in {1..22}
  do
    echo "Processing phenotype: ${Pheno}, chromosome: ${i}"
    /public/home/Vinceyang/wbs_data/software/plink2 \
    --bfile /public/mig_old_storage/home1/Vinceyang/UKB_gene_v3_imp_qc/UKB_gene_v3_imp_qc_chr${i} \
    --keep /public/home/Vinceyang/wbs_data/british_id.txt \
    --exclude /public/home/Vinceyang/wbs_data/snp/snp_chr${i}.txt \
    --glm hide-covar cols=chrom,pos,ax,a1freq,nobs,beta,se,tz,p \
    --geno 0.05 \
    --hwe 1e-6 \
    --maf 0.005 \
    --vif 1000 \
    --pheno /public/home/Vinceyang/lhy_data/GFAP_NFL/GWAS/data/gwas_pheno.tsv \
	  --pheno-name ${Pheno} \
    --covar /public/home/Vinceyang/lhy_data/GFAP_NFL/GWAS/data/gwas_cov.tsv \
    --out GWAS_${Pheno}_chr${i}
  done
done
