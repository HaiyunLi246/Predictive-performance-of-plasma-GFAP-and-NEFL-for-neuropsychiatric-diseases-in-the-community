#GFAP, same for NEFL
#!/bin/bash
for i in {1..22};do
    echo "Processing phenotype: GFAP, chromosome: ${i}"
    /public/home/Vinceyang/wbs_data/software/plink2 \
    --bfile /public/mig_old_storage/home1/Vinceyang/yl_data/WES/unrelated_0_084/ukb_wes_chr${i}_sample_qc_final_unrelated \
    --glm hide-covar cols=chrom,pos,ax,a1freq,nobs,beta,se,orbeta,ci,tz,p \
    --geno 0.05 \
    --mind 0.05 \
    --maf 0.01 \
    --hwe 1e-6 \
    --pheno /public/home/Vinceyang/lhy_data/GFAP_NFL/ExWAS/common_var/data/GFAP_common_saige.txt \
    --pheno-col-nums 3 \
    --covar /public/home/Vinceyang/lhy_data/GFAP_NFL/ExWAS/common_var/data/GFAP_common_saige.txt \
    --covar-col-nums 4-15 \
    --covar-variance-standardize \
    --out /public/home/Vinceyang/lhy_data/GFAP_NFL/ExWAS/common_var/result/GFAP/GFAP_wes_common_single_variant_chr${i}
done
