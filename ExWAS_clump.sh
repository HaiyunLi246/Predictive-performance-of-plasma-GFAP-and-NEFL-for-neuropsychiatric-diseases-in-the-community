# Clump the results of common variants
#!/bin/bash
cd /public/home/Vinceyang/lhy_data/GFAP_NFL/ExWAS/common_var/result
awk 'BEGIN{print "MarkerName\tAllele1\tAllele2\tEffect\tStdErr\tP-value\tDirection\tHetISq\tHetChiSq\tHetDf\tHetPVal"}''{OFS="\t"}{print $3,$4,$5,$8,$9,$11}' old_total_GFAP_wes_common.glm.linear > GFAP_ForClump.txt
sed -i '2,2d' GFAP_ForClump.txt

Pheno='GFAP'

for i in {1..22};do
/public/home/Vinceyang/dyt_data/software/plink \
--bfile /public/home2/nodecw_group/UKB_WES_data/qcstep5/unrelated_0_0884/ukb_wes_chr${i}_sample_qc_final_unrelated \
--clump ${Pheno}_ForClump.txt \
--clump-p1 1e-6 \
--clump-p2 1e-4 \
--clump-r2 0.1 \
--clump-kb 5000 \
--clump-field P-value \
--clump-snp-field MarkerName \
--out /public/home/Vinceyang/lhy_data/GFAP_NFL/ExWAS/common_var/result/${Pheno}_clump/clumped_chr.${i}.${Pheno}.wes_common.txt
done