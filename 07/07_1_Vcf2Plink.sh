#!/usr/bin/bash


cd /share/pub/likai/Star_protocols
software/plink/plink --vcf 06_variantannotation/cohort_combine_snpQC_sampleQC_fix_vep_combine.vcf --recode --make-bed --out 07_ExWAS/cohort_final
cat 07_ExWAS/cohort_final.bim | awk '{print $1"\t"$1":"$4":"$6":"$5"\t"$3"\t"$4"\t"$5"\t"$6}' > 07_ExWAS/cohort_final.bim1
mv 07_ExWAS/cohort_final.bim1 07_ExWAS/cohort_final.bim

