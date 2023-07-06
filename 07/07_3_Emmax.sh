#!/usr/bin/bash


cd /share/pub/likai/Star_protocols
software/plink/plink --bfile 07_ExWAS/cohort_final --recode12 --output-missing-genotype 0 --transpose --out 07_ExWAS/cohort_final_emmax
software/emmax-kin-intel64 -v -d 10 07_ExWAS/cohort_final_emmax
software/emmax-intel64 -v -d 10 -t 07_ExWAS/cohort_final_emmax -p 07_ExWAS/cohort_final_Ave.AD.pheno -k 07_ExWAS/cohort_final_emmax.aBN.kinf -o 07_ExWAS/cohort_final_emmax_Ave.AD_out
cat 07_ExWAS/cohort_final_emmax_Ave.AD_out.ps | cut -f1 | sed 's/:/\t/g' | paste - 07_ExWAS/cohort_final_emmax_Ave.AD_out.ps > 07_ExWAS/cohort_final_emmax_Ave.AD_out.txt
