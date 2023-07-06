#!/usr/bin/bash


cd /share/pub/likai/Star_protocols
software/plink/plink --vcf 04_sampleqc/cohort_combine_snpQC.vcf --recode --make-bed --out 04_sampleqc/cohort_combine_snpQC
cat 04_sampleqc/cohort_combine_snpQC.bim | awk '{print $1"\t"$1":"$4":"$6":"$5"\t"$3"\t"$4"\t"$5"\t"$6}' > 04_sampleqc/cohort_combine_snpQC.bim1
mv 04_sampleqc/cohort_combine_snpQC.bim1 04_sampleqc/cohort_combine_snpQC.bim
