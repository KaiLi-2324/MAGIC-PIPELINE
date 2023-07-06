#!/usr/bin/bash


cd /share/pub/likai/Star_protocols
software/plink/plink --vcf 04_sampleqc/cohort_combine_snpQC_sampleQC.recode.vcf --recode --make-bed --out 05_harmonization/cohort_combine_snpQC_sampleQC
cat 05_harmonization/cohort_combine_snpQC_sampleQC.bim | awk '{print $1"\t"$1":"$4":"$6":"$5"\t"$3"\t"$4"\t"$5"\t"$6}' > 05_harmonization/cohort_combine_snpQC_sampleQC.bim1
mv 05_harmonization/cohort_combine_snpQC_sampleQC.bim1 05_harmonization/cohort_combine_snpQC_sampleQC.bim
