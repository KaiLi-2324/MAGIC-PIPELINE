#!/usr/bin/bash


cd /share/pub/likai/Star_protocols
software/vcftools/bin/vcftools --vcf 04_sampleqc/cohort_combine_snpQC.vcf --extract-FORMAT-info GQ --out 04_sampleqc/cohort_combine_snpQC
cat 04_sampleqc/cohort_combine_snpQC.GQ.FORMAT | sed '1d' | awk '{for(i=1;i<=NF;i++) total[i]+=$i;} END {for(i=1;i<=NF;i++) printf "%f ",total[i]/NR ;}' | sed 's/ /\n/g' > 04_sampleqc/cohort_combine_snpQC.GQ
