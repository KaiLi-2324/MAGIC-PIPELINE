#!/usr/bin/bash


cd /share/pub/likai/Star_protocols/
software/bcftools-1.9/bin/bcftools stats -s - 04_sampleqc/cohort_combine_snpQC.vcf > 04_sampleqc/cohort_combine_snpQC_stats
cat 04_sampleqc/cohort_combine_snpQC_stats | grep "PSC" | grep -v "#" | awk '{print $3"\t"$7"\t"$8}' > 04_sampleqc/cohort_combine_snpQC_stats.tstv
cat 04_sampleqc/cohort_combine_snpQC_stats | grep "PSC" | grep -v "#" | awk '{print $3"\t"$5"\t"$6}' > 04_sampleqc/cohort_combine_snpQC_stats.homhet
cat 04_sampleqc/cohort_combine_snpQC_stats | grep "PSC" | grep -v "#" | awk '{print $3"\t"$5+$6"\t"$9}'> 04_sampleqc/cohort_combine_snpQC_stats.indelsnp
cat 04_sampleqc/cohort_combine_snpQC_stats | grep "PSC" | grep -v "#" | awk '{print $3"\t"$11}' > 04_sampleqc/cohort_combine_snpQC_stats.singleton
