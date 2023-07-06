#!/usr/bin/bash


cd /share/pub/likai/Star_protocols
source /share/pub/likai/Software/anaconda3/bin/activate /share/pub/likai/Software/anaconda3/envs/bio/
vcffixup 04_sampleqc/cohort_combine_snpQC_sampleQC.chr14.recode.vcf > 04_sampleqc/cohort_combine_snpQC_sampleQC_fix.chr14.recode.vcf

