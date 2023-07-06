#!/usr/bin/bash


cd /share/pub/likai/Star_protocols
software/vcftools/bin/vcftools --vcf 03_combine/cohort_after_VQSR.vcf --chr chrY --recode --recode-INFO-all --out 03_split/cohort_after_VQSR.chrY
