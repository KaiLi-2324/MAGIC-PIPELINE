#!/usr/bin/bash


software/vcftools/bin/vcftools --vcf 03_combine/cohort_after_VQSR.vcf --chr chr15 --recode --recode-INFO-all --out 03_split/cohort_after_VQSR.chr15
