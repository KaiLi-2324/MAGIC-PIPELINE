#!/usr/bin/bash


software/vcftools/bin/vcftools --vcf 03_split/cohort_after_VQSR.chr9.AB.VF.vcf --max-alleles 2 --hwe 0.000001 --max-missing 0.9 \
--recode --recode-INFO-all --out 03_split/cohort_after_VQSR.chr9.AB.GQ.DP.hwe.miss.vcf

