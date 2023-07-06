#!/usr/bin/bash


cd /share2/pub/lik/lik/Star_protocols
software/vcftools/bin/vcftools --gzvcf 02_vcf/cohort.chrY.vcf.gz --max-missing 0.9 --recode --recode-INFO-all --out 03_filter/cohort.chrY.mis0.9

