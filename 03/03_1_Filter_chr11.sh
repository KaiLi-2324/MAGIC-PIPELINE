#!/usr/bin/bash


cd /share2/pub/lik/lik/Star_protocols
software/vcftools/bin/vcftools --gzvcf 02_vcf/cohort.chr11.vcf.gz --max-missing 0.9 --recode --recode-INFO-all --out 03_filter/cohort.chr11.mis0.9

