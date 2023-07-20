#!/usr/bin/bash


software/vcftools/bin/vcftools --gzvcf 02_vcf/cohort.chr3.vcf.gz --max-missing 0.9 --recode --recode-INFO-all --out 03_filter/cohort.chr3.mis0.9

