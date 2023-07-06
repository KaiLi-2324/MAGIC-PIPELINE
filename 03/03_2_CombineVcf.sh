#!/usr/bin/bash


cd /share2/pub/lik/lik/Star_protocols/
software/vcftools/bin/vcf-concat \
03_filter/cohort.chr1.mis0.9.recode.vcf \
03_filter/cohort.chr2.mis0.9.recode.vcf \
03_filter/cohort.chr3.mis0.9.recode.vcf \
03_filter/cohort.chr4.mis0.9.recode.vcf \
03_filter/cohort.chr5.mis0.9.recode.vcf \
03_filter/cohort.chr6.mis0.9.recode.vcf \
03_filter/cohort.chr7.mis0.9.recode.vcf \
03_filter/cohort.chr8.mis0.9.recode.vcf \
03_filter/cohort.chr9.mis0.9.recode.vcf \
03_filter/cohort.chr10.mis0.9.recode.vcf \
03_filter/cohort.chr11.mis0.9.recode.vcf \
03_filter/cohort.chr12.mis0.9.recode.vcf \
03_filter/cohort.chr13.mis0.9.recode.vcf \
03_filter/cohort.chr14.mis0.9.recode.vcf \
03_filter/cohort.chr15.mis0.9.recode.vcf \
03_filter/cohort.chr16.mis0.9.recode.vcf \
03_filter/cohort.chr17.mis0.9.recode.vcf \
03_filter/cohort.chr18.mis0.9.recode.vcf \
03_filter/cohort.chr19.mis0.9.recode.vcf \
03_filter/cohort.chr20.mis0.9.recode.vcf \
03_filter/cohort.chr21.mis0.9.recode.vcf \
03_filter/cohort.chr22.mis0.9.recode.vcf \
03_filter/cohort.chrX.mis0.9.recode.vcf \
>03_combine/cohort.miss0.1.vcf
