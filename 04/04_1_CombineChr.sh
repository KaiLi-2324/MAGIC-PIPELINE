#!/usr/bin/bash


cd /share/pub/likai/Star_protocols
software/vcftools/bin/vcf-concat \
03_split/cohort_after_VQSR.chr1.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr2.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr3.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr4.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr5.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr6.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr7.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr8.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr9.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr10.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr11.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr12.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr13.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr14.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr15.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr16.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr17.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr18.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr19.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr20.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr21.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chr22.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
03_split/cohort_after_VQSR.chrX.AB.GQ.DP.hwe.miss.vcf.recode.vcf \
>04_sampleqc/cohort_combine_snpQC.vcf
