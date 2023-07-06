#!/usr/bin/bash


software/gatk-4.3.0.0/gatk --java-options "-Xmx248g"  SelectVariants -R ref/ucsc.no_hap.hg19.fa -V 03_combine/cohort_VQSR.snp.indel.recal.vcf \
--exclude-filtered --preserve-alleles -O 03_combine/cohort_after_VQSR.vcf
