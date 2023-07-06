#!/usr/bin/bash

software/gatk-4.3.0.0/gatk --java-options "-Xmx100g" VariantFiltration -R ref/ucsc.no_hap.hg19.fa -G-filter "GQ < 20.0" -G-filter-name lowGQ \
-G-filter "DP < 10.0" -G-filter-name lowDP --set-filtered-genotype-to-no-call -V 03_split/cohort_after_VQSR.chr18.AB.vcf -O 03_split/cohort_after_VQSR.chr18.AB.VF.vcf

