#!/usr/bin/bash


cd /share/pub/likai/Star_protocols/
software/gatk-4.3.0.0/gatk --java-options "-Xmx248g" ApplyVQSR -R ref/ucsc.no_hap.hg19.fa -V 03_combine/cohort.miss0.1.vcf \
--ts-filter-level 99.9 --tranches-file 03_combine/cohort_VQSR.snp.tranches --recal-file 03_combine/cohort_VQSR.snp.recal -mode SNP \
-O 03_combine/cohort_VQSR.snp.recal.vcf
