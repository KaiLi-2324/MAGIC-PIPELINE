#!/usr/bin/bash


cd /share2/pub/lik/lik/Star_protocols
software/gatk-4.0.11.0/gatk GenotypeGVCFs --java-options "-Xmx150g -Djava.io.tmpdir=temp/" \
--variant 02_gvcf/cohort.chr4.g.vcf -new-qual true -R ref/ucsc.no_hap.hg19.fa -O 02_vcf/cohort.chr4.vcf.gz

