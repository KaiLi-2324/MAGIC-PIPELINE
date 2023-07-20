#!/usr/bin/bash


software/gatk-4.0.11.0/gatk GenotypeGVCFs --java-options "-Xmx150g -Djava.io.tmpdir=temp/" \
--variant 02_gvcf/cohort.chrX.g.vcf -new-qual true -R ref/ucsc.no_hap.hg19.fa -O 02_vcf/cohort.chrX.vcf.gz

