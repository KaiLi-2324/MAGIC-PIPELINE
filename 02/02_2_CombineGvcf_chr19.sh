#!/usr/bin/bash


software/gatk-4.0.11.0/gatk CombineGVCFs --java-options "-Xmx150g -Djava.io.tmpdir=temp/" -R ref/ucsc.no_hap.hg19.fa --intervals chr19:1-59128983 \
--variant 02_gvcf/cohort.vcf.list -O 02_gvcf/cohort.chr19.g.vcf
