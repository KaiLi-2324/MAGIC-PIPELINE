#!/usr/bin/bash


cd /share2/pub/lik/lik/Star_protocols
software/gatk-4.0.11.0/gatk CombineGVCFs --java-options "-Xmx150g -Djava.io.tmpdir=temp/" -R ref/ucsc.no_hap.hg19.fa --intervals chr1:1-249250621 \
--variant 02_gvcf/cohort.vcf.list -O 02_gvcf/cohort.chr1.g.vcf
