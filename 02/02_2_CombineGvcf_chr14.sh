#!/usr/bin/bash


cd /share2/pub/lik/lik/Star_protocols
software/gatk-4.0.11.0/gatk CombineGVCFs --java-options "-Xmx150g -Djava.io.tmpdir=temp/" -R ref/ucsc.no_hap.hg19.fa --intervals chr14:1-107349540 \
--variant 02_gvcf/cohort.vcf.list -O 02_gvcf/cohort.chr14.g.vcf
