#!/usr/bin/bash


software/gatk-4.3.0.0/gatk CombineGVCFs --java-options "-Xmx150g -Djava.io.tmpdir=temp/" -R ref/ucsc.no_hap.hg19.fa --intervals chrX:1-155270560 \
--variant 02_gvcf/cohort.vcf.list -O 02_gvcf/cohort.chrX.g.vcf
