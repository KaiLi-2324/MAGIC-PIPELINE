#!/usr/bin/bash


export _JAVA_OPTIONS=-Djava.io.tmpdir=temp/

software/gatk-4.0.11.0/gatk --java-options "-Xmx8G" HaplotypeCaller -R ref/ucsc.no_hap.hg19.fa -I 01_realign/19BY0125.bam -L ref/Design_V2.hg19.bed -ERC BP_RESOLUTION \
-O 01_vcf/19BY0125_HaplotypeCaller.g.vcf && \
software/bgzip 01_vcf/19BY0125_HaplotypeCaller.g.vcf && software/tabix -p vcf 01_vcf/19BY0125_HaplotypeCaller.g.vcf.gz
