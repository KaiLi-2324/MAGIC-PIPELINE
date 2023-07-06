#!/usr/bin/bash


#!/usr/bin/bash


software/gatk-4.3.0.0/gatk --java-options "-Xmx248g" ApplyVQSR -R ref/ucsc.no_hap.hg19.fa -V 03_combine/cohort_VQSR.snp.recal.vcf --ts-filter-level 99.9 \
--tranches-file 03_combine/cohort_VQSR.indel.tranches --recal-file 03_combine/cohort_VQSR.indel.recal -mode INDEL -O 03_combine/cohort_VQSR.snp.indel.recal.vcf
