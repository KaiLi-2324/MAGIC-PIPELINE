#!/usr/bin/bash


software/gatk-4.3.0.0/gatk --java-options "-Xmx248g"  VariantRecalibrator -R ref/ucsc.no_hap.hg19.fa -V 03_combine/cohort_VQSR.snp.recal.vcf \
-resource:mills,known=true,training=true,truth=true,prior=12.0 ref/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ref/ncbi.hg19_dbsnp147.vcf -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode INDEL \
--max-gaussians 4 -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -O 03_combine/cohort_VQSR.indel.recal --tranches-file 03_combine/cohort_VQSR.indel.tranches \
--rscript-file 03_combine/cohort_VQSR.indel.plots.R
