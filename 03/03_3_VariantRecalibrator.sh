#!/usr/bin/bash


software/gatk-4.0.11.0/gatk --java-options "-Xmx248g" VariantRecalibrator -R ref/ucsc.no_hap.hg19.fa \
-V 03_combine/cohort.miss0.1.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ref/hapmap_3.3.hg19.sites.vcf \
-resource:omni,known=false,training=true,truth=true,prior=12.0 ref/1000G_omni2.5.hg19.sites.vcf \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 ref/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ref/ncbi.hg19_dbsnp147.vcf -an QD -an MQ -an MQRankSum -an ReadPosRankSum \
-an FS -an SOR -mode SNP --max-gaussians 4 -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 -O 03_combine/cohort_VQSR.snp.recal \
--tranches-file 03_combine/cohort_VQSR.snp.tranches --rscript-file 03_combine/cohort_VQSR.snp.plots.R
