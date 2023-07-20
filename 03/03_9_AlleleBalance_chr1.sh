#!/usr/bin/bash


/share/apps/R-4.2.1/bin/Rscript scripts/tools/AlleleBalanceBySample.R -v 03_split/cohort_after_VQSR.chr1.recode.vcf -o 03_split/cohort_after_VQSR.chr1.AB.vcf -i 0.2 -I 0.8
