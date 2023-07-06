#!/usr/bin/bash


software/vcftools/bin/vcftools --vcf 04_sampleqc/cohort_combine_snpQC.vcf --keep 04_sampleqc/keep_samples.txt \
--non-ref-ac-any 1 --hwe 0.000001 --max-missing 0.9 --recode --recode-INFO-all --out 04_sampleqc/cohort_combine_snpQC_sampleQC
