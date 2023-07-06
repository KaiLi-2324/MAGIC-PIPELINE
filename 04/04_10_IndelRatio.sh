#!/usr/bin/bash


software/vcftools/bin/vcftools --vcf 04_sampleqc/cohort_combine_snpQC.vcf --keep-only-indels --recode --recode-INFO-all  --out 04_sampleqc/cohort_combine_snpQC.indel
software/vcftools/bin/vcftools --vcf 04_sampleqc/cohort_combine_snpQC.indel.recode.vcf --extract-FORMAT-info GT --out 04_sampleqc/cohort_combine_snpQC.indel
python3 scripts/tools/calc_indel_ratio.py --path_indel_vcf 04_sampleqc/cohort_combine_snpQC.indel.recode.vcf --path_out 04_sampleqc/cohort_indel_ratio.txt --path_gt 04_sampleqc/cohort_combine_snpQC.indel.GT.FORMAT
