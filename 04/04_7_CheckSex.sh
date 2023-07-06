#!/usr/bin/bash


cd /share/pub/likai/Star_protocols
software/plink/plink --bfile 04_sampleqc/cohort_combine_snpQC --exclude ref/inversion.txt --range --indep-pairwise 50 5 0.2 --out 04_sampleqc/indepSNP
software/plink/plink --bfile 04_sampleqc/cohort_combine_snpQC --extract ref/indepSNP.prune.in --check-sex --set-hh-missing --out 04_sampleqc/cohort_combine_snpQC
