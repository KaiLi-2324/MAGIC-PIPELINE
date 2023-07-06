#!/usr/bin/bash


cd /share/pub/likai/Star_protocols
software/plink/plink --bfile 04_sampleqc/cohort_combine_snpQC --extract ref/indepSNP.prune.in --genome --out 04_sampleqc/cohort_combine_snpQC
