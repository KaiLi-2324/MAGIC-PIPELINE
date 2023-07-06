#!/usr/bin/bash


software/plink/plink --bfile 07_ExWAS/cohort_final --indep-pairwise 50 10 0.1 --out 07_ExWAS/cohort_final_indep
software/plink/plink --bfile 07_ExWAS/cohort_final --extract 07_ExWAS/cohort_final_indep.prune.in --maf 0.05 --snps-only --make-bed --out 07_ExWAS/cohort_final_maf0.05
software/plink/plink --bfile 07_ExWAS/cohort_final_maf0.05 --pca --autosome --out 07_ExWAS/cohort_final_PCA
