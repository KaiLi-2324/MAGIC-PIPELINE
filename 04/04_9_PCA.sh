#!/usr/bin/bash


software/plink/plink --bfile ref/Merge --geno 0.1 --mind 0.1 --maf 0.05 --allow-no-sex --make-bed --out 04_sampleqc/1kG_PCA
software/plink/plink --bfile 04_sampleqc/cohort_combine_snpQC --geno 0.1 --mind 0.1 --maf 0.05 --make-bed --out 04_sampleqc/cohort_combine_snpQC_maf0.05

awk '{print$2}' 04_sampleqc/cohort_combine_snpQC_maf0.05.bim > 04_sampleqc/cohort_combine_snpQC_maf0.05.txt
software/plink/plink --bfile 04_sampleqc/1kG_PCA --extract 04_sampleqc/cohort_combine_snpQC_maf0.05.txt --make-bed --out 04_sampleqc/1kG_PCA1
awk '{print $2}' 04_sampleqc/1kG_PCA1.bim > 04_sampleqc/1kG_PCA1.txt
software/plink/plink --bfile 04_sampleqc/cohort_combine_snpQC_maf0.05 --extract 04_sampleqc/1kG_PCA1.txt --recode --make-bed --out 04_sampleqc/cohort_combine_snpQC_maf0.05_1kg
awk '{print $2,$4}' 04_sampleqc/cohort_combine_snpQC_maf0.05_1kg.map > 04_sampleqc/buildhapmap.txt
software/plink/plink --bfile 04_sampleqc/1kG_PCA1 --update-map 04_sampleqc/buildhapmap.txt --make-bed --out 04_sampleqc/1kG_PCA2
software/plink/plink --bfile 04_sampleqc/cohort_combine_snpQC_maf0.05_1kg --bmerge 04_sampleqc/1kG_PCA2.bed 04_sampleqc/1kG_PCA2.bim 04_sampleqc/1kG_PCA2.fam --allow-no-sex --make-bed --out 04_sampleqc/cohort_combine_snpQC_maf0.05_1kg_pca
software/plink/plink --bfile 04_sampleqc/cohort_combine_snpQC_maf0.05_1kg_pca --pca --out 04_sampleqc/cohort_combine_snpQC_maf0.05_1kg_pca
software/plink/plink --bfile 04_sampleqc/cohort_combine_snpQC_maf0.05 --pca --out 04_sampleqc/cohort_combine_snpQC_maf0.05
