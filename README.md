# MAGIC PIPELINE

## Prepare Software
Create a directory to store the downloaded software packages
```shell
mkdir packages
cd packages
```
Install bcftools
```shell
wget https://sourceforge.net/projects/samtools/files/samtools/1.9/bcftools-1.9.tar.bz2
tar -xvf bcftools-1.9.tar.bz2
cd bcftools-1.9/
./configure --prefix=/share/pub/likai/Star_protocols/software/bcftools-1.9
make && make install
```
Install samtools
```shell
wget https://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2
tar -xvf samtools-1.9.tar.bz2
cd samtools-1.9/
./configure --prefix=/share/pub/likai/Star_protocols/software/samtools-1.9
make && make install
```

## I. Data Pre-processing
Mapping clean fastq files to the reference genome and mark duplicates in the bam file.
```shell
cd /share/pub/likai/Star_protocols/
bash scripts/01/01_Map.sh
```
Realign bam files with GATK BaseRecalibrator and ApplyBQSR.
```shell
cd /share2/pub/lik/lik/Star_protocols
bash scripts/01/01_2_Realign.sh
```
Run GATK HaplotypeCaller on the vcf file.
```shell
cd /share2/pub/lik/lik/Star_protocols
bash scripts/01/01_3_HaplotypeCaller.sh
```
Run quality control on the fast file. If you are using another probe, please replace ref/Design_V2.hg19.bed with your own probe.
```shell
cd /share2/pub/lik/lik/Star_protocols
bash scripts/01/01_4_QC.sh
```
## II. Variant Discovery
Generate a file with paths to all gvcf files, one gvcf per line.
```shell
cd /share2/pub/lik/lik/Star_protocols
bash scripts/02/02_1_GvcfList.sh
```
Combine gvcf files of all samples and split the combined gvcf by chromosomes since this step is time-consuming.
```shell
cd /share/pub/likai/Star_protocols/scripts/02
for file in 02_2_CombineGvcf_chr*.sh
do
  echo $file
  bash $file
done 
```
Perform joint genotyping on the gvcf files.
```shell
cd /share/pub/likai/Star_protocols/scripts/02
for file in 02_3_GenotypeGvcf_chr*.sh
do
  echo $file
  bash $file
done 
```
## III. Variant QC
Filter sites by vcftools with the parameter --max-missing 0.9.
```shell
cd /share/pub/likai/Star_protocols/scripts/03
for file in 03_1_Filter_chr*.sh
do
  echo $file
  bash $file
done 
```
Combine the filtered vcf files into one vcf file.
```shell
cd /share/pub/likai/Star_protocols/
bash scripts/03/03_2_CombineVcf.sh
```
Run variant recalibrator on the combined vcf file.
```shell
cd /share/pub/likai/Star_protocols/
bash 03_3_VariantRecalibrator.sh
```
Run Variant Quality Score Recalibration.
```shell
cd /share/pub/likai/Star_protocols/
bash 03_4_VQSR.sh
```
Run variant recalibrator on the indel regions.
```shell
cd /share/pub/likai/Star_protocols
bash 03_5_VariantRecalibratorIndel.sh
```
Run Variant Quality Score Recalibration on the indel regions.
```shell
cd /share/pub/likai/Star_protocols
bash scripts/03/03_6_VQSRIndel.sh
```
Run select variants.
```shell
cd /share/pub/likai/Star_protocols
bash scripts/03/03_7_SelectVariants.sh
```
Split vcf file.
```shell
cd /share/pub/likai/Star_protocols/scripts/03
for file in 03_8_SplitVcf_chr*.sh
do
  echo $file
  bash $file
done 
```
Allele Balance for vcf file. Path to Rscript MUST be provided in the script. And you must
install package stringr, magrittr, readr and getopt in R.
```shell
cd /share/pub/likai/Star_protocols/
for file in scripts/03/03_9_AlleleBalance_chr*.sh
do
  echo $file
  bash $file
done 
```
Variant filtration
```shell
cd /share/pub/likai/Star_protocols/
for file in scripts/03/03_10_VariantFiltration_chr*.sh
do
  echo $file
  bash $file
done 
```
Hard Filtering
```shell
cd /share/pub/likai/Star_protocols/
for file in scripts/03/03_11_HardFiltering_chr*.sh
do
  echo $file
  bash $file
done 
```
## IV. Sample QC
Merge the hard filtered vcf into one
```shell
cd /share/pub/likai/Star_protocols
bash scripts/04/04_1_CombineChr.sh
```
Missing Indv
```shell
cd /share/pub/likai/Star_protocols
bash scripts/04/04_2_MissingIndv.sh
```
Depth Indv
```shell
cd /share/pub/likai/Star_protocols
bash scripts/04/04_3_DepthIndv.sh
```
GQ Indv
```shell
cd /share/pub/likai/Star_protocols
bash scripts/04/04_4_GQIndv.sh
```
Bcftools stats
```shell
bash scripts/04/04_5_VcfStats.sh
```
Run plink
```shell
bash scripts/04/04_6_Vcf2Plink.sh
```
Check Sex
```shell
bash scripts/04/04_7_CheckSex.sh
```
IBD
```shell
bash scripts/04/04_8_IBD.sh
```
PCA
```shell
bash scripts/04/04_9_PCA.sh
```
Indel Ratio
```shell
bash scripts/04/04_10_IndelRatio.sh
```
Filter sample with keep samples
```shell
bash scripts/04/04_11_FilterSample.sh
```
## V. Harmonization
Convert vcf to plink format file
```shell
bash scripts/05/05_1_Vcf2Plink.sh
```
## VI. Variant Annotation
Split vcf
```shell
for file in scripts/06/06_1_SplitVcf_Chr*.sh
do
  echo $file
  bash $file
done
```
Vcf fix, you need to install [vcffixup](https://github.com/vcflib/vcflib#INSTALL)
```shell
for file in scripts/06/06_2_FixVcf_Chr*.sh
do
  echo $file
  bash $file
done 
```
VEP annotation, need to install [VEP](https://asia.ensembl.org/info/docs/tools/vep/index.html)
```shell
for file in scripts/06/06_3_VEP_Chr*.sh
do
  echo $file
  bash $file
done 
```
Merge VEP annotation files
```shell
bash 06_4_VEPMerge.sh
```
## VII. ExWAS 
Convert VEP annotated vcf file to plink format file
```shell
bash scripts/07/07_1_Vcf2Plink.sh
```
PCA
```shell
bash scripts/07/07_2_PCA.sh
```
Emmax
```shell
bash scripts/07/07_3_Emmax.sh
```



