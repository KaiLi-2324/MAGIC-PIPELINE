# MAGIC PIPELINE

## Resources download

**Timing: 2.0 h**

In this section, the required non-software resources (reference genome and annotation files) are retrieved. At the end of the section, these resources should be stored in the proper locations for further steps.

1.Set the directory for storing reference genome and genomic annotations intended for subsequent general utilization.
```shell
>PROJECT_PATH=/user/projects
>RESOURCES_PATH=$PROJECT_PATH/ref 
>mkdir -p $ RESOURCES_PATH
>SOFTWARE_PATH=$PROJECT_PATH/software
>mkdir -p $SOFTWARE_PATH
>git clone https://github.com/sulab-wmu/MAGIC-PIPELINE.git
>mv MAGIC-PIPELINE/scripts $SOFTWARE_PATH/
```
2.Download the hg19 version of the human reference genome.
```shell
>cd $RESOURCES_PATH
>wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
>../software/seqkit grep -i -r -p '^chr[\dXY]+$' hg19.fa.gz -o ucsc.no_hap.hg19.fa.gz
>gzip -d ucsc.no_hap.hg19.fa.gz
```
3.Download the GATK resource bundle for Variant Quality Score Recalibration.
- Known variants and rs(dbSNP) accessions: dbSNP15p1
- HapMap genotypes and sites VCFs
- OMNI 2.5 genotypes for 1000 Genomes samples, as well as sites, VCF
- The set from 1000G phase 1 of known-site information for local realignment
- The current best set of known indels to be used for local realignment
```shell
>cd $RESOURCES_PATH
>wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
>gunzip 00-All.vcf.gz
>mv 00-All.vcf ncbi.hg19_dbsnp151.vcf
>wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/hapmap_3.3.hg19.sites.vcf.gz
>gunzip hapmap_3.3.hg19.sites.vcf.gz
>wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.gz
>gunzip 1000G_omni2.5.hg19.sites.vcf.gz
>wget  -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
>gunzip 1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
>wget  -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
>gunzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
```
4.Download variant frequencies across population files which are required in rare variant association study according to their frequency in major population cohorts: gnomAD 2.1.1 and index.
```shell
>cd $RESOURCES_PATH
>wget https://storage.googleapis.com/gcp-public-data– gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi
>wget https://storage.googleapis.com/gcp-public-data– gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz
```
5.Download 1000 genome database files for population structure analysis.
```shell
>cd $RESOURCES_PATH
>prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr";
suffix=".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
for chr in {1..22}; do wget "${prefix}""${chr}""${suffix}" "${prefix}""${chr}""${suffix}".tbi ; done
>wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
```
6.Download the resources for variant annotation in VEP.
```shell
>cd $RESOURCES_PATH
>wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz
>wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/InDels.tsv.gz
>wget https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz
>wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/GERP_scores.final.sorted.txt.gz
```
7.Acquire RNA-seq data in the relevant tissue for your trait of interest.
Identify either an existing RNA-seq data or a co-expression network from a relevant cell type or tissue for your trait of interest. There may be a publicly available source of RNA-seq data or a previously published co-expression network that is suitable for your applications. In this protocol, the gene-level reads per kilobase million mapped reads (RPKM) values of the Genotype-Tissue Expression (GTEx) RNA-seq data(https://www.gtexportal.org/home/) were used for peripheral retina samples.
```shell
>cd $RESOURCES_PATH
>wget https://storage.googleapis.com/gtex_external_datasets/eyegex_data/annotations/EyeGEx_meta_combined_inferior_retina_summary_deidentified_geo_ids.csv
>wget https://storage.googleapis.com/gtex_external_datasets/eyegex_data/rna_seq_data/EyeGEx_retina_combined_genelevel_expectedcounts_byrid_nooutlier.tpm.matrix.gct
```
## Download and install software and tools

**Timing: 1.5 h**

The goal of this section is to download and install software and tools for executing our pipeline. In cases where users encounter restriction or lack the necessary permissions to install software/tools on the machines, it is advisable to reach out to the machine administrator for assistance. The following command line operations can be executed as provided in most Linux distributions. Of note, at the end of each code snippet, a final command is included to export to tool’s command to the filesystem environment, ensuring its availability in subsequent steps. 
8.Installing conda and modifying channels in conda config
```shell
>wget https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh
>bash Anaconda3-2023.09-0-Linux-x86_64.sh
>conda config --add channels conda-forge    
>conda config --add channels defaults   
>conda config --add channels r    
>conda config --add channels bioconda
```
9.Installing Python prerequisites. The scripts were primarily tested and are recommended to run using Python 3.8. Run command:
```shell
>pip=`which pip3`
>$pip install numpy==1.19.4
>$pip install tqdm==4.42.1
>$pip install scipy==1.3.3
>$pip install statsmodels==0.12.1
>$pip install rpy2==3.3.6
```
10.Installing R packages SKAT and ACAT for burden analysis, stringr, magrittr, readr, getopt and WGCNA
```shell
>R_PATH=`which R`
>$R_PATH 
>>install.packages("SKAT")
>>install.packages("devtools")
>>devtools::install_github("yaowuliu/ACAT")
>>install.packages("stringr")
>>install.packages("magrittr")
>>install.packages("readr")
>>install.packages("getopt")
>>>install.packages("BiocManager") 
>>BiocManager::install("WGCNA")
```
11.Download and install seqkit.
```shell
>cd $SOFTWARE_PATH
>wget wget https://github.com/shenwei356/seqkit/releases/download/v2.5.1/seqkit_linux_amd64.tar.gz
>tar -xvf seqkit_linux_amd64.tar.gz
```
12.Download and install bwa.
```shell
>cd $SOFTWARE_PATH
>mkdir packages
>cd packages
>wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2 
>tar -xvf bwa-0.7.12.tar.bz2 
>export BWA_PATH=$SOFTWARE_PATH/bwa-0.7.12
>cd $BWA_PATH  
>make 
```
13.Download and install GATK. 
```shell
>cd $SOFTWARE_PATH
>wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip
>unzip gatk-4.3.0.0.zip
```
14.Download and install vcflib.
```shell
>conda create -n vcflib
>conda install -c bioconda -n vcflib vcflib
```
15.Download and install VEP.
```shell
>cd $SOFTWARE_PATH
>git clone https://github.com/Ensembl/ensembl-vep.git
>cd ensembl-vep
>perl INSTALL.pl
```
16.Download and install plink.
```shell
>cd $SOFTWARE_PATH
>wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip
>unzip plink_linux_x86_64_20230116.zip
```
17.Download and install EMMAX.
```shell
>cd $SOFTWARE_PATH
>wget http://csg.sph.umich.edu//kang/emmax/download/emmax-intel-binary-20120210.tar.gz
>tar xvf emmax-intel-binary-20120210.tar.gz
```
18.Download and install GCTA.
```shell
>cd $SOFTWARE_PATH
>wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip
>unzip gcta-1.94.1-linux-kernel-3-x86_64.zip
>cp -r gcta-1.94.1-linux-kernel-3-x86_64 $SOFTWARE_PATH/gcta64
```
19.Download and install vcftools.
```shell
>cd $SOFTWARE_PATH
>wget https://sourceforge.net/projects/vcftools/files/latest/download/ vcftools_0.1.13.tar.gz
>tar -xzf vcftools_0.1.13.tar.gz
>cd vcftools_0.1.13
>./configure --prefix=$SOFTWARE_PATH/vcftools_0.1.13
>make && make install
```
20.Download and install samtools.
```shell
>cd $SOFTWARE_PATH
>wget https://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2
>tar -xvf samtools-1.9.tar.bz2
>cd samtools-1.9
>./configure --prefix=$SOFTWARE_PATH/samtools-1.9
>make && make install
```
21.Download and install bcftools. 
```shell
>cd $SOFTWARE_PATH
>wget https://sourceforge.net/projects/samtools/files/samtools/1.9/bcftools-1.9.tar.bz2
>tar -xvf bcftools-1.9.tar.bz2
>cd bcftools-1.9/
>./configure --prefix=$SOFTWARE_PATH/bcftools-1.9
>make && make install
```
22.Download and install sambamba
```shell
>cd $SOFTWARE_PATH
>wget https://github.com/biod/sambamba/releases/download/v0.6.6/sambamba_v0.6.6_linux.tar.bz2
>tar -xvf sambamba_v0.6.6_linux.tar.bz2
```
23.Download and install tabix and bgzip
```shell
>cd $SOFTWARE_PATH && mkdir htslib
>wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib1.10.2.tar.bz2
>tar -xvf htslib-1.10.2.tar.bz2
>cd htslib-1.10.2 && ./configure --prefix=$SOFTWARE_PATH/htslib
>make && make install
```
## Download the test dataset

**Timing: 1.5 h**

In this section, the raw data for the demonstration of the protocol are retrieved. At the end of the process, the data should be in the proper place for the continuation of the protocol. 

24.Set the directory where the raw data will be placed and download the raw WES data.
```shell
>cd $PROJECT_PATH
>mkdir rawdata
>cd rawdata
>wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099967/SRR099967_1.fastq.gz 
>wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099967/SRR099967_2.fastq.gz 
> wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099969/SRR099969_1.fastq.gz
>wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099969/SRR099969_2.fastq.gz
>wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099957/SRR099957_1.fastq.gz >wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099957/SRR099957_2.fastq.gz
```
25.Download test file for exome-wide association and burden test.
```shell
>cd $PROJECT_PATH
>mkdir 06_ExWAS/
>curl -o 06_ExWAS/ExWAS_test.fam \ https://zenodo.org/record/8189201/files/ExWAS_test.fam
>curl -o 06_ExWAS/ExWAS_test.bim \ https://zenodo.org/record/8189201/files/ExWAS_test.bim
>curl -o 06_ExWAS/ExWAS_test.bed \ https://zenodo.org/record/8189201/files/ExWAS_test.bed
>mkdir 07_collapsing/
>curl -o 07_collapsing/test.vcf \ https://zenodo.org/record/8189201/files/test.vcf
>curl -o 07_collapsing/test_class.txt \ https://zenodo.org/record/8189201/files/test_class.txt
>curl -o ref/Design_V2.hg19.bed \ https://zenodo.org/record/8189201/files/Design_V2.hg19.bed
```
26.Download RNAseq test file for WGCNA.
```shell
>mkdir 08_systems_genetics
>curl -o \
08_systems_genetics/datExpr1_combine_test.txt \
https://zenodo.org/record/8189201/files/datExpr1_combine_test.txt
>curl -o 08_systems_genetics/top50.txt \
https://zenodo.org/record/8189201/files/top50.txt
```
## Step-by-step method details
### Step 1: Data Pre-processing
**Timing: 4.0 h**

This section outlines the procedure of aligning the FASTQ pairs to the reference genome and collecting alignment statistics to ensure quality control. The output of this part comprises BAM files that are appropriate for the subsequent variant calling. At the end of the step, both BAM files and a report containing read alignment statistics are generated.
1.Index the reference genome.
```shell
>cd $RESOURCES_PATH
>$BWA_PATH/bwa index ucsc.no_hap.hg19.fa
>$SAMTOOLS_PATH/samtools faidx ucsc.no_hap.hg19.fa
>$SAMTOOLS_PATH/samtools dict ucsc.no_hap.hg19.fa > ucsc.no_hap.hg19.fa
```
2.Create a file with only one sample per line, then align each sample to the reference genome and mark duplicates in the BAM file. 
```shell
>cd $PROJECT_PATH
>mkdir 01_mapping/
>cat sample.list
SRR099967
SRR099969
SRR099957
>cat sample.list | while read line
do
	software/bwa-0.7.12/bwa mem -t 32 -M -R "@RG\tID:${line}\tPL:illumina\tSM:${line}\tCN:WY" ref/ucsc.no_hap.hg19.fa rawdata/${line}_1.fastq.gz rawdata/${line}_2.fastq.gz | software/samtools-1.9/bin/samtools view -uS -F 4 - -o 01_mapping/${line}.bam && software/sambamba-0.6.6/sambamba sort -t 4 -u --tmpdir=temp 01_mapping/${line}.bam && software/sambamba-0.6.6/sambamba markdup --sort-buffer-size 262144 --hash-table-size 60000 -t 4 -l 0 --overflow-list-size 400000 --tmpdir=temp 01_mapping/${line}.sorted.bam 01_mapping/${line}.sort.rmdup.bam
done
```
3.Base Quality Score Recalibration and application on BAM files.
```shell
>cd $PROJECT_PATH
>export _JAVA_OPTIONS=-Djava.io.tmpdir=temp/
>cat sample.list | while read line
do
	software/gatk-4.3.0.0/gatk --java-options \
"-Dsamjdk.use_async_io_read_samtools=true \
-Dsamjdk.use_async_io_write_samtools=true \
-Dsamjdk.use_async_io_write_tribble=false -Xmx8g" \
BaseRecalibrator -R ref/ucsc.no_hap.hg19.fa \
-I 01_mapping/${line}.sort.rmdup.bam --verbosity WARNING \
--known-sites ref/ncbi.hg19_dbsnp151.vcf -O 01_realign/${line}.table
software/gatk-4.3.0.0/gatk --java-options \
"-Dsamjdk.use_async_io_read_samtools=true \
-Dsamjdk.use_async_io_write_samtools=true \
-Dsamjdk.use_async_io_write_tribble=false -Xmx8g" \
ApplyBQSR -R ref/ucsc.no_hap.hg19.fa -I \
01_mapping/${line}.sort.rmdup.bam --verbosity WARNING -bqsr \
01_realign/${line}.table -O 01_realign/${line}.bam 
software/sambamba-0.6.6/sambamba index 01_realign/${line}.bam
done
```
4.Collection of alignment statistics. 
If you are using another WES probe, please replace ref/Design_V2.hg19.bed with your own probe. As previous WES analysis 18, we collect several statistics that pertain to the alignment process's quality control in this step. At the conclusion of the process, a text file containing the statistics should be generated.  
- Total sequenced reads.
- Aligned reads.
- Uniquely aligned reads (q>20).
- Reads overlapping targets.
- Total sequenced bases.
- Aligned bases.
- Uniquely aligned bases.
- Bases overlapping targets.
```shell
>cd $PROJECT_PATH
>RSCRIPT=`which Rscript`
>PERL=`which perl`
>cat sample.list | while read line
do
	$PERL software/QC_exome_mem.pl -R $RSCRIPT -i 01_mapping/${line}.sort.rmdup.bam -r ref/Design_V2.hg19.bed \
-c rawdata/${line}.R1.clean.fastq.gz -o 01_vcf/${line} -s software/samtools-1.9/bin/samtools -b software/bedtools2/bin/bedtools -p software/picard.jar -t PE -plot \
-ref ref/ucsc.no_hap.hg19.fa -d ref/ucsc.no_hap.hg19.dict -C ref/Design_V2.hg19.bed
done 
```
### Step 2: Variant Discovery	

**Timing: 2.0 h**

This section provides an overview of the variant calling procedure using GATK HaplotypeCaller. The BAM files generated during preprocessing are fed into GATK’s HaplotypeCaller to perform per-sample variant calling, outputting a genomic VCF (gVCF) for a single sample. Joint-calling and variant quality score recalibration (VQSR) generates a final, multi-sample VCFs files, serving as the starting point for downstream analysis.  
5.For each sample use GATK HaplotypeCaller to create a gVCF callset file.
```shell
>cd $PROJECT_PATH
>mkdir 02_gvcf && mkdir 02_vcf
>export _JAVA_OPTIONS=-Djava.io.tmpdir=temp/
>cat sample.list | while read line
do
	software/gatk-4.3.0.0/gatk --java-options "-Xmx8G" HaplotypeCaller -R ref/ucsc.no_hap.hg19.fa -I 01_realign/${line}.bam -L ref/Design_V2.hg19.bed -ERC BP_RESOLUTION -O 01_vcf/${line}_HaplotypeCaller.g.vcf && software/bgzip 01_vcf/${line}_HaplotypeCaller.g.vcf && software/tabix -p vcf 01_vcf/${line}_HaplotypeCaller.g.vcf.gz
done 
```
6.Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file.
```shell
>cd $PROJECT_PATH
>ls 01_vcf/*.gz > 02_gvcf/cohort.vcf.list
>software/gatk-4.3.0.0/gatk CombineGVCFs --java-options "-Xmx150g -Djava.io.tmpdir=temp/" -R ref/ucsc.no_hap.hg19.fa --variant 02_gvcf/cohort.vcf.list -O 02_gvcf/cohort.g.vcf
```
7.Perform joint-genotyping on the combined gVCF files.
```shell
>cd $PROJECT_PATH
>software/gatk-4.3.0.0/gatk GenotypeGVCFs --java-options "-Xmx150g -Djava.io.tmpdir=temp/" --variant 02_gvcf/cohort.g.vcf -new-qual true -R ref/ucsc.no_hap.hg19.fa -O 02_vcf/cohort.vcf.gz
```
8.Filter variants by hard-filtering before variant quality score recalibration (VQSR).
```shell
>cd $PROJECT_PATH
>software/vcftools/bin/vcftools --gzvcf 02_vcf/cohort.vcf.gz --max-missing 0.9 --recode --recode-INFO-all --out 03_filter/cohort.mis0.9
```
9.Calculate VQSLOD tranches for SNPs using VariantRecalibrator.
```shell
>cd $PROJECT_PATH
>software/gatk-4.3.0.0/gatk --java-options "-Xmx248g" VariantRecalibrator -R ref/ucsc.no_hap.hg19.fa -V 02_vcf/cohort.mis0.9.recode.vcf \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ref/hapmap_3.3.hg19.sites.vcf \
-resource:omni,known=false,training=true,truth=true,prior=12.0 ref/1000G_omni2.5.hg19.sites.vcf \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 ref/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ref/ncbi.hg19_dbsnp151.vcf \
-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP --max-gaussians 4 -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
-O 02_vcf/cohort_VQSR.snp.recal \
--tranches-file 02_vcf/cohort_VQSR.snp.tranches \
--rscript-file 02_vcf/cohort_VQSR.snp.plots.R
```
10.Filter SNPs on VQSLOD using ApplyVQSR.
```shell
>cd $PROJECT_PATH
>software/gatk-4.3.0.0/gatk --java-options "-Xmx248g" ApplyVQSR \
-R ref/ucsc.no_hap.hg19.fa -V 02_vcf/cohort.mis0.9.recode.vcf \
--ts-filter-level 99.9 
--tranches-file 02_vcf/cohort_VQSR.snp.tranches \
--recal-file 02_vcf/cohort_VQSR.snp.recal \
-mode SNP -O 02_vcf/cohort_VQSR.snp.recal.vcf
```
11.Calculate VQSLOD tranches for indels using VariantRecalibrator.
```shell
>cd $PROJECT_PATH
>software/gatk-4.3.0.0/gatk --java-options "-Xmx248g" VariantRecalibrator -R ref/ucsc.no_hap.hg19.fa -V 02_vcf/cohort_VQSR.snp.recal.vcf \
-resource:mills,known=true,training=true,truth=true,prior=12.0 ref/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ref/ncbi.hg19_dbsnp151.vcf \
-an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode INDEL \
--max-gaussians 4 -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
-O 02_vcf/cohort_VQSR.indel.recal \
--tranches-file 02_vcf/cohort_VQSR.indel.tranches \
--rscript-file 02_vcf/cohort_VQSR.indel.plots.R
```
12.Filter indels on VQSLOD using ApplyVQSR.
```shell
>cd $PROJECT_PATH
>software/gatk-4.3.0.0/gatk --java-options "-Xmx248g" ApplyVQSR \
-R ref/ucsc.no_hap.hg19.fa -V 02_vcf/cohort_VQSR.snp.recal.vcf \
--ts-filter-level 99.9 \
--tranches-file 02_vcf/cohort_VQSR.indel.tranches \
--recal-file 02_vcf/cohort_VQSR.indel.recal -mode INDEL \
-O 02_vcf/cohort_VQSR.snp.indel.recal.vcf
```
13.Subset to SNPs and indels callset with SelectVariants.
```shell
>cd $PROJECT_PATH
>software/gatk-4.3.0.0/gatk --java-options "-Xmx248g" SelectVariants -R ref/ucsc.no_hap.hg19.fa -V 02_vcf/cohort_VQSR.snp.indel.recal.vcf \
--exclude-filtered --preserve-alleles -O 02_vcf/cohort_after_VQSR.vcf
```
### Step 3: Quality control (QC) for variants
**Timing: 0.5 h**

In this section, we performed QC on genetic variants discovered through sequencing. Variants were prefiltered so that only those passed the GATK VQSR metric and those located outside of low-complexity regions were included into further analysis. Genotypes with a genotype depth (DP) < 10 and genotype quality (GQ) < 20 and heterozygous genotype calls with an allele balance >0.8 or <0.2 were determined as missing data. We excluded variants with a call rate<0.9, a case-control call rate difference >0.005, or a Hardy-Weinberg equilibrium (HWE) test P-value <10-6 on the basis of the combined case-control cohort.
14.Perform Allele balance on a VCF file from GATK. Path to Rscript must be provided in the script. And you must prioritize install package stringr, magrittr, readr and getopt in R.
```shell
>cd $PROJECT_PATH
>mkdir 03_variant_qc
>RSCRIPT=`which Rscript`
>$Rscript software/scripts/AlleleBalanceBySample.R -v 02_vcf/cohort_after_VQSR.vcf -o 03_variant_qc/cohort_after_VQSR.AB.vcf -i 0.2 -I 0.8
```
15.Hard filter a cohort callset with VariantFiltration.
```shell
>cd $PROJECT_PATH
>software/gatk-4.3.0.0/gatk --java-options "-Xmx100g" VariantFiltration -R ref/ucsc.no_hap.hg19.fa -G-filter "GQ < 20.0" -G-filter-name lowGQ \
-G-filter "DP < 10.0" -G-filter-name lowDP --set-filtered-genotype-to-no-call -V 03_variant_qc/cohort_after_VQSR.AB.vcf -O 03_variant_qc/cohort_after_VQSR.AB.VF.vcf
```
16.Apply final hard filters to a cohort callset.
```shell
>cd $PROJECT_PATH
>software/vcftools/bin/vcftools --vcf 03_variant_qc/cohort_after_VQSR.AB.VF.vcf --max-alleles 2 --hwe 0.000001 --max-missing 0.9 \
--recode --recode-INFO-all --out 03_variant_qc/cohort_snpQC
```
###  Step 4: Quality control for samples

**Timing: 0.5 h**

In this section, we performed sample QC. Samples were excluded if they showed a low average call rate (<0.9), low mean sequencing depth (<10), or low mean genotype quality (<65). Outliers (>4 SD from the mean) of the transition/transversion ratio, heterozygous/homozygous ratio, or insertion/deletion ratio within each cohort were further discarded. We also removed samples that were closely related to one another, had ambiguous sex status or had population outliers per principal component analysis (PCA).
17.Sample-based missing data reports.
```shell
>cd $PROJECT_PATH
>mkdir 04_sample_qc
>software/vcftools/bin/vcftools --vcf 03_variant_qc/cohort_snpQC.recode.vcf --missing-indv --out 04_sample_qc/cohort_snpQC
```
18.Calculate mean depth per sample.
```shell
>cd $PROJECT_PATH
>software/vcftools/bin/vcftools --vcf 03_variant_qc/cohort_snpQC.recode.vcf --depth --out 04_sample_qc/cohort_snpQC
```
19.Calculate mean genotype-quality per sample.
```shell
>cd $PROJECT_PATH
>software/vcftools/bin/vcftools --vcf 03_variant_qc/cohort_snpQC.recode.vcf --extract-FORMAT-info GQ --out 04_sample_qc/cohort_snpQC
cat 04_sample_qc/cohort_snpQC.GQ.FORMAT | sed '1d' | awk '{for(i=1;i<=NF;i++) total[i]+=$i;} END {for(i=1;i<=NF;i++) printf "%f ",total[i]/NR ;}' | sed 's/ /\n/g' > 04_sample_qc/cohort_snpQC.GQ
```
20.Calculate transition/transversion ratio, heterozygous/homozygous ratio, and insertion/deletion ratio from bcftools stats.
```shell
>cd $PROJECT_PATH
>software/bcftools-1.9/bin/bcftools stats -s - 03_variant_qc/cohort_snpQC.recode.vcf > 04_sample_qc/cohort_snpQC_stats
>cat 04_sample_qc/cohort_snpQC_stats | grep "PSC" | grep -v "#" | awk '{print $3"\t"$7"\t"$8}' > 04_sample_qc/cohort_snpQC_stats.tstv
>cat 04_sample_qc/cohort_snpQC_stats | grep "PSC" | grep -v "#" | awk '{print $3"\t"$5"\t"$6}' > 04_sample_qc/cohort_snpQC_stats.homhet
>cat 04_sample_qc/cohort_snpQC_stats | grep "PSC" | grep -v "#" | awk '{print $3"\t"$5+$6"\t"$9}'> 04_sample_qc/cohort_snpQC_stats.indelsnp
>cat 04_sample_qc/cohort_snpQC_stats | grep "PSC" | grep -v "#" | awk '{print $3"\t"$11}' > 04_sample_qc/cohort_snpQC_stats.singleton
>software/vcftools/bin/vcftools --vcf 03_variant_qc/cohort_snpQC.recode.vcf --keep-only-indels --recode --recode-INFO-all --out 04_sample_qc/cohort_snpQC.indel
>software/vcftools/bin/vcftools --vcf 04_sample_qc/cohort_snpQC.indel.recode.vcf --extract-FORMAT-info GT --out 04_sample_qc/cohort_snpQC.indel
>python3 software/scripts/calc_indel_ratio.py --path_indel_vcf 04_sample_qc/cohort_snpQC.indel.recode.vcf --path_out 04_sample_qc/cohort_indel_ratio.txt --path_gt 04_sample_qc/cohort_snpQC.indel.GT.FORMAT
```
21.Sex discrepancy. 
```shell
>cd $PROJECT_PATH
>software/plink/plink --vcf 03_variant_qc/cohort_snpQC.recode.vcf --recode --make-bed --out 04_sample_qc/cohort_snpQC
>cat 04_sample_qc/cohort_snpQC.bim | awk '{print $1"\t"$1":"$4":"$6":"$5"\t"$3"\t"$4"\t"$5"\t"$6}' > 04_sample_qc/cohort_snpQC.bim1
>mv 04_sample_qc/cohort_snpQC.bim1 04_sample_qc/cohort_snpQC.bim
>software/plink/plink --bfile 04_sample_qc/cohort_snpQC --exclude ref/inversion.txt --range --indep-pairwise 50 5 0.2 --out 04_sample_qc/indepSNP
>software/plink/plink --bfile 04_sample_qc/cohort_snpQC --extract ref/indepSNP.prune.in --check-sex --set-hh-missing --out 04_sample_qc/cohort_snpQC
```
22.Relatedness check by identity-by-descent
```shell
>cd $PROJECT_PATH
>software/plink/plink --bfile 04_sample_qc/cohort_snpQC --extract ref/indepSNP.prune.in --genome --out 04_sample_qc/cohort_snpQC
```
23.Estimating population structure and sample ancestry. Linkage Disequilibrium (LD) pruning on the variants were performed before running PCA.
```shell
>cd $PROJECT_PATH
>software/plink/plink --bfile ref/Merge --geno 0.1 --mind 0.1 --maf 0.05 --allow-no-sex --make-bed --out 04_sample_qc/1kG_PCA
>software/plink/plink --bfile 04_sample_qc/cohort_snpQC --geno 0.1 --mind 0.1 --maf 0.05 --make-bed --out 04_sample_qc/cohort_snpQC_maf0.05
>awk '{print$2}' 04_sample_qc/cohort_snpQC_maf0.05.bim > 04_sample_qc/cohort_snpQC_maf0.05.txt
>software/plink/plink --bfile 04_sample_qc/1kG_PCA --extract 04_sample_qc/cohort_snpQC_maf0.05.txt --make-bed --out 04_sample_qc/1kG_PCA1
>awk '{print $2}' 04_sample_qc/1kG_PCA1.bim > 04_sample_qc/1kG_PCA1.txt
>software/plink/plink --bfile 04_sample_qc/cohort_snpQC_maf0.05 --extract 04_sample_qc/1kG_PCA1.txt --recode --make-bed --out 04_sample_qc/cohort_snpQC_maf0.05_1kg
>awk '{print $2,$4}' 04_sample_qc/cohort_snpQC_maf0.05_1kg.map > 04_sample_qc/buildhapmap.txt
>software/plink/plink --bfile 04_sample_qc/1kG_PCA1 --update-map 04_sample_qc/buildhapmap.txt --make-bed --out 04_sample_qc/1kG_PCA2
>software/plink/plink --bfile 04_sample_qc/cohort_snpQC_maf0.05_1kg --bmerge 04_sample_qc/1kG_PCA2.bed 04_sample_qc/1kG_PCA2.bim 04_sample_qc/1kG_PCA2.fam --allow-no-sex --make-bed --out 04_sample_qc/cohort_snpQC_maf0.05_1kg_pca
>software/plink/plink --bfile 04_sample_qc/cohort_snpQC_maf0.05_1kg_pca --exclude ref/inversion.txt --range --indep-pairwise 50 5 0.2 --out 04_sample_qc/indepSNP
>software/plink/plink --bfile 04_sample_qc/cohort_snpQC_maf0.05_1kg_pca --extract ref/indepSNP.prune.in --pca --out 04_sample_qc/cohort_snpQC_maf0.05_1kg_pca
>software/plink/plink --bfile 04_sample_qc/cohort_snpQC_maf0.05 --exclude ref/inversion.txt --range --indep-pairwise 50 5 0.2 --out 04_sample_qc/indepSNP
>software/plink/plink --bfile 04_sample_qc/cohort_snpQC_maf0.05 --extract ref/indepSNP.prune.in --pca --out 04_sample_qc/cohort_snpQC_maf0.05
```
24.Retain samples that pass our quality control filters.
```shell
>cd $PROJECT_PATH
>software/vcftools/bin/vcftools --vcf 03_variant_qc/cohort_snpQC.recode.vcf --keep 04_sample_qc/keep_samples.txt \
--non-ref-ac-any 1 --hwe 0.000001 --max-missing 0.9 --recode --recode-INFO-all --out 04_sample_qc/cohort_snpQC_sampleQC
```
### Step 5: Variant Annotation

**Timing: 4.0 h**

In this section, the annotation of variants was performed with Ensembl’s VEP for human genome assembly GRCh37. We used the VEP, CADD, LOFTEE and SpliceAI plugins to generate additional bioinformatic predictions of variant deleteriousness.
25.Uses genotypes from the VCF file to correct AC (alternate allele count), AF (alternate allele frequency), NS (number of called), in the VCF file.
```shell
>cd $PROJECT_PATH
>mkdir 05_variantannot
>conda activate vcflib
>software/vcflib/vcffixup 04_sample_qc/cohort_snpQC_sampleQC.recode.vcf > 05_variantannot/cohort_snpQC_sampleQC_fix.recode.vcf
```
26.VEP annotation.
```shell
>cd $PROJECT_PATH
>software/ensembl-vep/vep -i 05_variantannot/cohort_snpQC_sampleQC_fix.recode.vcf --offline --cache --dir_cache ref/vep/vep_cache/ --plugin CADD,ref/vep/Plugins/whole_genome_SNVs.tsv.gz,ref/vep/Plugins/InDels.tsv.gz --plugin LoF,loftee_path:ref/vep/Plugins/loftee/loftee/,human_ancestor_fa:ref/vep/Plugins/loftee/human_ancestor.fa.gz,conservation_file:ref/vep/Plugins/loftee/phylocsf_gerp.sql,gerp_file:ref/vep/Plugins/loftee/GERP_scores.final.sorted.txt.gz --plugin SpliceAI,snv=ref/vep/vcf/spliceai_scores.raw.snv.hg19.vcf.gz,indel=ref/vep/vcf/spliceai_scores.raw.indel.hg19.vcf.gz --dir_plugins ref/vep/Plugins --custom ref/vep/vcf/gnomad.genomes.r2.1.1.sites.vcf.gz,gnomADg,vcf,exact,0,AC,AF,AF_afr,AF_amr,AF_asj,AF_eas,AF_fin,AF_nfe,AF_oth --fasta ref/ucsc.no_hap.hg19.fa --force_overwrite --sift b --polyphen b --regulatory --numbers --ccds --hgvs --symbol --xref_refseq --canonical --protein --biotype --af --af_1kg --af_esp --af_gnomad --max_af --pick --vcf --output_file 05_variantannot/cohort_snpQC_sampleQC_fix_vep.vcf
```
### Step 6: Exome-wide single-variant association analysis 

**Timing: 0.5 h**

In this section, we conducted two types of single-variant association analyses, including the EMMAX test and MLMA, on all samples. For a rapid test of this protocol, we provided test dataset from our large cohort.
27.Generate in-house principal component as a covariate after sample and variant QC. LD pruning on the variants were performed before running PCA.
```shell
>cd $PROJECT_PATH
>software/plink/plink --bfile 06_ExWAS/ExWAS_test --indep-pairwise 50 5 0.2 --out 06_ExWAS/ExWAS_test_indep
>software/plink/plink --bfile 06_ExWAS/ExWAS_test --extract 06_ExWAS/ExWAS_test_indep.prune.in --maf 0.05 --snps-only --make-bed --out 06_ExWAS/ExWAS_test_maf0.05
>software/plink/plink --bfile 06_ExWAS/ExWAS_test_maf0.05 --pca --autosome --out 06_ExWAS/ExWAS_test_PCA
```
28.Performed ExWAS using EMMAX.
```shell
>cd $PROJECT_PATH
>software/plink/plink --bfile 06_ExWAS/ExWAS_test --recode12 --output-missing-genotype 0 --transpose --out 06_ExWAS/ExWAS_test_emmax
>software/emmax-kin-intel64 -v -d 10 06_ExWAS/ExWAS_test_emmax
>software/emmax-intel64 -v -d 10 -t 06_ExWAS/ExWAS_test_emmax -p 06_ExWAS/ExWAS_test.pheno -k 06_ExWAS/ExWAS_test_emmax.aBN.kinf -o 06_ExWAS/ExWAS_test_emmax.out
>cat 06_ExWAS/ExWAS_test_emmax.out.ps | cut -f1 | sed 's/:/\t/g' | paste - 06_ExWAS/ExWAS_test_emmax.out.ps > 06_ExWAS/ExWAS_test_emmax.out.txt
```
29.Performed ExWAS using MLMA.
```shell
>cd $PROJECT_PATH
>software/gcta64 --bfile 06_ExWAS/ExWAS_test --autosome --make-grm --out 06_ExWAS/ExWAS_test --thread-num 5
>software/gcta64 --mlma --bfile 06_ExWAS/ExWAS_test --grm 06_ExWAS/ExWAS_test --pheno 06_ExWAS/ExWAS_test.pheno --out 06_ExWAS/ExWAS_test_MLMA --thread-num 5
>software/gcta64 --mlma-loco --bfile 06_ExWAS/ExWAS_test --grm 06_ExWAS/ExWAS_test --pheno 06_ExWAS/ExWAS_test.pheno --out 06_ExWAS/ExWAS_test_MLMALOCO --thread-num 5
```
### Step 7: Gene-based collapsing analysis

**Timing: 1.5 h**

To determine whether a single gene was enriched in or depleted of rare protein-coding variants in cases, we performed four gene-level association tests including Fisher’s exact test, burden, SKAT and SKAT-O.
30.Qualifying variants.
```shell
>cd $PROJECT_PATH
>software/vcftools/bin/vcftools --vcf 07_collapsing/test.vcf --mac 1 --max-maf 0.005 --recode --recode-INFO-all --out 07_collapsing/test_rare
>python3 software/scripts/vcf2mat.py -p 07_collapsing/test_rare.recode.vcf -o 07_collapsing/test_rare_lof.txt -e MAX_AF -F 0.005 -v lof -c 07_collapsing/test_rare_lof.out
```
31.Burden tests for rare variants using Fisher's method.
```shell
>cd $PROJECT_PATH
>python3 software/scripts/fisher_test_burden.py -p 07_collapsing/test_rare_lof.txt -c 07_collapsing/test_class.txt -o 07_collapsing/test_rare_lof_fisher.txt -l 40
```
32.Burden tests for rare variants using SKAT method.
```shell
>cd $PROJECT_PATH
>python3 software/scripts/skat_burden.py -p 07_collapsing/test_rare_lof.txt -c 07_collapsing/test_class.txt -o 07_collapsing/test_rare_lof_skat.txt
```
###  Step 8: Systems genetics

**Timing: 2.0 h**

In this section, we build a signed coexpression network for test dataset derived from peripheral retina by using weighted gene coexpression network analysis (WGCNA). Additionally, we calculated the odds ratios (ORs) and p-values of candidate risk genes enriched for module genes from WGCNA.
33.Weighted gene coexpression network analysis and module enrichment analysis in R.
```shell
>cd $PROJECT_PATH/08_systems_genetics
>RSCRIPT=`which Rscript`
>$RSCRIPT ../software/scripts/WGCNA.R ../rawdata/datExpr1_combine_test.txt
```



