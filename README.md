## Before you begin

Whole-exome sequencing (WES) is a research approach for identifying molecular defects in patients with suspected genetic disorders2-4. Large-scale WES analysis has effectively identified many novel susceptibility genes and pathways for various dieases5-7. This protocol was established to analyze the original WES data in large cohorts. For data pre-processing and variant discovery steps, the Genome Analysis Toolkit (GATK) is a recommended tool with workflows following community best practices8,9. Variant annotation and variant frequency addition from large population studies can be processed by referring to the instructions. The execution of this protocol requires additional resources beyond the basic software tools, which can be categorized into two types, i.e., reference file downloads (reference genome, genomic coordinates of the kit used for WES, and reference annotation databases) and additional software tools (for pre- and post-processing of the input and output data). Apart from data collection, the following steps can be performed once, as they involve online data and locally stored software resources. The users of this protocol were assumed to have already had a basic familiarity with the Unix/Linux command line as all commands are executed in the command line via a terminal.
Note: each line of executable code throughout this protocol is preceded by a greater-than-sign (>). 

## Resources download

Timing: 2.0 h

In this section, the required non-software resources (reference genome and annotation files) are retrieved. At the end of the section, these resources should be stored in the proper locations for further steps.

1. Set the directory for storing reference genome and genomic annotations intended for subsequent general utilization.
```
>PROJECT_PATH=/user/projects
>RESOURCES_PATH=$PROJECT_PATH/ref 
>mkdir -p $ RESOURCES_PATH
>SOFTWARE_PATH=$PROJECT_PATH/software
>mkdir -p $SOFTWARE_PATH
>git clone https://github.com/sulab-wmu/MAGIC-PIPELINE.git
>mv MAGIC-PIPELINE/scripts $SOFTWARE_PATH/
```
2. Download the hg19 version of the human reference genome.
```
>cd $RESOURCES_PATH
>wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
>../software/seqkit grep -i -r -p '^chr[\dXY]+$' hg19.fa.gz -o ucsc.no_hap.hg19.fa.gz
>gzip -d ucsc.no_hap.hg19.fa.gz
```
3. Download the GATK resource bundle for Variant Quality Score Recalibration.
a.Known variants and rs(dbSNP) accessions: dbSNP151
b.HapMap genotypes and sites VCFs
c.OMNI 2.5 genotypes for 1000 Genomes samples, as well as sites, VCF
d.The set from 1000G phase 1 of known-site information for local realignment
e.The current best set of known indels to be used for local realignment
```
>cd $RESOURCES_PATH
>wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
>gunzip 00-All.vcf.gz
>mv 00-All.vcf ncbi.hg19_dbsnp151.vcf
>wget ftp://ftp.broadinstitute.org/bundle/hg19/hapmap_3.3.hg19.sites.vcf.gz
>gunzip hapmap_3.3.hg19.sites.vcf.gz
>wget ftp://ftp.broadinstitute.org/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.gz
>gunzip 1000G_omni2.5.hg19.sites.vcf.gz
>wget >ftp://ftp.broadinstitute.org/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
>gunzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
>wget >ftp://ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
>gunzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
```
4. Download variant frequencies across population files which are required in rare variant association study according to their frequency in major population cohorts: gnomAD 2.1.1 and index.
```
>cd $RESOURCES_PATH
>wget https://storage.googleapis.com/gcp-public-data– gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi
>wget https://storage.googleapis.com/gcp-public-data– gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz
>wget https://storage.googleapis.com/gcp-public-data– gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi
>wget https://storage.googleapis.com/gcp-public-data– gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz
```
5. Download 1000 genome database files for population structure analysis.
```
>cd $RESOURCES_PATH
>prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr";
suffix=".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
for chr in {1..22}; do wget "${prefix}""${chr}""${suffix}" "${prefix}""${chr}""${suffix}".tbi ; done
>wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
```
6. Download the resources for variant annotation in VEP.
```
>cd $RESOURCES_PATH
>wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/whole_genome_SNVs.tsv.gz
>wget https://krishna.gs.washington.edu/download/CADD/v1.5/GRCh38/InDels.tsv.gz
>wget https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz
>wget https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/GERP_scores.final.sorted.txt.gz
```
Note: The annotations for all possible SNVs and indels in spliceAI were download from https://basespace.illumina.com/s/otSPW8hnhaZR.

7. Acquire RNA-seq data in the relevant tissue for your trait of interest.
Identify either an existing RNA-seq data or a co-expression network from a relevant cell type or tissue for your trait of interest. There may be a publicly available source of RNA-seq data or a previously published co-expression network that is suitable for your applications. In this protocol, the gene-level reads per kilobase million mapped reads (RPKM) values of the Genotype-Tissue Expression (GTEx) RNA-seq data(https://www.gtexportal.org/home/) were used for peripheral retina samples.
```
>cd $RESOURCES_PATH
>wget https://storage.googleapis.com/gtex_external_datasets/eyegex_data/annotations/EyeGEx_meta_combined_inferior_retina_summary_deidentified_geo_ids.csv
>wget https://storage.googleapis.com/gtex_external_datasets/eyegex_data/rna_seq_data/EyeGEx_retina_combined_genelevel_expectedcounts_byrid_nooutlier.tpm.matrix.gct
```
## Download and install software and tools

## Timing: 1.5 h

The goal of this section is to download and install software and tools for executing our pipeline. In cases where users encounter restriction or lack the necessary permissions to install software/tools on the machines, it is advisable to reach out to the machine administrator for assistance. The following command line operations can be executed as provided in most Linux distributions. We are using CentOS Linux release 7.9.2009. Of note, at the end of each code snippet, a final command is included to export to tool’s command to the filesystem environment, ensuring its availability in subsequent steps. 

8. Download and install vcflib.
```
>conda create -n vcflib
>conda install -c bioconda -n vcflib vcflib
```
9. Installing Python prerequisites. The scripts were primarily tested and are recommended to run using Python 3.8. Run command:
```
>pip=`which pip3`
>$pip install numpy==1.19.4
>$pip install tqdm==4.42.1
>$pip install scipy==1.3.3
>$pip install statsmodels==0.12.1
>$pip install rpy2==3.3.6
```
10. Installing R packages SKAT and ACAT for burden analysis.
```
>R_PATH=`which R`
>$R_PATH 
>install.packages("SKAT")
>install.packages("devtools")
>devtools::install_github("yaowuliu/ACAT")
```
11. Download and install bwa.
```
>cd $SOFTWARE_PATH
>mkdir packages
>cd packages
>wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.12.tar.bz2 
>tar -xvf bwa-0.7.12.tar.bz2 
>export BWA_PATH=$SOFTWARE_PATH/bwa-0.7.12
>cd $BWA_PATH  
>make
```
12. Download and install GATK. 
```
>cd $SOFTWARE_PATH
>wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip
>unzip gatk-4.3.0.0.zip
```
13. Download and install VEP.
```
>cd $SOFTWARE_PATH
>git clone https://github.com/Ensembl/ensembl-vep.git
>cd ensembl-vep
>perl INSTALL.pl
```
14. Download and install plink.
```
>cd $SOFTWARE_PATH
>wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip
>unzip plink_linux_x86_64_20230116.zip
```
15. Download and install EMMAX.
```
>cd $SOFTWARE_PATH
>wget http://csg.sph.umich.edu//kang/emmax/download/emmax-intel-binary-20120210.tar.gz
>tar xvf emmax-intel-binary-20120210.tar.gz
```
16. Download and install GCTA.
```
>cd $SOFTWARE_PATH
>wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip
>unzip gcta-1.94.1-linux-kernel-3-x86_64.zip
>cp gcta-1.94.1 $SOFTWARE_PATH/gcta64
```
17. Download and install vcftools.
```
>cd $SOFTWARE_PATH
>wget https://sourceforge.net/projects/vcftools/files/latest/download/ vcftools_0.1.13.tar.gz
>tar -xzf vcftools_0.1.13.tar.gz
>cd vcftools_0.1.13
>./configure --prefix=$SOFTWARE_PATH/vcftools_0.1.13
>make && make install
```
18. Download and install samtools.
```
>cd $SOFTWARE_PATH
>wget https://sourceforge.net/projects/samtools/files/samtools/1.9/samtools-1.9.tar.bz2
>tar -xvf samtools-1.9.tar.bz2
>cd samtools-1.9
>./configure --prefix=$SOFTWARE_PATH/samtools-1.9
>make && make install
```
19. Download and install bcftools. 
```
>cd $SOFTWARE_PATH
>wget https://sourceforge.net/projects/samtools/files/samtools/1.9/bcftools-1.9.tar.bz2
>tar -xvf bcftools-1.9.tar.bz2
>cd bcftools-1.9/
>./configure --prefix=$SOFTWARE_PATH/bcftools-1.9
>make && make install
```
20. Download and install sambamba
```
>cd $SOFTWARE_PATH
>wget https://github.com/biod/sambamba/releases/download/v0.6.6/sambamba_v0.6.6_linux.tar.bz2
>tar -xvf sambamba_v0.6.6_linux.tar.bz2
```
21. Installing WGCNA.
```
>R_PATH=`which R`
>$R_PATH
>install.packages("BiocManager") 
>BiocManager::install("WGCNA")
```
22. Download and install tabix and bgzip
```
>cd $SOFTWARE_PATH && mkdir htslib
>wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib1.10.2.tar.bz2
>tar -xvf htslib-1.10.2.tar.bz2
>cd htslib-1.10.2 && ./configure --prefix=$SOFTWARE_PATH/htslib
>make && make install
```
## Download the test dataset

Timing: 1.5 h

## In this section, the raw data for the demonstration of the protocol are retrieved. At the end of the process, the data should be in the proper place for the continuation of the protocol. 

23. Set the directory where the raw data will be placed and download the raw WES data.
```
>cd $PROJECT_PATH
>mkdir rawdata
>cd rawdata
>wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099967/SRR099967_1.fastq.gz 
>wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099967/SRR099967_2.fastq.gz 
> wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099969/SRR099969_1.fastq.gz
>wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099969/SRR099969_2.fastq.gz
>wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099957/SRR099957_1.fastq.gz >wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR099/SRR099957/SRR099957_2.fastq.gz
```
24. Download test file for exome-wide association and burden test.
```
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
25. Download RNAseq test file for WGCNA.
```
>mkdir 08_systems_genetics
>curl -o \
08_systems_genetics/datExpr1_combine_test.txt \
https://zenodo.org/record/8189201/files/datExpr1_combine_test.txt
> curl -o 08_systems_genetics/top50.txt \
https://zenodo.org/record/8189201/files/top50.txt
```
## Step-by-step method details
We prepared a GitHub repository containing all the necessary scripts for the method presented below. GitHub: https://github.com/sulab-wmu/MAGIC-PIPELINE
To follow the step-by-step instructions, users are encouraged to download the whole repository, including the example input files. In all the subsequent steps, the paths to the required software tools are the same as the “exported” paths in the respective command boxes under the “before you begin” section.
## Step 1: Data Pre-processing

Timing: hours to days; factors that affect timing include the number of samples and available computing resources

This section outlines the procedure of aligning the FASTQ pairs to the reference genome and collecting alignment statistics to ensure quality control. The output of this part comprises BAM files that are appropriate for the subsequent variant calling. At the end of the step, both BAM files and a report containing read alignment statistics are generated. 
1. Index the reference genome.
```
>cd $RESOURCES_PATH
>$BWA_PATH/bwa index ucsc.no_hap.hg19.fa
>$SAMTOOLS_PATH/samtools faidx ucsc.no_hap.hg19.fa
>$SAMTOOLS_PATH/samtools dict ucsc.no_hap.hg19.fa > ucsc.no_hap.hg19.fa
```
2. Alignment to the reference genome and mark duplicates in the BAM file.
```
>cd $PROJECT_PATH
>mkdir 01_mapping/
>software/bwa-0.7.12/bwa mem -t 32 -M -R "@RG\tID:SRR099967\tPL:illumina\tSM:SRR099967\tCN=WY" ref/ucsc.no_hap.hg19.fa rawdata/SRR099967_1.fastq.gz rawdata/SRR099967_2.fastq.gz | software/samtools-1.9/bin/samtools view -uS -F 4 - -o 01_mapping/SRR099967.bam && \
software/sambamba-0.6.6/sambamba sort -t 4 -u --tmpdir=temp 01_mapping/SRR099967.bam && \
software/sambamba-0.6.6/sambamba markdup --sort-buffer-size 262144 --hash-table-size 60000 -t 4 -l 0 --overflow-list-size 400000 --tmpdir=temp 01_mapping/SRR099967.sorted.bam  01_mapping/SRR099967.sort.rmdup.bam
```
3. Base Quality Score Recalibration and application on BAM files.
```
>cd $PROJECT_PATH
>export _JAVA_OPTIONS=-Djava.io.tmpdir=temp/
>software/gatk-4.3.0.0/gatk --java-options \
"-Dsamjdk.use_async_io_read_samtools=true \
-Dsamjdk.use_async_io_write_samtools=true \
-Dsamjdk.use_async_io_write_tribble=false -Xmx8g" \
BaseRecalibrator -R ref/ucsc.no_hap.hg19.fa \
-I 01_mapping/SRR099967.sort.rmdup.bam --verbosity WARNING \
--known-sites ref/ncbi.hg19_dbsnp151.vcf -O 01_realign/SRR099967.table
>software/gatk-4.3.0.0/gatk --java-options \
"-Dsamjdk.use_async_io_read_samtools=true \
-Dsamjdk.use_async_io_write_samtools=true \
-Dsamjdk.use_async_io_write_tribble=false -Xmx8g" \
ApplyBQSR -R ref/ucsc.no_hap.hg19.fa -I \
01_mapping/SRR099967.sort.rmdup.bam --verbosity WARNING -bqsr \
01_realign/SRR099967.table -O 01_realign/SRR099967.bam 
>software/sambamba-0.6.6/sambamba index 01_realign/SRR099967.bam
```
4.Collection of alignment statistics. 
If you are using another WES probe, please replace ref/Design_V2.hg19.bed with your own probe. In this step, several statistics related to the quality control of the alignment process are collected. A text file with statistics should be produced at the end of the process.
a. Total sequenced reads.
b. Aligned reads.
c. Uniquely aligned reads (q>20).
d. Reads overlapping targets.
e. Total sequenced bases.
f. Aligned bases.
g. Uniquely aligned bases.
h. Bases overlapping targets.
```
>cd $PROJECT_PATH
>RSCRIPT=`which Rscript`
>PERL=`which perl`
>$PERL software/QC_exome_mem.pl -R $RSCRIPT -i 01_mapping/SRR099967.sort.rmdup.bam -r ref/Design_V2.hg19.bed \
-c rawdata/SRR099967.R1.clean.fastq.gz -o 01_vcf/SRR099967 -s software/samtools-1.9/bin/samtools -b software/bedtools2/bin/bedtools -p software/picard.jar -t PE -plot \
-ref ref/ucsc.no_hap.hg19.fa -d ref/ucsc.no_hap.hg19.dict -C ref/Design_V2.hg19.bed
```
## Step 2: Variant Discovery	

Timing: hours to days

This section provides an overview of the variant calling procedure using GATK HaplotypeCaller. The BAM files generated during preprocessing are fed into GATK’s HaplotypeCaller to perform per-sample variant calling, outputting a genomic VCF (gVCF) for a single sample. Joint-calling and variant quality score recalibration (VQSR) generates a final, multi-sample VCFs files, serving as the starting point for downstream analysis.
5.For each sample use GATK HaplotypeCaller to create a gVCF callset file.
```
>cd $PROJECT_PATH
>mkdir 02_gvcf && mkdir 02_vcf
>export _JAVA_OPTIONS=-Djava.io.tmpdir=temp/
>software/gatk-4.3.0.0/gatk --java-options "-Xmx8G" HaplotypeCaller -R ref/ucsc.no_hap.hg19.fa -I 01_realign/SRR099967.bam -L ref/Design_V2.hg19.bed -ERC BP_RESOLUTION \
-O 01_vcf/SRR099967_HaplotypeCaller.g.vcf && \
software/bgzip 01_vcf/SRR099967_HaplotypeCaller.g.vcf && software/tabix -p vcf 01_vcf/SRR099967_HaplotypeCaller.g.vcf.gz
```
6.Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file.
```
>cd $PROJECT_PATH
>ls 01_vcf/*.gz > 02_gvcf/cohort.vcf.list
>software/gatk-4.3.0.0/gatk CombineGVCFs --java-options "-Xmx150g -Djava.io.tmpdir=temp/" -R ref/ucsc.no_hap.hg19.fa --variant 02_gvcf/cohort.vcf.list -O 02_gvcf/cohort.g.vcf
```
7.Perform joint-genotyping on the combined gVCF files.
```
>cd $PROJECT_PATH
>software/gatk-4.3.0.0/gatk GenotypeGVCFs --java-options "-Xmx150g -Djava.io.tmpdir=temp/" --variant 02_gvcf/cohort.g.vcf -new-qual true -R ref/ucsc.no_hap.hg19.fa -O 02_vcf/cohort.vcf.gz
```
8.Filter variants by hard-filtering before variant quality score recalibration (VQSR).
```
>cd $PROJECT_PATH
>software/vcftools/bin/vcftools --gzvcf 02_vcf/cohort.vcf.gz --max-missing 0.9 --recode --recode-INFO-all --out 03_filter/cohort.mis0.9
```










