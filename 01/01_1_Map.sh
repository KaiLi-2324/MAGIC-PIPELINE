#!/usr/bin/bash


software/bwa-0.7.12/bwa mem -t 32 -M -R "@RG\tID:19BY16411\tPL:illumina\tSM:19BY16411\tCN=WY" ref/ucsc.no_hap.hg19.fa \
rawdata/19BY0125.R1.clean.fastq.gz rawdata/19BY0125.R2.clean.fastq.gz | samtools-1.9/bin/samtools view -uS -F 4 - -o 01_mapping/19BY0125.bam && \
software/sambamba-0.6.6/sambamba sort -t 4 -u --tmpdir=temp 01_mapping/19BY0125.bam && \
software/sambamba-0.6.6/sambamba markdup --sort-buffer-size 262144 --hash-table-size 60000 -t 4 -l 0 --overflow-list-size 400000 --tmpdir=temp 01_mapping/19BY0125.sorted.bam  01_mapping/19BY0125.sort.rmdup.bam
