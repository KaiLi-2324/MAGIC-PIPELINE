#!/usr/bin/bash


export _JAVA_OPTIONS=-Djava.io.tmpdir=temp/
software/gatk-4.3.0.0/gatk --java-options "-Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Xmx8g" \
BaseRecalibrator -R ref/ucsc.no_hap.hg19.fa -I 01_mapping/19BY0125.sort.rmdup.bam --verbosity WARNING --known-sites ref/ncbi.hg19_dbsnp147.vcf -O 01_realign/19BY0125.table && \
software/gatk-4.3.0.0/gatk --java-options "-Dsamjdk.use_async_io_read_samtools=true -Dsamjdk.use_async_io_write_samtools=true -Dsamjdk.use_async_io_write_tribble=false -Xmx8g" \
ApplyBQSR -R ref/ucsc.no_hap.hg19.fa -I 01_mapping/19BY0125.sort.rmdup.bam --verbosity WARNING -bqsr 01_realign/19BY0125.table -O 01_realign/19BY0125.bam && \
software/sambamba-0.6.6/sambamba index 01_realign/19BY0125.bam
