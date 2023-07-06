#!/usr/bin/bash


/usr/bin/perl software/QC_exome_mem.pl -R /share/apps/R-3.5.1/bin/R -i 01_mapping/19BY0125.sort.rmdup.bam -r ref/Design_V2.hg19.bed \
-c rawdata/19BY0125.R1.clean.fastq.gz -o 01_vcf/19BY0125 -s software/samtools-1.9/bin/samtools -b software/bedtools2/bin/bedtools -p software/picard.jar -t PE -plot \
-ref ref/ucsc.no_hap.hg19.fa -d ref/ucsc.no_hap.hg19.dict

# delete the .Rout files
files=$(find . -name "*.Rout" 2> /dev/null | wc -l)
if [ "$files" != "0" ]
then
        for f in *.Rout
        do
                [[ -e "$f" ]] || break
                rm "$f"
        done
fi
