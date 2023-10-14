#!/usr/bin/perl -w
use warnings FATAL => 'all';
use Getopt::Long qw(:config no_ignore_case);
use Data::Dumper;
use File::Basename;

my ($R,$bamfile,$outdir,$regionfile,$Plot,$clean_data,$samtools,$bedtools,$type,$ref,$picard,$Core_gene,$ref_dict,$help);
GetOptions(
	"i:s"=>\$bamfile,
	"o:s"=>\$outdir,
	"r:s"=>\$regionfile,
	"c:s"=>\$clean_data,
	"s:s"=>\$samtools,
	"b:s"=>\$bedtools,
	"t:s"=>\$type,
	"R:s"=>\$R,
	"C:s"=>\$Core_gene,
	"ref:s"=>\$ref,
	"d:s"=>\$ref_dict,
	"p:s"=>\$picard,
	"plot"=>\$Plot,
	"h"=>\$help,
);

my $usage=<<USAGE;
usage:perl $0
		-i <bamfile>
		-r <region file>
		-c <clean_data>
		-o <outdir>
		-R <R path>
		-s <samtools path>
		-b <bedtools path>
		-t <'PE' or 'SE'>
                -p <picard path>
		-plot
                -ref <reference fasta path>
                -d <reference fasta dict>
		-h help
USAGE

$samtools ||= "/home/leiyoubing/software/samtools-0.1.19/samtools";
$bedtools ||= "/home/leiyoubing/software/bedtools2-2.19.1/bin";
$type ||= 'PE';
$ref ||= "/home/hongqy/workdir/database/ref/ucsc.no_hap.hg19.fa";
#$picard = "/share/public/software/picard-tools-1.57/picard-1.57.jar";
$picard ||= "/share/public/software/picard-2.10.7/picard.jar";
$Core_gene ||= "/home/hongqy/workdir/database/bed/Circle.V2.2.hg19.bed";
$ref_dict ||= "/home/hongqy/workdir/database/genome/hg19/ucsc.no_hap.hg19.dict";
$R ||= "/share/public/software/R/3.4.0/bin/R";
die $usage if (!$bamfile || $help || !$outdir || !$regionfile || !$clean_data);

mkdir ("$outdir", 0755) unless (-d $outdir);
mkdir ("$outdir/QC", 0755) unless(-d "$outdir/QC");
my $sample_name = basename $outdir;
my $Initial_bases_on_target=0;
my $Initial_bases_near_target=0;
my $Initial_bases_on_or_near_target=0;
my $Total_effective_yield=0;
my $Effective_sequences_on_target=0;
my $Effective_sequences_near_target=0;
my $Uniquely_mapped_to_target=0;
my $Effective_sequences_on_or_near_target=0;
my $Uniquely_mapped_to_genome=0;
my $Fraction_of_effective_bases_on_target=0;
my $Fraction_of_uniq_mapping_on_target=0;
my $Average_sequencing_depth_on_target=0;
my $Average_sequencing_depth_near_target=0;
my $Fraction_of_effective_bases_on_or_near_target=0;
my $Mismatch_rate_in_target_region=0;
my $Mismatch_rate_in_all_effective_sequence=0;
my $Base_covered_on_target=0;
my $Coverage_of_target_region=0;
my $Base_covered_near_target=0;
my $Coverage_of_flanking_region=0;
my $Fraction_of_target_covered_with_at_least_50x=0;
my $Fraction_of_target_covered_with_at_least_20x=0;
my $Fraction_of_target_covered_with_at_least_10x=0;
my $Fraction_of_target_covered_with_at_least_4x=0;
my $Mismatch_base_in_target_region=0;
my $Mismatch_base_in_all_effective_sequence=0;
my $clean_reads=0;
my $mapped_reads=0;
my $mapping_rate=0;
my $dup_reads=0;
my $unique_reads=0;
my $unique_rate=0;
my $dup_rate=0;
my $clean_bases=0;
#get clean reads number
my @r1 = split /,/,$clean_data;
foreach my $r1(@r1){
	if(-B $r1 ){
	    	my $tmp = `gzip -dc $r1|wc -l `;
		    $clean_reads += $tmp;
	}else{
		my $tmp = `wc -l $r1`;
		$tmp = $1 if ($tmp =~ /^\s*(\d+)\s*/);
		$clean_reads += $tmp;
	}
}
if ($type =~ /PE/i){
	$clean_reads=$clean_reads/2;
}elsif($type =~ /SE/i){
	$clean_reads=$clean_reads/4;
}
$clean_bases = $clean_reads * 150;


#get unique reads number of uniquely mapped to genome
my $key = (split /\./, basename $bamfile)[0];
open BAM,"$samtools view  -F 0x0004 $bamfile | " or die $!;
while(my $info = <BAM>){
	chomp($info);
	my @info=split /\s+/,$info;
	unless($info[1] & 0x800){ #0x800   SUPPLEMENTARY
		unless($info[1] & 0x100){ #0x100   SECONDARY
			unless($info[1] & 0x4){ #0x4     UNMAP
				$mapped_reads++;
				if($info[1] & 0x400){
					$dup_reads++;
					next;
				}
				if($info[4] >= 5)
				{
					$Uniquely_mapped_to_genome++;
				}
				if($info=~/NM:i:(\d+)/)
				{
					my $temp_m = $1;
					$temp_m -= $1 if($info[5] =~ /(\d+)[ID]/);
					$Mismatch_base_in_all_effective_sequence+=$temp_m;
				}
			}
		}
	}
}
close BAM;

$unique_reads = $mapped_reads - $dup_reads;
$unique_rate = sprintf("%.2f%%", 100*$unique_reads/$mapped_reads);
$mapping_rate=$mapped_reads/$clean_reads;
$mapping_rate=sprintf("%.2f%%",100*$mapping_rate);
$dup_rate=$dup_reads/$mapped_reads;
$dup_rate=sprintf("%.2f%%",100*$dup_rate);

open REG,"$regionfile" or die $!;
open FREG,">$outdir/QC/freg_tmp.txt" or die $!;
while(<REG>)
{
	chomp;
	next if(/^\s*$/ || /^#/);
	my @info=split (/\t/,$_);
	next if (@info < 3);
	$Initial_bases_on_target+=$info[2]-$info[1];
	my $pri_beg_pos=$info[1]-200;
	my $pri_end_pos=$info[1];
	my $next_beg_pos=$info[2];
	my $next_end_pos=$info[2]+200;
	$pri_beg_pos=0 if($pri_beg_pos < 0);
	if ($info[1] <= 1)
	{
		print FREG "$info[0]\t$next_beg_pos\t$next_end_pos\n";
	}
	else
	{
		print FREG "$info[0]\t$pri_beg_pos\t$pri_end_pos\n";
		print FREG "$info[0]\t$next_beg_pos\t$next_end_pos\n";
	}
}
close(FREG);
close(REG);

`for i in \$(seq 1 22) X Y;do grep -w chr\$i $outdir/QC/freg_tmp.txt|sort -k2n >>$outdir/QC/tmp.txt;done`;
`mv $outdir/QC/tmp.txt $outdir/QC/freg_tmp.txt`;
`$bedtools  merge  -i $outdir/QC/freg_tmp.txt >$outdir/QC/tmp.txt`;
`$bedtools subtract -a $outdir/QC/tmp.txt -b $regionfile >$outdir/QC/freg.txt`;
`rm $outdir/QC/freg_tmp.txt $outdir/QC/tmp.txt`;

$Initial_bases_near_target=`awk '{total+=\$3-\$2};END{print total}' $outdir/QC/freg.txt`;
chomp($Initial_bases_near_target);

`$samtools depth -q 20 -b  $regionfile $bamfile >$outdir/QC/target.depth`;
`$samtools depth -q 20 -b $outdir/QC/freg.txt $bamfile >$outdir/QC/flanking.depth`;

$Effective_sequences_on_target=`awk '{total+=\$3};\$3 != 0 {cov+=1};END{print total"\t"cov}' $outdir/QC/target.depth`;
chomp($Effective_sequences_on_target);
($Effective_sequences_on_target,$Base_covered_on_target)=split(/\t/,$Effective_sequences_on_target);

$Effective_sequences_near_target=`awk '{total+=\$3};\$3 != 0 {cov+=1};END{print total"\t"cov}' $outdir/QC/flanking.depth`;
chomp($Effective_sequences_near_target);
if($Effective_sequences_near_target!~/\d/){$Effective_sequences_near_target="0\t0";}
($Effective_sequences_near_target,$Base_covered_near_target)=split (/\t/,$Effective_sequences_near_target);

$Initial_bases_on_or_near_target=$Initial_bases_on_target+$Initial_bases_near_target;

`$samtools depth -q 20 $bamfile >$outdir/QC/whole_genome.depth`;
$Total_effective_yield=`awk '{total+=\$3};END{print total}' $outdir/QC/whole_genome.depth`;
chomp $Total_effective_yield;


$Fraction_of_effective_bases_on_target=$Effective_sequences_on_target/$Total_effective_yield;
$Effective_sequences_on_or_near_target=$Effective_sequences_on_target + $Effective_sequences_near_target;
$Fraction_of_effective_bases_on_or_near_target=$Effective_sequences_on_or_near_target/$Total_effective_yield;

open TMP,"$samtools view  -F 0x4 -L $regionfile $bamfile | " or die $!;
while(<TMP>)
{
	chomp;
	my @arry = split /\t/, $_;
	my $fflag=$arry[1];
	unless($fflag & 0x800){ #0x800   SUPPLEMENTARY
		unless($fflag & 0x100){ #0x100   SECONDARY
			unless($fflag & 0x400){ #0x400     DUPLICATE
				if($arry[4] >= 5)
				{
					$Uniquely_mapped_to_target++;
				}
				my $temp_m = 0;
				if($_=~/NM:i:(\d+)/){
					$temp_m = $1;
					$temp_m -= $1 if($arry[5] =~ /(\d+)[ID]/);
				}
				$Mismatch_base_in_target_region+=$temp_m;
			}
		}
	}
}
close(TMP);

$Mismatch_rate_in_target_region=$Mismatch_base_in_target_region/$Effective_sequences_on_target;
$Mismatch_rate_in_all_effective_sequence=$Mismatch_base_in_all_effective_sequence/$Total_effective_yield;
$Fraction_of_uniq_mapping_on_target = $Uniquely_mapped_to_target/$Uniquely_mapped_to_genome;

$Average_sequencing_depth_on_target=$Effective_sequences_on_target/$Initial_bases_on_target;
$Average_sequencing_depth_near_target=$Effective_sequences_near_target/$Initial_bases_near_target;
	
$Coverage_of_target_region=$Base_covered_on_target/$Initial_bases_on_target;
$Coverage_of_flanking_region=$Base_covered_near_target/$Initial_bases_near_target;

my $tmp1=`awk '\$3 >=20 {total1++};\$3 >=10 {total2++};\$3 >=4 {total3++};\$3 >=50 {total4++};END{print total1"\t"total2"\t"total3"\t"total4}' $outdir/QC/target.depth`;
chomp($tmp1);
my @info1;
@info1=split /\t/,$tmp1;
$info1[0]=0 unless($info1[0] =~ /\d+/);
$info1[1]=0 unless($info1[1] =~ /\d+/);
$info1[2]=0 unless($info1[2] =~ /\d+/);
$info1[3]=0 unless($info1[3] =~ /\d+/);
$Fraction_of_target_covered_with_at_least_50x=$info1[3]/$Initial_bases_on_target;
$Fraction_of_target_covered_with_at_least_20x=$info1[0]/$Initial_bases_on_target;
$Fraction_of_target_covered_with_at_least_10x=$info1[1]/$Initial_bases_on_target;
$Fraction_of_target_covered_with_at_least_4x=$info1[2]/$Initial_bases_on_target;
my $name=basename($bamfile);
$name=~s/\.sort\.rmdup\.bam//g;
my $sample=$name;
`$samtools view  -f66 $bamfile |cut -f 9|sed 's/^-//' >$outdir/QC/$sample.insert.txt`;
Plot_insert_size($outdir,"$outdir/QC/$sample.insert.txt",$sample);
my $mean_insert_size=` samtools view  -f66 $bamfile |awk '{if (\$9 > 0 && \$9 < 500) {S+=\$9; T+=1}}END{print  S/T}' `;


open STAT,">$outdir/QC/$key\_QC.xls" or die $!;
print STAT "Sample\t$sample\n";
print STAT "Clean bases\t$clean_bases\n";
print STAT "Clean reads\t$clean_reads\n";
print STAT "Mapped reads\t$mapped_reads($mapping_rate)\n";
print STAT "Duplicate reads\t$dup_reads($dup_rate)\n";
print STAT "Unique reads\t$unique_reads($unique_rate)\n";
print STAT "Reads uniquely mapped to target\t$Uniquely_mapped_to_target\n";
print STAT "Reads uniquely mapped to genome\t$Uniquely_mapped_to_genome\n";
print STAT "Total bases on target\t$Initial_bases_on_target\n";
print STAT "Total bases near target\t$Initial_bases_near_target\n";
printf STAT "Total sequences on target(Mb)\t%.2f\n",$Effective_sequences_on_target/1000000;
printf STAT "Total sequences near target(Mb)\t%.2f\n",$Effective_sequences_near_target/1000000;
printf STAT "Total effective yield(Mb)\t%.2f\n",$Total_effective_yield/1000000;
printf STAT "Fraction of effective bases on target\t%.2f%%\n",100*$Fraction_of_effective_bases_on_target;
printf STAT "Fraction of effective bases on or near target\t%.2f%%\n",100*$Fraction_of_effective_bases_on_or_near_target;
printf STAT "Fraction of uniquely mapped on target\t%.2f%%\n",100*$Fraction_of_uniq_mapping_on_target;
printf STAT "Average sequencing depth on target\t%.2f\n",$Average_sequencing_depth_on_target;
printf STAT "Average sequencing depth near target\t%.2f\n",$Average_sequencing_depth_near_target;
printf STAT "Average insert size of the library\t%.2f\n",$mean_insert_size;

print STAT "Base covered on target\t$Base_covered_on_target\n";
printf STAT "Coverage of target region\t%.2f%%\n",100*$Coverage_of_target_region;
print STAT "Base covered near target\t$Base_covered_near_target\n";
printf STAT "Coverage of near target\t%.1f%%\n",100*$Coverage_of_flanking_region;
printf STAT "Mismatch rate in target\t%.2f%%\n",100*$Mismatch_rate_in_target_region;
printf STAT "Mismatch rate in total\t%.2f%%\n",100*$Mismatch_rate_in_all_effective_sequence;
printf STAT "Fraction of target covered at least 4x\t%.2f%%\n",100*$Fraction_of_target_covered_with_at_least_4x;
printf STAT "Fraction of target covered at least 10x\t%.2f%%\n",100*$Fraction_of_target_covered_with_at_least_10x;
printf STAT "Fraction of target covered at least 20x\t%.2f%%\n",100*$Fraction_of_target_covered_with_at_least_20x;
printf STAT "Fraction of target covered at least 50x\t%.2f%%\n",100*$Fraction_of_target_covered_with_at_least_50x;
close STAT;

open DB,"<$outdir/QC/target.depth" or die $!;
open DF,">$outdir/QC/depth_frequency.xls" or die $!;
my %depth=();
while(<DB>)
{
	chomp;
	my @tmp=split;
	$depth{$tmp[2]}++;
}
close(DB);
my $maxCov=0;
#$depth{"0"}+=$Initial_bases_on_target-$Base_covered_on_target;
	
foreach my $depth (sort {$a<=>$b} keys %depth)
{
	next if($depth==0);
	my $per=$depth{$depth}/$Initial_bases_on_target;
	$maxCov = $per if($per > $maxCov);
	print DF "$depth\t$per\t$depth{$depth}\n";
}
close(DF);
	
open CU,">$outdir/QC/cumu.xls" or die $!;
print CU "Depth\tTRPercent\n";
my @depth= sort {$a<=>$b} keys %depth;
foreach my $depth1 (sort {$a<=>$b} keys %depth)
{
	my $tmp=0;
	next if($depth1==0);
	foreach my $depth2 (@depth)
	{
		if($depth2 >= $depth1)
		{
			$tmp+=$depth{$depth2};
		}
	}
	$tmp = $tmp/$Initial_bases_on_target;
	print CU "$depth1\t$tmp\n";
}
close(CU);

if($Plot)
{
	my $ylim = 100*$maxCov;
	my ($xbin,$ybin);
	$ylim= int($ylim) + 1;
	if($ylim <= 3)
	{
		$ybin = 0.5;
	}else{
		$ybin=1;
	}
	my $xlim=0;
	if($Average_sequencing_depth_on_target<30)
	{
		$xlim=100;
		$xbin=20;
	}elsif($Average_sequencing_depth_on_target < 50)
	{
		$xlim=160;
		$xbin=20;
	}elsif($Average_sequencing_depth_on_target < 120)
	{
		$xlim=250;
		$xbin=50;
	}else{
		$xlim=600;
		$xbin=100;
	}
	histPlot($outdir,"$outdir/QC/depth_frequency.xls",$ylim,$ybin,$xlim,$xbin, $key);
	cumuPlot($outdir,"$outdir/QC/cumu.xls",$xlim,$xbin, $key);
}

my %regionfile_depth;
open IN,"$outdir/QC/target.depth" or die "Can't open the $outdir/QC/target.depth:$!";
while (<IN>)
{
	chomp;
	next if (/^\s*$/);
	my @tem = split /\t/;
	next unless ($tem[0] =~/^chr\d+$/i || $tem[0]=~/^chrX$/i || $tem[0]=~/^chrY$/i);
	$regionfile_depth{$tem[0]."\t".$tem[1]}= $tem[2];
}
close IN;

my (%chr_depth,%chr_cover,%chr_len) = ();
open RE,$regionfile || die "$!";
while (my $line = <RE>){
	chomp($line);
	next if ($line =~ /^\s*$/);
	my @tem = split /\t/,$line;

	my $regionfile_chr="$tem[0]";
	my $regionfile_len = $tem[2]-$tem[1]+1;
	my ($regionfile_depth,$regionfile_num,$regionfile_dep) = (0,0,0);
	for (my $i=$tem[1];$i<=$tem[2];$i++){
		my $search = $tem[0]."\t".$i;
		if (exists $regionfile_depth{$search}){
			$regionfile_dep = $regionfile_depth{$search};
			$regionfile_num++;
		}else{
			$regionfile_dep=0;
		}
		$regionfile_depth += $regionfile_dep;
		push @{$regionfile_chr},$regionfile_dep;
	}
	$chr_depth{$regionfile_chr} += $regionfile_depth;
	$chr_cover{$regionfile_chr} += $regionfile_num;
	$chr_len{$regionfile_chr} += $regionfile_len;
}
close RE;

open OUT,">$outdir/QC/dist.xls" or die $!;

my (@Cover,$Cover_median);
foreach my $index (1..22, 'X','Y'){
	my $chr = "chr$index";
	my $cover = $chr_depth{$chr};
	my $len = $chr_len{$chr};
	$cover = sprintf ("%.2f",($cover/$len));
	push @Cover,$cover;
}
$Cover_median = Median(@Cover);

foreach my $tmp (1..22,'X','Y'){
	my $chr = "chr$tmp";
	my $cover = $chr_cover{$chr};
	my $Depth = $chr_depth{$chr};
	my $len = $chr_len{$chr};
	$cover = sprintf ("%.2f",($cover*100/$len));
	$Depth = sprintf ("%.2f",($Depth/$len));
	my $ratio = $Depth/$Cover_median;
	my $ploidy = $ratio * 2;
	my $ploidy_norm;
	if (abs($ploidy -2) <= 0.8 ){
		$ploidy_norm = 2;
	}
	elsif (($ploidy -2) > 0.8 ){
		$ploidy_norm = 3;
	}
	elsif(($ploidy -2) < -0.8 && ($ploidy -2) > -1.7){
		$ploidy_norm = 1;
	}	
	elsif(($ploidy -2) <= -1.7){
		$ploidy_norm = 0;
	}
	
	@{$chr} = sort {$a<=>$b} @{$chr};
	my ($number,$number1,$number2,$median);
	if (@{$chr}%2==1){
		$number = (@{$chr}-1)/2;
		$median = sprintf("%.2f", @{$chr}[$number]);
	}else{
		$number1 = (@{$chr}/2)-1;
		$number2 = @{$chr}/2;
		$median = sprintf ("%.2f",(@{$chr}[$number1]+@{$chr}[$number2])/2);
   	}
	print OUT "$chr\t$Depth\t$cover%\t$median\t$ploidy_norm\t$ploidy\n";
}
close OUT;

sub Median{
	my (@Cover) = @_;
	my @vals = sort {$a <=> $b} @Cover;
	my $len = @vals;
	if($len%2){
		return $vals[int($len/2)];
	}
	else{
		return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
	}
}

distPlot($outdir,"$outdir/QC/dist.xls",$sample);

GC_bias($sample);

Gene_QC($sample);

#`rm $outdir/QC/freg.txt $outdir/QC/whole_genome.depth $outdir/QC/target.depth $outdir/QC/flanking.depth`; //annot by hqy 2017-11-11

sub cumuPlot {
        my ($outdir, $dataFile, $xlim, $xbin, $sample) = @_;
        my $figFile = "$outdir/QC/$sample\_cumuPlot.pdf";
        my $Rline=<<Rline;
        pdf(file="$figFile",w=8,h=6)
        rt <- read.table("$dataFile",header=T)
        opar <- par()
        x <- rt\$Depth[1:($xlim+1)]
        y <- 100*rt\$TRPercent[1:($xlim+1)]
        par(mar=c(4.5, 4.5, 2.5, 2.5))
        plot(x,y,col="red",type='l', lwd=2, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
        xpos <- seq(0,$xlim,by=$xbin)
        ypos <- seq(0,100,by=20)
        axis(side=1, xpos, tcl=0.2, labels=FALSE)
        axis(side=2, ypos, tcl=0.2, labels=FALSE)
        mtext("Cumulative sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
        mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
        mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.4)
        mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.4)
        par(opar)
        dev.off()
        png(filename="$outdir/QC/$sample\_cumuPlot.png",width = 480, height = 360)
        par(mar=c(4.5, 4.5, 2.5, 2.5))
        plot(x,y,col="red",type='l', lwd=3, bty="l",xaxt="n",yaxt="n", xlab="", ylab="", ylim=c(0, 100))
        xpos <- seq(0,$xlim,by=$xbin)
        ypos <- seq(0,100,by=20)
        axis(side=1, xpos, tcl=0.2, labels=FALSE)
        axis(side=2, ypos, tcl=0.2, labels=FALSE)
        mtext("Cumulative sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
        mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
        mtext(xpos, side=1, las=1, at=xpos, line=0.3, cex=1.5)
        mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.5)
        par(opar)
        dev.off()

Rline
        open (ROUT,">$figFile.R");
        print ROUT $Rline;
        close(ROUT);

        system("$R CMD BATCH $figFile.R");
}


sub histPlot {
	my ($outdir, $dataFile, $ylim, $ybin, $xlim, $xbin, $sample) = @_;
	my $figFile = "$outdir/QC/$sample\_histPlot.pdf";
	my $Rline=<<Rline; 
	pdf(file="$figFile",w=8,h=6)
	rt <- read.table("$dataFile")
	t=sum(rt\$V2[($xlim+1):length(rt\$V2)])
	y=c(rt\$V2[1:$xlim],t)
	y <- y*100
	x <- rt\$V1[1:($xlim+1)]
	#add by 2018.05.10
	#complement percent of seqeuncing depth (0).  
	y_sum <- sum(y)
	y_compare <- 100-y_sum
	x <- c(0,x)
	y <- c(y_compare,y) 
	data <- data.frame(x,y)
	#retrieve 25% and 75% index
	index_low <- function(y){
		per <- 0
		for(i in 1:length(y)){
			per <- per + y[i]
			if(per >= 25){
			return(i)
			break
			}
		}
	}
	index_high <- function(y){
		per <- 0
		for(i in 1:length(y)){
			per <- per + y[i]
			if(per >= 75){
			return(i)
			break
			}
		}
	}
	###full width at half maximum:FWHM
	half_low <- function(y,y_max){
		per <- 0
		for(i in 2:index_max){
			if(y[i] >= (y_max/2)){
			return(i)
			break
			}
		}
	}

	half_high <- function(y,y_max){
		per <- 0
		for(i in index_max:$xlim){
			if(y[i] <= (y_max/2)){
			return(i)
			break
			}
		}
	}
	############################################################	
	X_low <- index_low(y)
	X_high <- index_high(y) - 1
	y_max <- max(y[X_low:X_high])
	index_max <- which(y == y_max)
	###75%-25% (Peak/Width)
	ratio <- 100*(y_max)/(X_high-X_low)
	x_a <- X_high - X_low
	#FWHM
	x_half_low <- half_low(y,y_max)
	x_half_high <- half_high(y,y_max) -1
	x_b <- x_half_high-x_half_low
	############################################################
	opar <- par()
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="gray",type='h', lwd=1.5, xaxt="n",yaxt="n", xlab="", ylab="", bty="l",ylim=c(0,$ylim),xlim=c(0,$xlim))
	with(data, polygon(x=c(x[c(X_low,X_low:X_high,X_high)]), y= c(0, y[X_low:X_high], 0), col="dodgerblue",border =NA))
	with(data, segments(x_half_low,y_max/2,x_half_high,y_max/2,col="firebrick1",lty="dashed",lwd="2"))
	#tag <- paste("FWHM",x_b,sep = "=")
	tag <- paste("IQR",x_a,sep = "=")
	text((x_half_low+x_half_high)/2,y_max/2 + 0.1,tag)
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,$ylim,by=$ybin)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	end <- length(xpos)-1
	mtext(c(xpos[1:end],"$xlim+"), side=1, las=1, at=xpos, line=0.3, cex=1.4)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.4)
	par(opar)
	dev.off()
	png(filename="$outdir/QC/$sample\_histPlot.png",width = 480, height = 360)
	par(mar=c(4.5, 4.5, 2.5, 2.5))
	plot(x,y,col="gray",type='h', lwd=1.5, xaxt="n",yaxt="n", xlab="", ylab="", bty="l",ylim=c(0,$ylim),xlim=c(0,$xlim))
	with(data, polygon(x=c(x[c(X_low,X_low:X_high,X_high)]), y= c(0, y[X_low:X_high], 0), col="dodgerblue",border =NA))
	with(data, segments(x_half_low,y_max/2,x_half_high,y_max/2,col="firebrick1",lty="dashed",lwd="2"))
	#tag <- paste("FWHM",x_b,sep = "=")
	tag <- paste("IQR",x_a,sep = "=")
	text((x_half_low+x_half_high)/2,y_max/2 + 0.1,tag)
	xpos <- seq(0,$xlim,by=$xbin)
	ypos <- seq(0,$ylim,by=$ybin)
	axis(side=1, xpos, tcl=0.2, labels=FALSE)
	axis(side=2, ypos, tcl=0.2, labels=FALSE)
	mtext("Sequencing depth",side=1, line=2, at=median(xpos), cex=1.5 )
	mtext("Fraction of target bases (%)",side=2, line=3, at=median(ypos), cex=1.5 )
	end <- length(xpos)-1
	mtext(c(xpos[1:end],"$xlim+"), side=1, las=1, at=xpos, line=0.3, cex=1.5)
	mtext(ypos, side=2, las=1, at=ypos, line=0.3, cex=1.5)
	par(opar)
	dev.off()
	Ratio <- y_max/x_b
	IQR <- x_a
	FWHM <- x_b
	info <- rbind(IQR,FWHM,Ratio,ratio,index_max)
	rownames(info) <- c("IQR","FWHM","Peak/FWHM","Peak/IQR Width","Max_Depth")
	write.table(info,"$outdir/QC/$sample\_FWHM.stat.txt",col.names = FALSE,row.names=TRUE)
	
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);
	system("$R CMD BATCH $figFile.R");
	#system("rm  $figFile.R  $figFile.Rout");
}

sub Plot_insert_size {
	my ($outdir, $dataFile,$sample) = @_;
	my $figFile = "$outdir/QC/$sample\_insert_size.pdf";
	my $Rline=<<Rline; 
	data <- scan("$dataFile")
	pdf("$figFile",width=8)
	if(range(data)[2] > 1000){
        	plot(table(data), xlab="Insert Size", ylab="Reads Number",xlim=c(0,1000))
        	data <- data[data<1000]
    		hist(data, xlim=c(0,1000))
	}else{
        	plot(table(data), xlab="Insert Size", ylab="Reads Number")
        	data <- data[data<1000]
		hist(data)
	}
	dev.off()
	png(filename="$outdir/QC/$sample\_insert_size.png",width = 480, height = 360)
	if(range(data)[2] > 1000){
		plot(table(data), xlab="Insert Size", ylab="Reads Number",xlim=c(0,1000))
		data <- data[data<1000]
		hist(data, xlim=c(0,1000))
	}else{
		plot(table(data), xlab="Insert Size", ylab="Reads Number")
		data <- data[data<1000]
		hist(data)
	}
	dev.off()
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);
	system("$R CMD BATCH $figFile.R");
}

sub distPlot
{
	my ($outdir, $dataFile,$sample) = @_;
	my $figFile = "$outdir/QC/$sample\_reads_distribution.png";
	my $main = "$sample Read Distribution; Karyotype";
	my $Rline=<<Rline;
	library("plotrix")
	dat<-read.table("$dataFile",header=F,row.names=1)
	ploidy=sum(dat[,4])
	if(dat[23,4] == 2){
		sex.x="XX"
	}
	if(dat[23,4] == 1){
		sex.x="X"
	}
	if(dat[23,4] == 0){
		sex.x=""
	}
	
	if(dat[24,4] == 1){
		sex.y="Y"
	}
	if(dat[24,4] == 0){
		sex.y=""
	}
	sex = paste(sex.x,sex.y,sep="")
	dat[,2]=gsub("%","",dat[,2])
	rownames(dat)=gsub("chr","",rownames(dat))
	index <-as.numeric(rownames(dat))[!is.na(as.numeric(rownames(dat)))]
	d_s1<-dat[pmatch(sort(index),rownames(dat)),]
	index <-which(is.na(as.numeric(rownames(dat)))) 
	d_s2<-dat[pmatch(sort(rownames(dat)[index]),rownames(dat)),]
	dat <- rbind(d_s1,d_s2) 
	rownames(dat)<-c(rownames(d_s1),rownames(d_s2))
	options(bitmapType='cairo')
	png(filename="$figFile", width = 850, height = 600)
	par(mar=c(5.1,4.1,4.1,4.1))
	xpos=seq(1,nrow(dat))
	b <- min(as.numeric(rownames(dat)))
	xxpos <- paste("chr",rownames(dat),sep="")
	ymax=max(c(dat[,1], dat[,3]))

	karyo=paste(ploidy,sex,sep=",")
	main=paste("$main",karyo,sep=":")
	twoord.plot(xpos, as.numeric(dat[,1]), xpos, as.numeric(dat[,2]), lylim=c(0, ymax*1.4), rylim=c(0,130), type=c("bar", "b"), lcol="blue", halfwidth=0.2, main=main, xticklab="")
	points(xpos,dat[,3],type="o",pch=20,col="purple")
	mtext("Depth",side=2,line=2,font=2)
	mtext("Coverage(%)",side=4,line=2,font=2)
	text(xpos,-(ymax/50),labels=xxpos,srt=45,adj=1,xpd=T)
	legend("topleft",legend=c("Mean depth","Coverge","Median Depth"),col=c("blue","red","purple"),lty=1,box.lty=0)
	dev.off()
Rline
	open (ROUT,">$figFile.R");
	print ROUT $Rline;
	close(ROUT);
	system("$R CMD BATCH $figFile.R");
}

sub GC_bias{
	my ($sample) = @_;
	my $GC =
		"java -Xmx6G -jar $picard CollectGcBiasMetrics " .
		"R=$ref " .
		"I=$bamfile " . 	
		"SCAN_WINDOW_SIZE=100 " .
		"O=$outdir/QC/$sample.gc_bias_metrics.txt " .
		"CHART=$outdir/QC/$sample.gc_bias_metrics.pdf " .
		"S=$outdir/QC/$sample.summary_metrics.txt";
	system($GC) == 0 or die();

}

sub Gene_QC{
	my ($sample) = @_;
	my $BedToInterval =
		"java -Xmx6G -jar $picard BedToIntervalList " .
		"I=$Core_gene " .
		"O=$outdir/QC/Core.gene.interval " .
		"SD=$ref_dict";
	system($BedToInterval) == 0 or die();

	my $OnTargetCore = 
		"java -Xmx6G -jar $picard CollectHsMetrics " .
		"R=$ref " .
		"I=$bamfile " . 	
		"O=$outdir/QC/$sample.Core.gene.coverage.txt " .
		"BI=$outdir/QC/Core.gene.interval " .
		"TI=$outdir/QC/Core.gene.interval";
	system($OnTargetCore) == 0 or die();

	my $filter = 
		"egrep \"^BAIT|^Core\" $outdir/QC/$sample.Core.gene.coverage.txt | awk \'{for(i=1;i<=NF;i++)if(NR == 1){a[i]=\$i;}else{a[i]=a[i]\"\t\"\$i;{print a[i]}}}\' - | grep \"PCT_TARGET\" > $outdir/QC/$sample.Core.gene.xls";
	system($filter) == 0 or die();
}


