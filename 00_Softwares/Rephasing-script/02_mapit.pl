#!/usr/bin/perl

use strict;
use Getopt::Std;
our($opt_1,$opt_2,$opt_f,$opt_b,$opt_p,$opt_o);
getopts('1:2:f:b:p:o:');

my $left_reads = $opt_1;
my $right_reads = $opt_2;
my $reference = $opt_f;
my $bwa_index = $opt_b;
my $threads = (defined $opt_p)?$opt_p:6;
my $output = $opt_o;

my $Usage = "\n$0 -1 -2 -f <options>
        -1 <FILE>  left reads
        -2 <FILE>  right reads
        -f <FILE>  reference genome
        -b <INDEX> BWA index
        -p <INT>   number of threads [6]
        -o <STR>   output prefix
        \n";
die $Usage unless ($opt_1 && $opt_2 && $opt_f && $opt_b && $opt_o);

#######################################################
# my software #
my $bwa = "bwa";
my $java = "java";
my $picard = "picard.jar";
my $samtools = "samtools";


# checking the index of reference genome #
if (!-e "$reference.fai") {
        !system "$samtools faidx $reference" or die;
}
my $prefix = (split /\.fa/,$reference)[0];
if (!-e "$prefix.dict") {
        !system "$java -jar $picard CreateSequenceDictionary R=$reference O=$prefix.dict" or die;
}


# get read group #
chomp(my $header = `zcat $left_reads |head -1`);
my ($id1,$id2) = (split /:/,$header)[2,3];
my $id = "$id1.$id2";

my $sample_id = (split /_/,(split /\//,$output)[-1])[0];
!system "$bwa mem -t $threads -M -R \"\@RG\\tID:bwa\\tPL:Illumina\\tLB:1\\tSM:$sample_id\\tPU:bwa\" $bwa_index $left_reads $right_reads | samtools view -@ $threads -o $output.bam -" or die "Error with bwa mem:$!";
print "$bwa mem -t $threads -M -R \"\@RG\tID:bwa\tPL:Illumina\tLB:1\tSM:$sample_id\tPU:bwa\" $bwa_index $left_reads $right_reads | samtools view -@ $threads -o $output.bam -";

!system "export JAVA_TOOL_OPTIONS=-Xmx30g";
## picard ##
# save the raw bam file #
!system "$java -jar $picard SortSam I=$output.bam O=$output.sorted.bam SORT_ORDER=coordinate" or die "Error with SortSam:$!";
unlink "$output.bam";

!system "$java -jar $picard MarkDuplicates I=$output.sorted.bam O=$output.sorted.marked.duplicates.bam M=$output.marked_dup_metrics.txt REMOVE_DUPLICATES=true" or die "Error with MarkDuplicates:$!";
unlink "$output.sorted.bam";

## index bam ##
!system "$samtools index $output.sorted.marked.duplicates.bam" or die "ERROR with index:$!";
