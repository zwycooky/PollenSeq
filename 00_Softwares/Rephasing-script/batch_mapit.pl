#!/usr/bin/perl

use strict;
use Cwd 'abs_path';
use Parallel::ForkManager;
use Getopt::Std;
our($opt_i,$opt_o,$opt_f,$opt_p,$opt_t);
getopts('i:o:f:p:t:');

my $inputdir = $opt_i;
my $outputdir = $opt_o;
my $ref_genome = $opt_f;
my $mapping_cpu = $opt_p;
my $cpu = $opt_t;

my $Usage = "\n\t$0 
	-i <input fastq dir>
	-o <out bam dir>
	-f <genome>
	-p <mapping threads>
	-t <threads>
\n";
die $Usage unless ($opt_i && $opt_o && $opt_f && $opt_p && $opt_t);

chomp(my $pwd = `pwd`);

if (!-e $outputdir) {
	mkdir $outputdir;
}

$inputdir = abs_path($inputdir);
$outputdir = abs_path($outputdir);
my $mapping_script = "02_mapit.pl";

## ForkManager setup ##
my $pm = Parallel::ForkManager->new($cpu);


chomp(my @fq = `find $inputdir -name "*.gz"`);
my $count = 0;
my ($fq1,$fq2);
foreach (sort @fq) {
	$count ++;
	if ($count < 2) {
		$fq1 = abs_path($_);
	}else{
		$fq2 = abs_path($_);
		my $id = (split /\./,(split /\//,$fq2)[-1])[0];
		$count = 0;
		$pm->start and next;
        	&mapping($fq1,$fq2,$id,$outputdir,$ref_genome,$mapping_cpu);
        	$pm->finish;
	}
}
$pm->wait_all_children();

sub mapping {
	my ($fq1,$fq2,$id,$outdir,$genome,$mapping_cpu) = @_;
	!system "$mapping_script -1 $fq1 -2 $fq2 -f $genome -b $genome -p $mapping_cpu -o $outdir/$id" or die "ERROR with fastp:$!";
}

