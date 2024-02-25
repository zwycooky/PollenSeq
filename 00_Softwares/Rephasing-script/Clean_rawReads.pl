#!/usr/bin/perl

use strict;
use Cwd 'abs_path';
use Parallel::ForkManager;

my ($inputdir,$outputdir,$fastp_cpu, $cpu) = @ARGV[0,1,2,3];
my $Usage = "\n\t$0 <input fastq dir> <out clean fastq dir> <fastp threads> <threads>
\n";
die $Usage unless (@ARGV == 4);

chomp(my $pwd = `pwd`);

if (!-e $outputdir) {
	mkdir $outputdir;
}

$inputdir = abs_path($inputdir);
$outputdir = abs_path($outputdir);

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
        &fastp($fq1,$fq2,$id,$outputdir,$fastp_cpu);
        $pm->finish;
	}
}
$pm->wait_all_children();

sub fastp {
	my ($fq1,$fq2,$id,$outdir,$fastp_cpu) = @_;
	!system "fastp -i $fq1 -I $fq2 -o $outdir/$id.R1.clean.fq.gz -O $outdir/$id.R2.clean.fq.gz --thread=$fastp_cpu" or die "ERROR with fastp:$!";
}
