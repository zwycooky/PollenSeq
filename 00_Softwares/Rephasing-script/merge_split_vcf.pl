#!/usr/bin/perl
#
use strict;

my ($input,$output) = @ARGV[0,1];
my $Usage = "\n\t$0 <Raw_split_vcf.list> <output file>
\n";
die $Usage unless (@ARGV == 2);


my (%ex);
open OUT,'>',"$output" or die;
my $count = 0;
open IN,'<',"$input" or die;
while (<IN>) {
	chomp;
	my $path = $_;
	if ($count == 0){
		open VCF,'<',"$path" or die "Cannot open $path:$!";
		while (<VCF>) {
			chomp;
			print OUT "$_\n";
			unless (/\A#/) {
				my ($chr,$pos) = (split)[0,1];
				my $key = "$chr\t$pos";
				$ex{$key} = 1;
			}
		}
		close VCF;
		$count = 1;
	}else{
		open VCF,'<',"$path" or die "Cannot open $path:$!";
                while (<VCF>) {
                        chomp;
			if (/\A#/) { next };
			my ($chr,$pos) = (split)[0,1];
                        my $key = "$chr\t$pos";
			if (!exists $ex{$key}) {
                        	print OUT "$_\n";
				$ex{$key} = 1;
			}
                }
                close VCF;
	}
}
close IN;

close OUT;
