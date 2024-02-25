#!/usr/bin/perl
#
use strict;

my ($vcf,$out) = @ARGV[0,1];
my $Usage = "\n\t$0 <VCF> <OUT>
\n\tGet SNPs with PASS flag\n";
die $Usage unless (@ARGV == 2);

open OUT,'>',"$out" or die;
open VCF,'<',"$vcf" or die;
while (<VCF>) {
	chomp;
	if (/\A#/) {
		print OUT "$_\n";
	}else{
		my $filter = (split)[6];
		if ($filter eq 'PASS') {
			print OUT "$_\n";
		}
	}
}
close VCF;
close OUT;

