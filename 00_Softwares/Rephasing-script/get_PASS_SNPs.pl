#!/usr/bin/perl
#
use strict;

my $vcf = @ARGV[0];
my $Usage = "\n\t$0 <VCF>
\n\tGet SNPs with PASS flag\n";
die $Usage unless (@ARGV == 1);

open VCF,'<',"$vcf" or die;
while (<VCF>) {
	chomp;
	if (/\A#/) {
		print "$_\n";
	}else{
		my $filter = (split)[6];
		if ($filter eq 'PASS') {
			print "$_\n";
		}
	}
}
close VCF;
