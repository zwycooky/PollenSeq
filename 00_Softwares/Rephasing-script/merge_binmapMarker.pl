#!/usr/bin/perl

use strict;

my ($inputdir,$postfix) = @ARGV[0,1];
my $Usage = "\n\t$0 <inputdir> <postfix>
\n";
die $Usage unless (@ARGV == 2);

my $header;
my @merged;
chomp(my @file = `find $inputdir -name "*.$postfix"`);
foreach (@file) {
	my $fir = 0;
	open IN,'<',$_ or die;
	while (<IN>) {
		chomp;
		if ($fir == 0) {
			$header = $_;
			$fir ++;
		}else{
			push @merged, $_;
		}
	}
	close IN;
}

print "$header\n";
foreach (@merged) {
	print "$_\n";
}
