#!/usr/bin/perl

use strict;

my $hap_dir = "/public/home/hjqiu/people/zwy/Fuding_population/04phasing_FD/03phasing/hapOut";
my @hap_file = glob "$hap_dir/*.hap.txt";

my @haps;
foreach (@hap_file) {
	my $fir = 0;
	open IN,'<',"$_" or die;
	while (<IN>) {
		chomp;
		if ($fir == 0) {
			$fir = 1;
			next;
		}
		my ($chr,$pos,$hap1,$hap2,$ratio,$confi) = (split)[1,2,-5,-4,-2,-1];
		if ($ratio > 0.8 && $confi ne 'L') {
			push @haps,"$chr\t$pos\t$hap1\t$hap2";
		}
	}
	close IN;
}

foreach (@haps) {
	print "$_\n";
}

