#!/usr/bin/perl

use strict;
my $contig_len = "/public/home/hjqiu/people/zwy/Fuding_population/04phasing_FD/03phasing/find_bins/FD_hifi.asm.bp.p_ctg.length.txt";
my $inputdir = @ARGV[0];
my $Usage = "\n\t$0 <inputdir>
\n";
die $Usage unless (@ARGV == 1);

my ($total_len,%contig_len);
open IN,'<',"$contig_len" or die;
while (<IN>) {
	chomp;
	my ($chr,$len) = (split)[0,1];
	$total_len += $len;
	$contig_len{$chr} = $len;
}
close IN;

my $contig_phasing_len;
chomp(my @file = `find ./ -name "*.binmap.txt"`);
foreach (@file) {
	my $name = (split /\./,(split /\//,$_)[-1])[0];
	$contig_phasing_len += $contig_len{$name};
}

my $phasing_per = sprintf ("%.2f", ($contig_phasing_len / $total_len) * 100);
print "$phasing_per%\n";

