#!/usr/bin/perl

use strict;
use Getopt::Std;
our($opt_i);
getopts('i:');

my $input = $opt_i;

my $Usage ="\n$0 -i <INPUT>
	\n";
die $Usage unless ($opt_i);
#############################
open IN,'<',"$input" or die "Cannot open input file:$!";
my ($name,%seq);
while (<IN>) {
	chomp;
	if (/>/) {
		$name = (split />/,(split)[0])[1];
		$seq{$name} = '';
	}else{
		$seq{$name} .= $_;
	}
}

my @seq;
foreach ( keys %seq ) {
	my $len = length($seq{$_});
	#print "$_\t$len\n";
	push @seq,"$_\t$len";
}

foreach (sort {(split /\t/,$b)[1] <=> (split /\t/,$a)[1]} @seq) {
	print "$_\n";
}
