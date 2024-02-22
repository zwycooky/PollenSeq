#!/usr/bin/perl

use strict;

my ($bins,$binpos,$wrong,$haps_dir,$linkage_group_dir, $output) = @ARGV[0,1,2,3,4,5];
my $Usage = "\n\t$0 <re-pased bins> <binmap.txt> <wrong phasing bins> <haps dir> <linkage group dir> <output>
\n";
die $Usage unless (@ARGV == 6);

my %keep_bins;
open IN,'<',"$bins" or die;
while (<IN>) {
	chomp;
	if (/\Alocus_name/) { next };
	my $id = (split)[0];
	#if ($id eq $exclude) { next };
	$keep_bins{$id} = 1;
}
close IN;

my $binmap;
open IN,'<',"$binpos" or die;
while (<IN>) {
	chomp;
	if (/\Achr/) { next };
	my ($id,$chr,$s,$e) = split;
	if (exists $keep_bins{$id}) {
		push @{$binmap->{$chr}}, "$id\t$s\t$e";
	}
}
close IN;

my %wrong_bin;
open IN,'<',"$wrong" or die;
while (<IN>) {
	chomp;
	$wrong_bin{$_} = 1;
}
close IN;

## @{$binmap->{$chr}}, %wrong_bin ##
my @hapfile = glob "$haps_dir/*.hap.txt";
my $correct_hap;
foreach (@hapfile) {
	print "$_ start\n";
	open HAP,'<',"$_" or die;
	while (<HAP>) {
		chomp;
		if (/\Achr/) { next };
		my ($chr,$pos,$hap1,$hap2,$rate,$frame) = (split)[1,2,-5,-4,-2,-1];
		if ($rate <= 0.8 || $frame eq 'L') { next };
  		if (!exists $binmap->{$chr}) { next };
		my @tmp = @{$binmap->{$chr}};
		foreach (@tmp) {
			my ($id,$s,$e) = split;
			if ($pos >= $s && $pos <= $e) {
				
				my $chr_num = (split /chr/,$chr)[1];
				my $snp_id = "$chr_num\_$pos";
				
				if (exists $wrong_bin{$id}) {
					push @{$correct_hap->{$id}}, "$snp_id\t$hap2\t$hap1";
				}else{
					push @{$correct_hap->{$id}}, "$snp_id\t$hap1\t$hap2";
				}
				
				last;
			}
		}
	}
	close HAP;
}

open OUT,'>',$output;
my @linkage_file = glob "$linkage_group_dir/*.txt";
foreach (@linkage_file) {
	my $chr = (split /\./,(split /\//,$_)[-1])[0];
	open LINK,'<',"$_" or die;
	
	while (<LINK>) {
		chomp;
		my ($id,$cM) = split;
		my @tmp = @{$correct_hap->{$id}};
		foreach (@tmp) {
			print OUT "$id\t$chr\t$_\n";
		}
	}
	close LINK;
	
}
close OUT;
