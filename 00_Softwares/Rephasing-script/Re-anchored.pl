#!/usr/bin/perl

use strict;

my ($linkage_dir,$binmap) = (@ARGV)[0,1];
my $Usage = "\n\t$0 <linkage group dir> <merged_binmap>
\n";
die $Usage unless (@ARGV == 2);

my %bins;
open BIN,'<',"$binmap" or die;
while (<BIN>) {
	chomp;
	if (/\Achr/) { next };
	my ($mid,$chr,$s,$e) = split;
	$bins{$mid} = "$chr\t$s\t$e";	
}
close BIN;

## merge continues bins ##
my @file = glob "$linkage_dir/*.txt";
foreach (@file) {
	my $chr_prefix = (split /\./,(split /\//,$_)[-1])[0];
	open IN,'<',"$_" or die;
	my ($merged,$fir,$pre_id,$pre_chr);
	while (<IN>) {
		chomp;
		my $bin_id = (split)[0];
		my ($chr,$numid) = (split /_/,$bin_id)[1,2];
		if ($fir == 0) {
			$fir ++;
			$pre_id = $numid;
			$pre_chr = $chr;
			push @{$merged->{$fir}},$bin_id;
		}elsif ( abs($pre_id - $numid) == 1 && $pre_chr eq $chr) {
			$pre_id = $numid;
			$pre_chr = $chr;
			push @{$merged->{$fir}},$bin_id;
		}else{
			$fir ++;
			$pre_id = $numid;
			$pre_chr = $chr;
			push @{$merged->{$fir}},$bin_id;
		}
	}
	close IN;
	
	my $pre_end;
	my @assembly;
	
	foreach (1..$fir) {
		if ($_ == 1) {
			$pre_end = 0;
		}
		my @tmp = @{$merged->{$_}};
		if (@tmp == 1) {
			my $bin_id = $tmp[0];
			my ($chr,$s,$e) = (split /\t/,$bins{$bin_id});
			my $len = $e - $s + 1;
			my $format_s = $pre_end;
			my $format_e = $pre_end + $len;
			push @assembly,"$bin_id\t$format_s\t$format_e\t+";
			$format_e += 500;
			$pre_end = $format_e;
			next;
		}
		
		my $num1 = (split /_/,$tmp[0])[-1];
		my $num2 = (split /_/,$tmp[-1])[-1];
		
		if ($num1 < $num2) {
			
			my $fir1 = 0;
			my $pree = 0;
			foreach (@tmp) {
				my ($chr,$s,$e) = (split /\t/,$bins{$_});
				my $len = $e - $s + 1;
				if ($fir1 == 0) {
					my $format_s = $pre_end;
					my $format_e = $pre_end + $len;
					push @assembly,"$_\t$format_s\t$format_e\t+";
					$pre_end = $format_e;
					$pree = $e;
					$fir1 = 1;
				}else{
					my $format_s = $pre_end + ($s - $pree);
					my $format_e = $format_s -1 + $len;
					push @assembly,"$_\t$format_s\t$format_e\t+";
					$pre_end = $format_e;
					$pree = $e;
				}
			}
			
			$pre_end += 500;
			
		}else{
		
			my $fir1 = 0;
			my $pres = 0;
			foreach (@tmp) {
				my ($chr,$s,$e) = (split /\t/,$bins{$_});
				my $len = $e - $s + 1;
				if ($fir1 == 0) {
					my $format_s = $pre_end;
					my $format_e = $pre_end + $len;
					push @assembly,"$_\t$format_s\t$format_e\t-";
					$pre_end = $format_e;
					$pres = $s;
					$fir1 = 1;
				}else{
					my $format_s = $pre_end + ($pres - $e);
					my $format_e = $format_s -1 + $len;
					push @assembly,"$_\t$format_s\t$format_e\t-";
					$pre_end = $format_e;
					$pres = $s;
				}
			}
			
			$pre_end += 500;
			
		}
	}
	
	foreach (@assembly) {
			print "$chr_prefix\t$_\n";
	}
	#exit;
	
}


sub overlap_contig {
	my $s1 = shift @_;
	my $e1 = shift @_;
	my @contig = @_;
	my @res;
	foreach (@contig) {
		my $line = $_;
		my ($s2,$e2) = (split)[1,2];
		my $over = &overlap($s1,$e1,$s2,$e2);
		if ($over == 1) {
			push @res, $line;
		}
	}
	return @res;
}

sub overlap {
	my ($s1,$e1,$s2,$e2) = @_;
	my $sub1 = abs($s1-$e1);
	my $sub2 = abs($s2-$e2);
	my @tmp = sort {$a <=> $b} @_;
	my $sub_all = abs($tmp[0]-$tmp[-1]);
	if ($sub_all < ($sub1+$sub2)) {
		return(1)
	}else{
		return(0)
	}
}
