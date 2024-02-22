#!/usr/bin/perl

use strict;
my ($pvcf,$svcf) = @ARGV[0,1];
my $Usage = "\n\t$0 <parent vcf> <single sperm vcf>
\tOnly heterozygotes in 3W-4 are keeped\n";
die $Usage unless (@ARGV == 2);

my %her;
open PVCF,'<',"$pvcf" or die;
while(<PVCF>) {
        chomp;
        if (/\A#/) {
                next;
        }else{
                my ($chr,$pos) = (split)[0,1];
                my $geno = (split /:/,(split)[9])[0];
                if ($geno eq '0/1') {
                        my $key = "$chr\t$pos";
                        $her{$key} = 1;
                }
        }
}
close PVCF;

open IN,'<',"$svcf" or die "Cannot open Single sperm vcf file:$!";
while (<IN>) {
        chomp;
        if (/\A#/) {
                print "$_\n";
        }else{
                my $line = $_;
                my ($chr,$pos) = (split)[0,1];
                my $key = "$chr\t$pos";
                if (exists $her{$key}) {
                        print "$line\n";
                }
        }
}
close IN;
