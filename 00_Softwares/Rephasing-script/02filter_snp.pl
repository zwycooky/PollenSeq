#!/usr/bin/perl

use strict;

my ($vcf,$exclude) = @ARGV[0,1];
my $Usage = "\n\t$0 <vcf> <exclude samples>
\tfilter variants with missing rate greater than 0.8, heterozygous and non-bialleic
\n";
die $Usage unless (@ARGV == 2);

my %ex;
open EX,'<',"$exclude" or die;
while (<EX>) {
        chomp;
        $ex{$_} = 1;
}
close EX;

open OUT1,'>',"02Filter.report.txt";
open OUT2,'>',"02HetSNPs_matrix.txt";
print OUT1 "chr\tposition\tRef_allele\tAlt_allele\tHeterozygoes\tMutiallelic\tTotal_SNPs\n";
open VCF,'<',"$vcf" or die "Cannot open vcf file:$!";
my @ex_sam;
while (<VCF>) {
        chomp;
        if (/\A##/) {
                print "$_\n";
        }elsif (/\A#CHROM/) {
                my @line = split;
                my $ex_count = 0;
                foreach (@line) {
                        if (!exists $ex{$_}) {
                                push @ex_sam,$ex_count;
                        }
                        $ex_count ++;
                }
                @line = @line[@ex_sam];
                my $line = join("\t",@line);
                print "$line\n";
        }else{
                my @line = split;
                @line = @line[@ex_sam];
                my ($chr,$pos,$refbase,$altbase) = (split)[0,1,3,4];

                if (length($refbase) > 1 || length($altbase) > 1) { next };

                $line[7] = '.';
                $line[8] = 'GT';
                my ($het,$mut,$ref,$alt,$miss) = (0,0,0,0,0);
                my @het_line;
                foreach (9..(@line-1)) {
                        my $gt = (split /:/,$line[$_])[0];
                        if ($gt eq '0/1') {
                                $line[$_] = './.';
                                $het ++;
                                push @het_line,1;
                        }elsif ($gt eq '1/1') {
                                $line[$_] = '1/1';
                                $alt ++;
                                push @het_line,0;
                        }elsif ($gt eq '0/0') {
                                $line[$_] = '0/0';
                                $ref ++;
                                push @het_line,0;
                        }elsif ($gt eq '1/2' || $gt =~ '2/2') {
                                $line[$_] = './.';
                                $mut ++;
                                push @het_line,0;
                        }else{
                                $line[$_] = './.';
                                $miss ++;
                                push @het_line,0;
                        }
                }

                my $het_line = join ("\t",@het_line);
                if ($het > 1) {
                        #my $chr_numbers = (split /chr/,$chr)[1];
                        print OUT2 "$chr\t$pos\t$het_line\n";
                }

                my $sum = $het + $alt + $ref + $mut + $miss;

                if ($ref == 0 || $alt == 0) {
                        print OUT1 "$chr\t$pos\t$ref\t$alt\t$het\t$mut\t$sum\n";
                        next;
                }

                my $ratio = $ref / ($ref + $alt);

                if ($miss / $sum > 0.8 || $mut > 0 || $het > 1 || $ratio > 0.75 || $ratio < 0.25) {
                        print OUT1 "$chr\t$pos\t$ref\t$alt\t$het\t$mut\t$sum\n";
                        next;
                }

                my $line = join("\t",@line);
                print "$line\n";

        }
}
close OUT1;
close OUT2;
