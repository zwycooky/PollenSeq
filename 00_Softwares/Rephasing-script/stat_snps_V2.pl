#!/usr/bin/perl

use strict;

my $Usage = "\n\t$0 <vcf> <Row_snp_stat.txt> <Sample_snp_stat.txt>\n";
die $Usage unless (@ARGV == 3);

my $vcf = $ARGV[0];
my $row_snp_stat_file = $ARGV[1];
my $sample_snp_stat_file = $ARGV[2];

my @sample;
my $col = 0;
my $row = 0;

open OUT1,'>', $row_snp_stat_file or die "Cannot open Row_snp_stat.txt: $!";
print OUT1 "miss\thet\tallele_frequency\tratio_others\tsum\n";

open VCF,'<', $vcf or die "Cannot open vcf file: $!";

while (<VCF>) {
    chomp;
    if (/\A##/) {next};
    if (/#CHROM/) {
        my @line = split;
        @line = @line[9..(@line-1)];
        $col = 0;
        foreach (@line) {
            $sample[$col][0] = $_;
            $sample[$col][1] = 0; ## miss ##
            $sample[$col][2] = 0; ## het ##
            $sample[$col][3] = 0; ## ATL allel ##
            $sample[$col][4] = 0; ## others ##
            $col ++;
        }
    } else {
        $row ++;
        my @line = split;
        @line = @line[9..(@line-1)];
        $col = 0;
        foreach (@line) {
            if (/\A\.\/\./) {
                $sample[$col][1] ++;
            } elsif (/\A0\/1/) {
                $sample[$col][2] ++;
            } elsif (/\A1\/1/) {
                $sample[$col][3] ++;
            } elsif (/\A0\/0/) {
            } else {
                $sample[$col][4] ++;
            }
            $col ++;
        }
        my $stat_row = &stat_snp(@line);
        print OUT1 "$stat_row\n";
    }
}
close OUT1;

open OUT2,'>', $sample_snp_stat_file or die "Cannot open Sample_snp_stat.txt: $!";
print OUT2 "miss\thet\tallele_frequency\tratio_others\tsum\n";

foreach (1..$col) {
    my $sample_name = $sample[$_-1][0];
    my $miss = $sample[$_-1][1];
    my $het = $sample[$_-1][2];
    my $var = $sample[$_-1][3];
    my $others = $sample[$_-1][4];

    if ($row == $miss) {
        print OUT2 "$sample_name\t1.0000\tNA\tNA\tNA\t$row\n";
    } else {
        my $ra_miss = sprintf("%.4f", ($miss/$row) );
        my $ra_het = sprintf("%.4f", ($het/ ($row-$miss)) );
        my $ra_allele = sprintf("%.4f",( ($het + $var*2) / (($row-$miss)*2) ));
        my $ra_others = sprintf("%.4f", ($others/ ($row-$miss)) );
        if ($ra_allele > 0.5) {
            $ra_allele = 1 - $ra_allele;
        }
        print OUT2 "$sample_name\t$ra_miss\t$ra_het\t$ra_allele\t$ra_others\t$row\n";
    }
}
close OUT2;

## -- sub -- ##
sub stat_snp {
    my $miss = 0;
    my $gen = 0;
    my $var = 0;
    my $her = 0;
    my $others = 0;
    my @snp = @_;
    my $n = @snp;
    foreach (@snp) {
        if (/\A\.\/\./) {
            $miss ++;
        } elsif (/\A0\/0/) {
            $gen ++;
        } elsif (/\A1\/1/) {
            $var ++;
        } elsif (/\A0\/1/) {
            $her ++;
        } else {
            $others ++;
        }
    }
    my $ra_miss = sprintf("%.4f", ($miss/$n) );
    my $ra_her = sprintf("%.4f", ($her/($n-$miss)) );
    my $ra_allele = sprintf("%.4f",( ($her + $gen*2) / (($n-$miss)*2) ));
    my $ra_others = sprintf("%.4f", ($others/($n-$miss)) );
    if ($ra_allele > 0.5) {
        $ra_allele = 1 - $ra_allele
    }
    return("$ra_miss\t$ra_her\t$ra_allele\t$ra_others\t$n");
}
