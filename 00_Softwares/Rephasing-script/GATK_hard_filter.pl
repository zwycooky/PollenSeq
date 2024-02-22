#!/usr/bin/perl
#
use strict;
my ($vcf,$output) = @ARGV[0,1];
my $Usage = "\n\t$0 <VCF> <OUTPUT>
\n";
die $Usage unless (@ARGV == 2);

print "START GATK\n";
my $gatk = "gatk";

print "Selecting SNPs\n";
#Select SNPs
!system "MALLOC_ARENA_MAX=2 $gatk SelectVariants -select-type SNP -V $vcf --verbosity ERROR -O $vcf.snp.vcf" or die "ERROR with GATK:$!";

print "Filter SNPs\n";
#Filter SNPs
!system "MALLOC_ARENA_MAX=2 $gatk VariantFiltration -V $vcf.snp.vcf --verbosity ERROR --filter-expression \"QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filter-name \"Filter\" -O $vcf.snp.filter.vcf" or die "ERROT with GATK:$!";

print "Select InDels\n";
#Select INDELs
!system "MALLOC_ARENA_MAX=2 $gatk SelectVariants -select-type INDEL --verbosity ERROR -V $vcf -O $vcf.indel.vcf" or die "ERROR with GATK:$!";

print "Filter InDels\n";
#Filter INDELs
!system "MALLOC_ARENA_MAX=2 $gatk VariantFiltration -V  $vcf.indel.vcf --verbosity ERROR --filter-expression \"QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -20.0\" --filter-name \"Filter\" -O $vcf.indel.filter.vcf" or die "ERROR with GATK:$!";

print "Merge SNPs and InDels\n";
#Merge VCF
!system "MALLOC_ARENA_MAX=2 $gatk MergeVcfs -I $vcf.snp.filter.vcf -I $vcf.indel.filter.vcf -O $output" or die "ERROR with GATK:$!";

print "Remove intermidiate files\n";
`rm $vcf.snp.vcf $vcf.snp.vcf.idx $vcf.snp.filter.vcf $vcf.snp.filter.vcf.idx $vcf.indel.vcf $vcf.indel.vcf.idx $vcf.indel.filter.vcf $vcf.indel.filter.vcf.idx`;

print "DONE!\n";
