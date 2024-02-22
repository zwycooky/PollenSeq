#!/usr/bin/perl

use strict;

my ($vcf_list, $dir, $thread) = @ARGV[0, 1, 2];
my $Usage = "\n\t$0 <vcf list> <output dir> <thread>
\n";
die $Usage unless (@ARGV == 3);

my @file = (split /,/,$vcf_list);


foreach (@file) {
        my $contig = (split /\./,(split /\//,$_)[-1])[0];
        !system "Hapi_step_by_step_parallel.R $_ $dir/$contig.hap.txt $thread" or warn "warnning with hapi phasing $contig:$!";
        print "Hapi_step_by_step_parallel.R $_ $dir/$contig.hap.txt $thread\n";
        print "$contig DONE!\n";
}
