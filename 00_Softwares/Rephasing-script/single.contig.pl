#!/usr/bin/perl

use strict;
use Cwd;
use Cwd 'abs_path';

my $dir = @ARGV[0];
my $dir2 = @ARGV[1];
my $dir3 = @ARGV[2];
my $thread = @ARGV[3];
my $Usage = "\n\t$0 <input vcf dir> <script folder> <output folder> <thread>\n";
die $Usage unless (@ARGV == 4);

$dir = abs_path($dir);

my $phasing_script = "contig_Hapi_phasing.pl";

mkdir $dir3;
my $dir3 = abs_path($dir3);

mkdir $dir2;
chdir $dir2;

my @file = glob "$dir/*.vcf";

my ($counter,@tmp,$jobs);

foreach (@file) {
        $counter ++;
        push @tmp,$_;
        if ($counter == 1) {
                my $vcf_list = join(",",@tmp);
                $jobs ++;
                $counter = 0;
                &phasing_sge($vcf_list,$jobs);
                @tmp = ();
        }
}


sub phasing_sge {
        my ($vcf_list,$job_id) = @_;
        my $sge = "

$phasing_script $vcf_list $dir3 $thread
";
        open SGE,'>',"$job_id.sh";
        print SGE "$sge";
        close SGE;
        `bash $job_id.sh`;
}
