#!/usr/bin/perl

use strict;
use Cwd;
use Cwd 'abs_path';
use Parallel::ForkManager;

my $dir = @ARGV[0];
my $dir2 = @ARGV[1];
my $dir3 = @ARGV[2];
my $hapi_thread = @ARGV[3];
my $thread = @ARGV[4];
my $Usage = "\n\t$0 <input vcf dir> <script folder> <output folder> <hapi thread> <thread>\n";
die $Usage unless (@ARGV == 5);

$dir = abs_path($dir);

my $phasing_script = "contig_Hapi_phasing.pl";

mkdir $dir3;
my $dir3 = abs_path($dir3);

mkdir $dir2;
chdir $dir2;

## ForkManager setup ##
my $pm = Parallel::ForkManager->new($thread);

my @file = glob "$dir/*.vcf";

my ($counter,@tmp,$jobs);

foreach (@file) {
        $counter ++;
        push @tmp,$_;
        if ($counter == 1) {
                my $vcf_list = join(",",@tmp);
                $jobs ++;
                $counter = 0;
				@tmp = ();
				$pm->start and next;
				&phasing_sge($vcf_list,$jobs);
				$pm->finish;
				@tmp = ();
        }
}
$pm->wait_all_children();


sub phasing_sge {
        my ($vcf_list,$job_id) = @_;
        my $sge = "

$phasing_script $vcf_list $dir3 $hapi_thread
";
        open SGE,'>',"$job_id.sh";
        print SGE "$sge";
        close SGE;
        `bash $job_id.sh`;
}
