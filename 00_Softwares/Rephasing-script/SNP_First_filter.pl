#!/usr/bin/perl

use strict;
use Cwd;
use Cwd 'abs_path';
use Parallel::ForkManager;

my ($raw_vcf_dir, $output, $num_threads) = @ARGV[0, 1, 2];
my $Usage = "\n\t$0 <raw vcf dir> <output> <num threads>\n";
die $Usage unless (@ARGV == 3);

!system "export JAVA_TOOL_OPTIONS=-Xmx5g";

$raw_vcf_dir = abs_path($raw_vcf_dir);
$output = abs_path($output);
my $gatk_filter = "GATK_hard_filter.pl";
my $get_pass = "get_PASS_SNPs.pl";

chomp(my @raw_vcf = `find $raw_vcf_dir -name "*.raw.vcf"`);

if (! -e "filter_scripts") {
        mkdir "filter_scripts";
}

my $dir = abs_path("filter_scripts");
chdir $dir;

## ForkManager setup ##
my $pm = Parallel::ForkManager->new($num_threads);

foreach my $job (1..@raw_vcf) {
    $pm->start and next;

    my $vcf_file = abs_path($raw_vcf[$job - 1]);
    my $file_name = (split /\.raw\.vcf/,$vcf_file)[0];
    my $out_name = $file_name . ".filtered.vcf";
    my $out_PASS = $file_name . ".filtered.PASS.vcf";

    open SGE,'>',"$job.sh";
    my $sh = "
cd $dir

$gatk_filter $vcf_file $out_name
$get_pass $out_name > $out_PASS
";

    open my $SGE, '>', "$file_name.sh";
    print $SGE $sh;
    close $SGE;

    system("bash $file_name.sh");

    $pm->finish;
}

$pm->wait_all_children();
