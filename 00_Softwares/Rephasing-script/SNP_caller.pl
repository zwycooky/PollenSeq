#!/usr/bin/perl

use strict;
use Cwd;
use Cwd 'abs_path';
use Parallel::ForkManager;

my ($input_dir, $tmp_output, $output, $genome, $num_threads) = @ARGV[0, 1, 2, 3, 4];
my $Usage = "\n\t$0 <input bam dir> <tmp output dir> <output dir> <genome> <num threads>\n";
die $Usage unless (@ARGV == 5);

!system "export JAVA_TOOL_OPTIONS=-Xmx5g";

$input_dir = abs_path($input_dir);
$tmp_output = abs_path($tmp_output);
$output = abs_path($output);
$genome = abs_path($genome);
my $gatk = "gatk";
my $postfix = "sorted.marked.duplicates.bam";

## find bam file ##
chomp(my @bam = `find $input_dir -name "*.$postfix"`);
if (@bam == 0) {
    die "No bam file found:$!";
}

mkdir $tmp_output;
mkdir $output;

## get chromosome name and length ##
chomp(my @sca_id = `grep '>' $genome`);
open GENOME, '<', "$genome" or die;
my ($gname, %glen);
while (<GENOME>) {
    chomp;
    if (/>/) {
        $gname = (split />/, (split)[0])[1];
        $glen{$gname} = '';
    } else {
        $glen{$gname} += length($_);
    }
}
close GENOME;

## get bam file input command ##
my @bam_com;
foreach (@bam) {
    my $tmp = "-I $_";
    push @bam_com, $tmp;
}
my $bam_com = join(" ", @bam_com);

## ForkManager setup ##
my $pm = Parallel::ForkManager->new($num_threads);


chomp(my $pwd = `pwd`);
foreach (@sca_id) {
    my $id = (split />/, (split)[0])[1];
    chdir "$pwd";
    mkdir "$tmp_output/$id";
    chdir "$tmp_output/$id";
    my $glen = $glen{$id};
    my $window = 10000000;
    my $step = 9900000;
    my $num = int($glen / $step) + 1;

    foreach (1..$num) {
        my $postfix = $_;
        my $split_s = 1 + ($_-1) * $step;
        my $split_e = $split_s + $window;
        if ($split_e >= $glen) {
                $split_e = $glen;
        }

        my $lsf_script = &lsf($id, $split_s, $split_e, $postfix);

        my $lsf_filename = "$tmp_output/$id/$id.$postfix.sh";
        open LSF, '>', $lsf_filename;
        print LSF "$lsf_script";
        close LSF;

        $pm->start and next;
        system("bash $lsf_filename");
        $pm->finish;

    }
}

$pm->wait_all_children();

print("All processes have been finished\n");


sub lsf {
    my ($id, $s, $e, $postfix) = @_;
    my $lsf;

    if ($postfix != 0) {
        $lsf = <<LSF;

cd $tmp_output/$id

MALLOC_ARENA_MAX=2 gatk HaplotypeCaller -R $genome $bam_com --verbosity ERROR -L $id:$s-$e -O $id.$postfix.raw.vcf
LSF
    } else {
        $lsf = <<LSF;

cd $tmp_output/$id

MALLOC_ARENA_MAX=2 gatk HaplotypeCaller -R $genome $bam_com --verbosity ERROR -L $id -O $id.raw.vcf
LSF
    }
    return ($lsf);
}
