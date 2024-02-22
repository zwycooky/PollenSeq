#!/usr/bin/perl

use strict;
use Cwd 'abs_path';

my ($input_dir, $len_file, $outdir) = @ARGV[0, 1, 2];
my $Usage ="\n\t $0 <hapOut> <scaffold lenght file> <Out Dir>
\n";
die $Usage unless (@ARGV == 3);

$len_file = abs_path($len_file);

my @co_file = glob "$input_dir/*co.txt";
my @hap_file = glob "$input_dir/*hap.txt";

my %co_file;
foreach (@co_file) {
        my $file_path = abs_path($_);
        my $contig_id = (split /\./, (split /\//,$file_path)[-1])[0];
        $co_file{$contig_id} = $file_path;
}

my %hap_file;
foreach (@hap_file) {
        my $file_path = abs_path($_);
        my $contig_id = (split /\./, (split /\//,$file_path)[-1])[0];
        $hap_file{$contig_id} = $file_path;
}

my $find_bin_script = "find_bin.R";
my $find_bin_script2 = "get_co_from_sca.R";

if (!-e $outdir) {
        mkdir $outdir;
        chdir $outdir;
}else{
        chdir $outdir;
}

foreach (sort keys %hap_file) {
        my $contig_id = $_;
        if (exists $co_file{$_}) {
                !system "$find_bin_script $co_file{$_} $hap_file{$_}" or die "Error with find bins:$!";
                print "$find_bin_script $co_file{$_} $hap_file{$_}\n";
        }else{
                !system "$find_bin_script2 $hap_file{$_} $len_file" or die "Error with find bins2:$!";
        }
        #exit;
}
