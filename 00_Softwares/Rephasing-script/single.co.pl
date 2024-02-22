#!/usr/bin/perl

use strict;
use Cwd;
use Cwd 'abs_path';

my $dir = @ARGV[0];
my $Usage = "\n\t$0 <hapOut dir>\n";
die $Usage unless (@ARGV == 1);

$dir = abs_path($dir);

mkdir "co_scripts";
chdir "co_scripts";

my @sca = glob "$dir/*hap.txt";

foreach (@sca) {
my $prefix = (split /\./,(split /\//,$_)[-1])[0];
my $sge = "

Hapi_co.R $_ $dir/$prefix.co.txt

";
    open SGE,'>',"$prefix.co.sh";
    print SGE "$sge\n";
    close SGE;
   `bash $prefix.co.sh`;
   print "$prefix.co completed!\n";

}

