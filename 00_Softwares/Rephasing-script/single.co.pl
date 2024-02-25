#!/usr/bin/perl

use strict;
use Cwd;
use Cwd 'abs_path';
use Parallel::ForkManager;

my ($dir,$cpu) = @ARGV[0,1];
my $Usage = "\n\t$0 <hapOut dir> <cpu>
\n";
die $Usage unless (@ARGV == 2);

$dir = abs_path($dir);

## ForkManager setup ##
my $pm = Parallel::ForkManager->new($cpu);

my @sca = glob "$dir/*hap.txt";

foreach (@sca) {
	my $prefix = (split /\./,(split /\//,$_)[-1])[0];
	$pm->start and next;
	&Hapi_co($_,"$dir/$prefix.co.txt");
	$pm->finish;
}
$pm->wait_all_children();

sub Hapi_co {
	my ($hap,$out) = @_;
	`Hapi_co.R $hap $out`;
}

