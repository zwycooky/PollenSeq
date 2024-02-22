#!/usr/bin/perl

use strict;

chomp(my @list = <>);

@list = sort {(split /\./,(split /\//,$a)[-1])[0] cmp (split /\./,(split /\//,$b)[-1])[0] or (split /\./,(split /\//,$a)[-1])[1] <=> (split /\./,(split /\//,$b)[-1])[1] } @list;

foreach (@list) {
	print "$_\n";
}
