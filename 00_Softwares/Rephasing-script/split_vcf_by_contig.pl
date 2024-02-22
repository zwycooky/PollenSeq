#!/usr/bin/perl
#
use strict;
use Cwd 'abs_path';

my ($input,$output) = @ARGV[0,1];
my $Usage = "\n\t$0 <input vcf> <output dir>
\n";
die $Usage unless (@ARGV == 2);

mkdir $output;
$output = abs_path($output);
chomp(my $pwd = `pwd`);
#$output = "$pwd/$output";

open IN,'<',"$input" or die;
my ($sample,%chr,$counter,@vcf_dir,%pos_ex);
while (<IN>) {
	chomp;
	if (/\A##/) {
        	next;
    	}elsif (/\A#CHROM/) {
        	my @sample = split;
        	@sample = @sample[9..(@sample-1)];
        	$sample = join("\t",@sample);
        	$sample = "chr\tpos\tref\talt\t$sample";
	}else{
		my ($chr,$pos,$ref,$alt) = (split)[0,1,3,4];
    #    if ($chr =~ /sca/) { next };
		
        	my @snp = split;
        	@snp = @snp[9..(@snp-1)];
        	my @snp2;
        	foreach (@snp) {
            		if ($_ eq './.') {
                		push @snp2, "NA";
            		}elsif ($_ eq '0/0') {
                		push @snp2, $ref;
            		}elsif ($_ eq '1/1') {
                		push @snp2, $alt;
            		}else{
                		die ":$!";
            		}
        	}
        	my $snp = join("\t",@snp2);
        	@snp = (); @snp2 = ();

		if ($counter == 0) {
			open OUT,'>',"$output/$chr.vcf";
            		print OUT "$sample\n";
            		print OUT "$pos\t$chr\t$pos\t$ref\t$alt\t$snp\n";
			$chr{$chr} = 1;
            		$counter = 1;
			$pos_ex{$pos} = 1;
			next;
		}elsif (!exists $chr{$chr}) {
			close OUT;
			%pos_ex = ();

            		open OUT,'>',"$output/$chr.vcf";
			print OUT "$sample\n";
            		print OUT "$pos\t$chr\t$pos\t$ref\t$alt\t$snp\n";
			$chr{$chr} = 1;
			$pos_ex{$pos} = 1;			

		}else{
			if (!exists $pos_ex{$pos}) {
				print OUT "$pos\t$chr\t$pos\t$ref\t$alt\t$snp\n";
				$pos_ex{$pos} = 1;
			}
		}
	}
}
close IN;
close OUT;





