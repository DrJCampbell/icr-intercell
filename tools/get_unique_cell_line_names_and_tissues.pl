#!/usr/bin/perl -w
# script used to get unique gene names from COSMIC cell line names and 
# tissues extracted from CosmicCLP_CompleteExport_v68.tsv.gz
# jamesc@icr.ac.uk, 13th May 2014

use strict;

open IN, "< CosmicCLP_CompleteExport_v68.cell_line_names_and_tissues.txt" or die "can't open input: $!\n";

my %names_tissues;

while(<IN>){
	my ($cell_line, $tissue) = split(/\t/);
	chomp($tissue);
	$names_tissues{$cell_line} = $tissue;
}

close IN;

open OUT, "> CosmicCLP_CompleteExport_v68.cell_line_names_and_tissues.unique.txt" or die "Can't open output: $!\n";

while(my ($key, $value) = each %names_tissues){
	print OUT "$key\t$value\n";
}

close OUT;


