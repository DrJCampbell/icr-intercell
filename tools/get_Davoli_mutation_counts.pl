#!/usr/bin/perl -w
# extract and count gene*mutations from
# the mutation dataset published by
# Davoli et al. (2012) and available
# from 
# http://www.cell.com/cell/supplemental/S0092-8674%2813%2901287-7

use strict;

my $gene_name_dictionary_file = "../resources/gene_name_dictionary.txt";
my $mutations_file = "../resources/Davoli_mutation_data/Mutation_Dataset.txt";


# format of the data is:
#Gene    Genome.position.hg19    Reference       Mutation        Protein_Change  Mutation_Type   Tumor_Sample    Tumor_Type
#A1BG    19:58861774-58861774    C       -       p.G385fs        Indel Frameshift        TCGA-18-3406-01A        Lung Squamous Cell Carcinoma
#A1BG    19:58864380-58864381    CC      -       p.G85fs Indel Frameshift        TCGA-G2-A3IE-01A        Bladder Carcinoma
#A1BG    19:58862784-58862784    C       T       p.A295T Missense        pfg019T Stomach Adenocarcinoma

# Note that the Gene column contains official HGNC symbols (eg KMT2B, not MLL2 or MLL4...)


# read in the gene name dictionary
my %gene_names;
open GN, "< $gene_name_dictionary_file" or die "unable to open gene name dictionary $gene_name_dictionary_file: $!\n";
while(<GN>){
  next if /^#/;
  my ($alias, $gene_name) = split(/\t/);
  chomp($gene_name);
  $gene_names{$alias} = $gene_name;
}
close GN;


open MUTS, "< $mutations_file" or die "Can't read mutations input file $mutations_file: $!\n";
my %gene_pos_mut_counts;
while(<MUTS>){

  my @fields = split(/\t/);
  my $gene = $fields[0];
  my $prot_pos = $fields[4];
  my $effect = $fields[5];

  my $standard_gene = '';
  if(exists($gene_names{$gene})){
  	$standard_gene = $gene_names{$gene};
  }
  else{
  	warn "$gene not found in the gene name dictionary - using the original name.\n";
  	$standard_gene = $gene;
  }
  
  $gene_pos_mut_counts{$standard_gene . "\t" . $prot_pos . "\t" . $effect} ++;
  
}
close MUTS;

my $outfile = $mutations_file . '.mut_counts.txt';

open PROT_POS_COUNTS, "> $outfile" or die "Can't write protein position counts: $!\n";
while(my($prot_pos, $count) = each %gene_pos_mut_counts){
  next if $count < 2;
  print PROT_POS_COUNTS "$prot_pos\t$count\n";
}
close PROT_POS_COUNTS;


