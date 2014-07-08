#!/usr/bin/perl -w
# ===================================================== #
# *** Depricated - renamed to process mutation data ***
# ===================================================== #
# process_CCLE_mutations.pl
# read data from CCLE, count mutations by type
# for each gene and cell line
# jamesc@icr.ac.uk, 14th May 2014
# ===================================================== #

use strict;
use Getopt::Long;

my ($help, $ccle_data, $mut_freqs, $genes, $cell_lines, $output);

GetOptions (
  "ccle_data=s" => \$ccle_data,
  "mut_freqs=s" => \$mut_freqs,
  "genes=s" => \$genes,
  "cell_lines=s" => \$cell_lines,
  "output=s" => \$output,
  "help" => \$help,
 );

# print usage message if requested or no args supplied
if(defined($help)) {
  &usage;
  exit(0);
}

$output = $ccle_data . ".table" unless defined $output;


# read cell line name dictionary
open CL, "< $cell_lines" or die "Can't read cell line name dictionary $cell_lines: $!\n";

my %cell_lines;

while(<CL>){
  next if /^#/;
  my ($cell_line, $standard_cell_line, $source) = split(/\t/);
  $cell_lines{$cell_line} = $standard_cell_line;
}

close CL;

# read gene name dictionary
open GN, "< $genes" or die "Can't read gene name dictionary $genes: $!\n";

my %genes;

while(<GN>){
  next if /^#/;
  my ($gene_alias, $standard_gene_name) = split(/\t/);
  chomp($standard_gene_name);
  $genes{$gene_alias} = $standard_gene_name;
}

close GN;

# read in the mutation frequency data
open MFRQ, "< $mut_freqs" or die "Can't read mutation frequency information file $mut_freqs: $!\n";

my %mut_freqs;
while(<MFRQ>){
  my ($gene_name, $prot_pos, $var_type, $freq) = split(/\t/);
  chomp($freq);
  my $key = "$gene_name\t$prot_pos";
  $mut_freqs{$key} = $freq;
}

close MFRQ;

# types of mutation consequences present in file:
my %mutation_consequences = (
	"De_novo_Start_InFrame" => "aa_sub",
	"In_Frame_Del" => "aa_sub",
	"In_Frame_Ins" => "aa_sub",
	"Missense_Mutation" => "aa_sub",
	"Nonstop_Mutation" => "aa_sub",
	"Stop_Codon_DNP" => "aa_sub",
	"Stop_Codon_Ins" => "aa_sub",
	"3'UTR" => "none",
	"5'Flank" => "none",
	"5'UTR" => "none",
	"Intron" => "none",
	"Silent" => "none",
	"De_novo_Start_OutOfFrame" => "trunc",
	"Frame_Shift_Del" => "trunc",
	"Frame_Shift_Ins" => "trunc",
	"Nonsense_Mutation" => "trunc",
	"Splice_Site_Del" => "trunc",
	"Splice_Site_DNP" => "trunc",
	"Splice_Site_Ins" => "trunc",
	"Splice_Site_SNP" => "trunc",
	"Start_Codon_Del" => "trunc"
	);
	

# ====================================================== #
# Want a table with one row per cell line and one column
# per (gene*consequence). The consequences could be like
# truncating (fs,*,ess) / recurrent missense / other
# 
# How can we represent more complex data such as position
# of mutations in a protein and z-scores? As a binned 
# image, scaled to protein length.
# ====================================================== #

my %ccle_truncs;		# store truncating counts for cell*gene
my %ccle_rec_mis;		# store recurrent missense counts for cell*gene
my %ccle_other;			# store other counts for cell*gene
my %cell_lines_seen;	# store standard cell line names for all seen
my %genes_seen;			# store standard gene names for all seen
my %mutations_seen;		# store cell line * genome mutation to prevent counting repeats

# Open and read mutations file
open MUTS, "< $ccle_data" or die "Can't read mutations file $ccle_data: $!\n";

my $header = <MUTS>;

while(<MUTS>){

  # get gene, cell line, prot. mut., consequence etc.
  # [1]		Entrez_Gene_Id
  # [8]		Variant_Classification	// Missense_Mutation, Intron, Frame_Shift_Del
  # [9]		Variant_Type			// DEL, SNP
  # [15]	Tumor_Sample_Barcode	// JHUEM2_ENDOMETRIUM
  # [32]	Genome_Change
  # [37]	Protein_Change

  my @fields = split(/\t/);
  my $entrez_gene = $fields[0];
  my $var_class = $fields[8];
  my $var_type = $fields[9];
  my $sample = $fields[15];
  my $genome_change = $fields[32];
  my $prot_change = $fields[37];
  
  # lookup gene and cell line name in dictionaries
  my $standard_gene = $entrez_gene;
  if(exists($genes{$entrez_gene})){
    $standard_gene = $genes{$entrez_gene};
  }
  my $standard_cell_line = $sample;
  if(exists($cell_lines{$sample})){
    $standard_cell_line = $cell_lines{$sample};
  }

  # store a hash key for every cell line * genome_change * consequence seen
  # skip if already counted...
  my $sample_genome_change_key = $standard_cell_line . "_" . $genome_change . "_" . $var_class;
  if(exists($mutations_seen{$sample_genome_change_key})){
    print "already seen $sample_genome_change_key - skipping.\n";
    next;
  }
  else{
    $mutations_seen{$sample_genome_change_key} = 1;
  }
  
  my $ccle_key = "$standard_cell_line\t$standard_gene";
  
  $cell_lines_seen{$standard_cell_line} = 1;
  $genes_seen{$standard_gene} = 1;
  
  if($mutation_consequences{$var_class} eq "trunc"){	# if mutation is truncating
    $ccle_truncs{$ccle_key} ++;							# increment cellline+gene truncating count
  }
  elsif($mutation_consequences{$var_class} eq "aa_sub"){
    # check if the mutation is recurrent as per Davoli
    my $davoli_key = "$standard_gene\t$prot_change";
    if(exists($mut_freqs{$davoli_key})){
      $ccle_rec_mis{$ccle_key} ++;
    }
    else{
      $ccle_other{$ccle_key} ++;
    }
  }
} # finished reading file, close
close MUTS;


my @cell_lines_seen = keys %cell_lines_seen;
my @genes_seen = keys %genes_seen;

my $output_data = '';

my $output_header = "cell_line";
foreach my $seen_gene (@genes_seen){
  $output_header .= "\t$seen_gene (trunc)\t$seen_gene (rec_mis)\t$seen_gene (other)";
}
$output_header .= "\n";

foreach my $seen_cell_line (@cell_lines_seen){
  $output_data .= $seen_cell_line;
  foreach my $seen_gene (@genes_seen){
    my $ccle_key = "$seen_cell_line\t$seen_gene";
    if(exists($ccle_truncs{$ccle_key})){
      $output_data .= "\t$ccle_truncs{$ccle_key}";
    }
    else{
      $output_data .= "\t0";
    }
    if(exists($ccle_rec_mis{$ccle_key})){
      $output_data .= "\t$ccle_rec_mis{$ccle_key}";
    }
    else{
      $output_data .= "\t0";
    }
    if(exists($ccle_other{$ccle_key})){
      $output_data .= "\t$ccle_other{$ccle_key}";
    }
    else{
      $output_data .= "\t0";
    }
  }
  $output_data .= "\n";
}
open OUT, "> $output" or die "Can't write to output file $output: $!\n";
print OUT "$output_header$output_data";
close OUT;


sub usage() {
  my $usage =<<END;
#-----------------------------#
#  process_CCLE_mutations.pl  #
#-----------------------------#

James Campbell (jamesc\@icr.ac.uk)

Usage:
perl process_CCLE_mutations.pl [options]

Options
--help			Display this message and quit
--ccle_data		Path to the CCLE maf file [required]
--mut_freqs		Path to the processed Davoli data [required]
--genes			Path to the gene name dictionary [required]
--cell_lines		Path to the cell line name dictionary [required]
--output		Path to output. Defaults to ccle_data.proc [optional]
END

  print $usage;
}







