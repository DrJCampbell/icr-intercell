#!/usr/bin/perl -w
# ========================================= #
# prep_expression_data.pl
# select a subset of genes from the
# Sanger expression data and output with
# standardised cell line and gene names
# cell lines should be sorted by tissue
#
# jamesc@icr.ac.uk, 29th Sept 2014
# ========================================= #

use strict;
use Getopt::Long;

my ($help, $exprn_data, $genes, $cell_lines, $output, $output_genes_list);

GetOptions (
  "exprn_data=s" => \$exprn_data,
  "genes=s" => \$genes,
  "cell_lines=s" => \$cell_lines,
  "output=s" => \$output,
  "output_genes_list=s" => \$output_genes_list,
  "help" => \$help,
 );


# print usage message if requested
if(defined($help)) {
  &usage;
  exit(0);
}

$output = $exprn_data . ".subset.txt" unless defined $output;

# ============================== #
# read cell line name dictionary #
# ============================== #

open CL, "< $cell_lines" or die "Can't read cell line name dictionary $cell_lines: $!\n";

my %cell_lines;

while(<CL>){
  next if /^#/;
  my ($cell_line, $standard_cell_line, $source) = split(/\t/);
  $cell_lines{$cell_line} = $standard_cell_line;
}

close CL;

# ========================= #
# read gene name dictionary #
# ========================= #

open GN, "< $genes" or die "Can't read gene name dictionary $genes: $!\n";

my %genes;

while(<GN>){
  next if /^#/;
  my ($gene_alias, $standard_gene_name) = split(/\t/);
  chomp($standard_gene_name);
  $genes{$gene_alias} = $standard_gene_name;
}

close GN;


# ==================================== #
# read in the list of genes we want to #
# output at the end and standardise    #
# ==================================== #

# we now have upto three identifiers for each gene
# 	the gene symbol (e.g. EIF1AX)
# 	the EntrezGene ID (e.g. 101060318)
# 	the ensembl gene ID(s) (e.g. ENSG00000173674_ENSG00000198692)
# this is followed by the cancer gene classification
# (TSG or OG) at the end.
# These need to be made into hashes so that any of the ID types
# can be used to look up the full set of identifiers...

open OUTGENES, "< $output_genes_list" or die "Can't read list of genes to output from file $output_genes_list: $!\n";
my %output_genes;
my %symbol_to_output_genes;
my %entrez_to_output_genes;
my %ensembl_to_output_genes;

while(<OUTGENES>){
  next if /^#/;
  next if /^[\r\n]/; # skip blank lines
  my ($symbol, $entrez, $ensembl, $type) = split /\t/;
  chomp($type);
  #my $standard_gene = $gene;
  #if(exists($genes{$gene})){
  #  $standard_gene = $genes{$gene};
  #}
  #$output_genes{$standard_gene} = $type;
  my @ensembls = split(/_/, $ensembl);
  $output_genes{"$symbol\t$entrez\t$ensembl"} = $type;
  $symbol_to_output_genes{$symbol} = "$symbol\t$entrez\t$ensembl";
  $entrez_to_output_genes{$entrez} = "$symbol\t$entrez\t$ensembl";
  foreach my $this_ensembl (@ensembls){
    $ensembl_to_output_genes{$this_ensembl} = "$symbol\t$entrez\t$ensembl";
  }

}
close OUTGENES;


# ============================= #
# process the cell lines header #
# ============================= #

open IN, "< $exprn_data" or die "Can't read expression data: $! \n";

# store cell line names in array
# Note that there is no header for the gene column,
# cell line names start in first position
my $cell_line_row = <IN>;
chomp($cell_line_row);
my @cell_lines = split(/,/, $cell_line_row);

# convert cell lines names to standard cell line
# names. e.g. A-549 => A549_LUNG. Use @cell_lines
# to look up the standard cell line name for each
# column of expression data
my %standard_cell_line_header;
my %standard_cell_line_tissues;
foreach my $cell_line (@cell_lines){
  $standard_cell_line_header{$cell_line} = $cell_lines{$cell_line};
  my $tissue = $cell_lines{$cell_line};
  $tissue =~ s/^[^_]+_//;
  $standard_cell_line_tissues{$cell_line} = $tissue;
}

# sort the standardised cell lines by tissue
my @standard_cell_lines_by_tissue = sort { $standard_cell_line_tissues{$a} cmp $standard_cell_line_tissues{$b} } keys %standard_cell_line_tissues;
# Note that the keys are the cell line names from the header, not the standardised cell line names



# ================== #
# set up the outputs #
# ================== #

# open the output file
open OUT, "> $output" or die "Can't write output file $output: $!\n";

my $header = "GeneSymbol" . '_' . "EntrezGeneID" . '_' . "EnsemblGeneID";
foreach my $standard_cell_lines_by_tissue (@standard_cell_lines_by_tissue){
  $header .= "\t$standard_cell_line_header{$standard_cell_lines_by_tissue}";
}
$header .= "\n";


# print a header for the output. First three cols are for gene IDs
print OUT $header;


# =========================== #
# process the expression data #
# =========================== #

# process the expression data line by line
while(<IN>){
  
  my %cell_line_gene_expression = ();
  
  my $line = $_;
  chomp($line);
  
  my @expression_values = split(/,/, $line);
  
  my $gene = shift(@expression_values);

  # look up gene symbol in dictionary and store standardised name
  # we may need to clean up the symbols. Many have /[ ]+\/\/\/.+/
  #$gene =~ s/[ ]+\/\/\/.+//;
#  if($gene =~ /^[^ ]+[ ]+\/\/\/[^ ]+(.+)/){
#    $gene = $1;
#  }
  
  
  
  my $standard_gene_symbol = undef;
  if(exists $genes{$gene}){
    $standard_gene_symbol = $genes{$gene};
  }
  else{
###    warn "$gene is not in the gene symbol dictionary... Skipping.\n";
    next;
  }
  
  my $output_gene = undef;
  if(exists $symbol_to_output_genes{$standard_gene_symbol}){
    $output_gene = $symbol_to_output_genes{$standard_gene_symbol};
  }
  else{
#    print "skipping gene $standard_gene_symbol - not in the output genes list\n";
    next;
  }
    
  # save the expression values to a hash with cellline+gene
  # as the keys. We need to assume that @expression_values is
  # in the same order as @cell_lines so we can associate them
  # in the %expression_values hash. Check they have the same
  # number of elements...
  
  die "number of expression values and cell lines not equal for $gene\n" if $#expression_values != $#cell_lines;
  
  # store expression value with standard cell line + standard gene as key
  for(my $i = 0; $i <= $#expression_values; $i ++){
    $cell_line_gene_expression{$cell_lines[$i] . "\t" . $output_gene} = $expression_values[$i];
  }
  
  # write out the expression values sorted by tissue...
  my $gene_for_row = $output_gene;
  $gene_for_row =~ s/\t/_/g;
  
  print OUT $gene_for_row;
  
  foreach my $cell_line_to_print (@standard_cell_lines_by_tissue){
    my $cell_line_gene_expression_key = $cell_line_to_print . "\t" . $output_gene;
    print OUT "\t$cell_line_gene_expression{$cell_line_gene_expression_key}";
  }
  print OUT "\n";
    
}

close IN;
close OUT;


sub usage() {
  my $usage =<<END;
#---------------------------#
#  prep_expression_data.pl  #
#---------------------------#

James Campbell (jamesc\@icr.ac.uk)

Usage:
perl prep_expression_data.pl [options]

Options
--help			Display this message and quit
--exprn_data	Path to the Sanger expression data (EN input) [required]
--genes			Path to the gene name dictionary [required]
--cell_lines	Path to the cell line name dictionary [required]
--output		Path to output. Defaults to ccle_data.proc [optional]
END

  print $usage;
}



