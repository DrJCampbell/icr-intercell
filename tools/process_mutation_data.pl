#!/usr/bin/perl -w
# ============================================= #
# process_mutation_data.pl
# read data from CCLE etc, count mutations by type
# for each gene and cell line
# jamesc@icr.ac.uk, 20th May 2014
# ============================================= #

use strict;
use Getopt::Long;

my ($help, $ccle_data, $cosmic_data, $mut_freqs, $genes, $cell_lines, $output, $output_genes_list, $vcf_data);

GetOptions (
  "ccle_data=s" => \$ccle_data,
  "cosmic_data=s" => \$cosmic_data,
  "vcf_data=s" => \$vcf_data,
  "mut_freqs=s" => \$mut_freqs,
  "genes=s" => \$genes,
  "cell_lines=s" => \$cell_lines,
  "output=s" => \$output,
  "output_genes_list=s" => \$output_genes_list,
  "help" => \$help,
 );

# what to do with $output_genes_list?
# give a list of gene symbols,
# convert to the standardised list
# include in the list if the gene
# is a TSG or an OG. Use this info
# to determine what to count when
# outputing matrices



# print usage message if requested or no args supplied
if(defined($help)) {
  &usage;
  exit(0);
}

$output = "combined_mutation_data_table.txt" unless defined $output;


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

# read in the list of genes we want to output at the end and standardise
open OUTGENES, "< $output_genes_list" or die "Can't read list of genes to output from file $output_genes_list: $!\n";
my %output_genes;
while(<OUTGENES>){
  next if /^#/;
  my ($gene, $type) = split /\t/;
  chomp($type);
  my $standard_gene = $gene;
  if(exists($genes{$gene})){
    $standard_gene = $genes{$gene};
  }
  $output_genes{$standard_gene} = $type;
}
close OUTGENES;


# types of mutation consequences present in file:
my %mutation_consequences = (
	"De_novo_Start_InFrame" => "aa_sub", # These cover CCLE consequences
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
	"Start_Codon_Del" => "trunc",
	"Complex - deletion inframe" => "aa_sub", # These cover COSMIC consequences
	"Complex - insertion inframe" => "aa_sub",
	"Deletion - In frame" => "aa_sub",
	"Insertion - In frame" => "aa_sub",
	"Substitution - Missense" => "aa_sub",
	"Nonstop extension" => "aa_sub",
	"Complex - frameshift" => "trunc",
	"Deletion - Frameshift" => "trunc",
	"Insertion - Frameshift" => "trunc",
	"No detectable mRNA/protein" => "trunc",
	"Substitution - Nonsense" => "trunc",
	"Substitution - coding silent" => "unknown",
	"Unknown" => "unknown",
	"aa_sub" => "aa_sub",
	"mirna" => "mirna",
	"non_coding" => "non_coding",
	"silent" => "silent",
	"trunc" => "trunc",
	"unknown" => "unknown"
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


# these will store all cell lines and genes seen across
# the various data sets being integrated
my %master_cell_lines_seen;
my %master_genes_seen;



# ===================== #
# Process the CCLE data #
# ===================== #


my %cell_lines_seen;	# store standard cell line names for all seen
my %genes_seen;			# store standard gene names for all seen
my %mutations_seen;		# store cell line * genome mutation to prevent counting repeats



my %ccle_truncs;		# store truncating counts for cell*gene
my %ccle_rec_mis;		# store recurrent missense counts for cell*gene
my %ccle_other;			# store other counts for cell*gene

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
#  my $var_type = $fields[9];
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
#    print "already seen $sample_genome_change_key - skipping.\n";
    next;
  }
  else{
    $mutations_seen{$sample_genome_change_key} = 1;
  }
  
  my $ccle_key = "$standard_cell_line\t$standard_gene";
  
  $cell_lines_seen{$standard_cell_line} = 1;
  $genes_seen{$standard_gene} = 1;

  # update the master record of all cell lines and genes seen
  $master_cell_lines_seen{$standard_cell_line} = 1;
  $master_genes_seen{$standard_gene} = 1;
  
  
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


# ======================= #
# Process the COSMIC data #
# ======================= #


# reset these so we can reuse for COSMIC
%cell_lines_seen = ();
%genes_seen = ();
%mutations_seen = ();

my %cos_truncs;		# store truncating counts for cell*gene
my %cos_rec_mis;	# store recurrent missense counts for cell*gene
my %cos_other;		# store other counts for cell*gene

# Open and read mutations file
open COS, "< $cosmic_data" or die "Can't read mutations file $cosmic_data: $!\n";

$header = <COS>;

while(<COS>){

  # get gene, cell line, prot. mut., consequence etc.
  # [1]		Entrez_Gene_Id
  # [8]		Variant_Classification	// Missense_Mutation, Intron, Frame_Shift_Del
  # [9]		Variant_Type			// DEL, SNP
  # [15]	Tumor_Sample_Barcode	// JHUEM2_ENDOMETRIUM
  # [32]	Genome_Change
  # [37]	Protein_Change
	
  # column indexes for COSMIC
  # [0]   Gene name
  # [3]   Sample name
  # [13]  Mutation AA
  # [14]  Mutation Description
  # [15]  Mutation zygosity
  # [18]  Mutation GRCh37 genome position
  my @fields = split(/\t/);
  my $entrez_gene = $fields[0];
  my $var_class = $fields[14];
#  my $var_type = $fields[];
  my $sample = $fields[3];
  my $genome_change = $fields[18];
  my $prot_change = $fields[13];
  my $mutation_zygosity = $fields[15];
  
  # sometimes var_class is blank and we can do nothing with it... skip
  next if $var_class eq '';
  
  # lookup gene and cell line name in dictionaries
  my $standard_gene = $entrez_gene;
  # The COSMIC data often has the gene ID concatenated to the 
  # transcript ID. e.g. CEACAM18_ENST00000451626
  # strip out anything after and '_' from entrez_gene
  $standard_gene =~ s/_.+//;
  
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
#    print "already seen $sample_genome_change_key - skipping.\n";
    next;
  }
  else{
    $mutations_seen{$sample_genome_change_key} = 1;
  }
  
  my $cos_key = "$standard_cell_line\t$standard_gene";
  
  $cell_lines_seen{$standard_cell_line} = 1;
  $genes_seen{$standard_gene} = 1;
  
  # update the master record of all cell lines and genes seen
  $master_cell_lines_seen{$standard_cell_line} = 1;
  $master_genes_seen{$standard_gene} = 1;

#  
#  print "VARCLASS = $var_class\n";
#  
  
  if($mutation_consequences{$var_class} eq "trunc"){	# if mutation is truncating
    $cos_truncs{$cos_key} ++;							# increment cellline+gene truncating count
  }
  elsif($mutation_consequences{$var_class} eq "aa_sub"){
    # check if the mutation is recurrent as per Davoli
    my $davoli_key = "$standard_gene\t$prot_change";
    if(exists($mut_freqs{$davoli_key})){
      $cos_rec_mis{$cos_key} ++;
    }
    else{
      $cos_other{$cos_key} ++;
    }
  }
} # finished reading file, close
close COS;



# ============================ #
# Process the ICR VCF+VEP data #
# ============================ #



# [0] original_file
# [1] cell_line
# [2] CHROM
# [3] POS
# [4] ID
# [5] REF
# [6] ALT
# [7] QUAL
# [8] FILTER
# [9] INFO
# [10] FORMAT
# [10] ES2
# [11] JapHapMap
# [12] Uploaded_variation
# [13] Location
# [14] Allele
# [15] Gene
# [16] Feature
# [17] Feature_type
# [18] Consequence
# [19] cDNA_position
# [20] CDS_position
# [20] Protein_position
# [21] Amino_acids
# [22] Codons
# [23] Existing_variation
# [24] Extra
# [25] most_severe_consequence
# [26] protein_mutation







# reset these so we can reuse for VCF+VEP data
%cell_lines_seen = ();
%genes_seen = ();
%mutations_seen = ();

my %icr_truncs;		# store truncating counts for cell*gene
my %icr_rec_mis;	# store recurrent missense counts for cell*gene
my %icr_other;		# store other counts for cell*gene

# Open and read mutations file
open ICR, "< $vcf_data" or die "Can't read mutations file $vcf_data: $!\n";

$header = <ICR>;

while(<ICR>){

  # get gene, cell line, prot. mut., consequence etc.
  # [1]		Entrez_Gene_Id
  # [8]		Variant_Classification	// Missense_Mutation, Intron, Frame_Shift_Del
  # [9]		Variant_Type			// DEL, SNP
  # [15]	Tumor_Sample_Barcode	// JHUEM2_ENDOMETRIUM
  # [32]	Genome_Change
  # [37]	Protein_Change

  # column indexes for ICR VCF+VEP data
  # [24] Extra
  # [25] most_severe_consequence
  # [1] cell_line
  # [12] Uploaded_variation
  # [26] protein_mutation


  my @fields = split(/\t/);
  my $entrez_gene = $fields[24];
  my $var_class = $fields[25];
#  my $var_type = $fields[];
  my $sample = $fields[1];
  my $genome_change = $fields[12];
  my $prot_change = $fields[26];
#  my $mutation_zygosity = $fields[15];
  
  # sometimes var_class is blank and we can do nothing with it... skip
  next if $var_class eq '';
  
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
#    print "already seen $sample_genome_change_key - skipping.\n";
    next;
  }
  else{
    $mutations_seen{$sample_genome_change_key} = 1;
  }
  
  my $icr_key = "$standard_cell_line\t$standard_gene";
  
  $cell_lines_seen{$standard_cell_line} = 1;
  $genes_seen{$standard_gene} = 1;
  
  # update the master record of all cell lines and genes seen
  $master_cell_lines_seen{$standard_cell_line} = 1;
  $master_genes_seen{$standard_gene} = 1;

  if($mutation_consequences{$var_class} eq "trunc"){	# if mutation is truncating
    $icr_truncs{$icr_key} ++;							# increment cellline+gene truncating count
  }
  elsif($mutation_consequences{$var_class} eq "aa_sub"){
    # check if the mutation is recurrent as per Davoli
    my $davoli_key = "$standard_gene\t$prot_change";
    if(exists($mut_freqs{$davoli_key})){
      $icr_rec_mis{$icr_key} ++;
    }
    else{
      $icr_other{$icr_key} ++;
    }
  }
} # finished reading file, close
close ICR;



# ==================================== #
# Process the hashes to format outputs
# ==================================== #


# get lists of all genes and cell lines seen in any data set
my @cell_lines_seen = keys %master_cell_lines_seen;
my @genes_seen = keys %master_genes_seen;


my $output_data = '';

my $output_header = "cell_line";
foreach my $seen_gene (@genes_seen){
  $output_header .= "\t$seen_gene (trunc)\t$seen_gene (rec_mis)\t$seen_gene (other)";
}
$output_header .= "\n";


foreach my $seen_cell_line (@cell_lines_seen){
  $output_data .= $seen_cell_line;
  foreach my $seen_gene (@genes_seen){
    my $hash_key = "$seen_cell_line\t$seen_gene";
    if(exists($ccle_truncs{$hash_key}) || exists($cos_truncs{$hash_key})){
      $output_data .= "\t1";
    }
    else{
      $output_data .= "\t0";
    }
    if(exists($ccle_rec_mis{$hash_key}) || exists($cos_rec_mis{$hash_key})){
      $output_data .= "\t1";
    }
    else{
      $output_data .= "\t0";
    }
    if(exists($ccle_other{$hash_key}) || exists($cos_other{$hash_key})){
      $output_data .= "\t1";
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
#----------------------------#
#  process_mutation_data.pl  #
#----------------------------#

James Campbell (jamesc\@icr.ac.uk)

Usage:
perl process_mutation_data.pl [options]

Options
--help			Display this message and quit
--ccle_data		Path to the CCLE maf file [required]
--cosmic_data	Path to the COSMIC mutations file [required]
--vcf_data		Path to the ICR mutation data [required]
--mut_freqs		Path to the processed Davoli data [required]
--genes			Path to the gene name dictionary [required]
--cell_lines		Path to the cell line name dictionary [required]
--output		Path to output. Defaults to ccle_data.proc [optional]
END

  print $usage;
}







