#!/usr/bin/perl -w
# ========================================= #
# process_mutation_data.pl
# read data from CCLE etc, count mutations
# by type for each gene and cell line
# jamesc@icr.ac.uk, 25th July 2014
# ========================================= #

use strict;
use Getopt::Long;

my ($help, $ccle_mut_data, $cosmic_mut_data, $biankin_mut_data, $mut_freqs, $genes, $cell_lines, $output, $output_genes_list, $vcf_mut_data, $ccle_cna_data, $wtsi_mut_data, $wtsi_expression_z_data);

GetOptions (
  "ccle_muts=s" => \$ccle_mut_data,
  "cosmic_muts=s" => \$cosmic_mut_data,
  "vcf_muts=s" => \$vcf_mut_data,
  "ccle_cna=s" => \$ccle_cna_data,
  "wtsi_muts=s" => \$wtsi_mut_data,
#  "cosmic_cna=s" => \$cosmic_cna_data,
  "biankin_muts=s" => \$biankin_mut_data,
  "wtsi_expression_z_data=s" => \$wtsi_expression_z_data,
  "mut_freqs=s" => \$mut_freqs,
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

$output = "combined_mutation_data_table.txt" unless defined $output;


# ================================================================ #
# this will be a hash of hashes to store the mutation details for
# all gene_celllines in each data set. Primary keys are cell_line
# and gene. Secondary keys are data sources (eg. comsic_mut)
# ================================================================ #
 my %details_table = ();


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

# =================================== #
# read in the mutation frequency data #
# =================================== #

open MFRQ, "< $mut_freqs" or die "Can't read mutation frequency information file $mut_freqs: $!\n";
my %mut_freqs;
while(<MFRQ>){
  my ($gene_name, $prot_pos, $freq, $gene_seq_count) = split(/\t/);
  chomp($gene_seq_count);
  my $three_percent_of_seq_counts = $gene_seq_count * 0.03;
  my $threshold = (3, $three_percent_of_seq_counts)[3 < $three_percent_of_seq_counts]; # find max of 3 or $three_percent_of_seq_counts
  my $key = "$gene_name\t$prot_pos";
  if($freq >= $threshold){
    $mut_freqs{$key} = $freq;
  }
}
close MFRQ;


# ==================================== #
# read in the list of genes we want to #
# output at the end and standardise    #
# ==================================== #

# we have three identifiers for each gene
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
  my @ensembls = split(/_/, $ensembl);
  $output_genes{"$symbol\t$entrez\t$ensembl"} = $type;
  $symbol_to_output_genes{$symbol} = "$symbol\t$entrez\t$ensembl";
  $entrez_to_output_genes{$entrez} = "$symbol\t$entrez\t$ensembl";
  foreach my $this_ensembl (@ensembls){
    $ensembl_to_output_genes{$this_ensembl} = "$symbol\t$entrez\t$ensembl";
  }
}
close OUTGENES;


# ================================================ #
# types of mutation consequences present in files: #
# ================================================ #
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
	"aa_sub" => "aa_sub", # These are used in the modifed VCF+VEP calls for ovarian data
	"mirna" => "mirna",
	"non_coding" => "non_coding",
	"silent" => "silent",
	"trunc" => "trunc",
	"unknown" => "unknown",
	"initiator_codon_variant" => "aa_sub",		# These are used for the Biankin VEP data for the TKCC pancreatic samples
	"initiator_codon_variant,splice_region_variant" => "aa_sub",
	"missense_variant" => "aa_sub",
	"missense_variant,NMD_transcript_variant" => "aa_sub",
	"missense_variant,splice_region_variant" => "aa_sub",
	"missense_variant,splice_region_variant,NMD_transcript_variant" => "aa_sub",
	"splice_acceptor_variant" => "trunc",
	"splice_acceptor_variant,nc_transcript_variant" => "trunc",
	"splice_acceptor_variant,NMD_transcript_variant" => "trunc",
	"splice_donor_variant" => "trunc",
	"splice_donor_variant,nc_transcript_variant" => "trunc",
	"splice_donor_variant,NMD_transcript_variant" => "trunc",
	"splice_donor_variant,non_coding_exon_variant,nc_transcript_variant" => "trunc",
	"stop_gained" => "trunc",
	"stop_gained,NMD_transcript_variant" => "trunc",
	"stop_gained,splice_region_variant" => "trunc",
	"stop_gained,splice_region_variant,NMD_transcript_variant" => "trunc",
	"stop_lost" => "aa_sub",
	"3_prime_UTR_variant" => "non_coding",
	"3_prime_UTR_variant,NMD_transcript_variant" => "non_coding",
	"5_prime_UTR_variant" => "non_coding",
	"5_prime_UTR_variant,NMD_transcript_variant" => "non_coding",
	"coding_sequence_variant" => "non_coding",
	"downstream_gene_variant" => "non_coding",
	"incomplete_terminal_codon_variant,coding_sequence_variant" => "non_coding",
	"intergenic_variant" => "non_coding",
	"intron_variant" => "non_coding",
	"intron_variant,nc_transcript_variant" => "non_coding",
	"intron_variant,NMD_transcript_variant" => "non_coding",
	"non_coding_exon_variant,nc_transcript_variant" => "non_coding",
	"regulatory_region_variant" => "non_coding",
	"splice_region_variant,3_prime_UTR_variant" => "non_coding",
	"splice_region_variant,3_prime_UTR_variant,NMD_transcript_variant" => "non_coding",
	"splice_region_variant,5_prime_UTR_variant" => "non_coding",
	"splice_region_variant,5_prime_UTR_variant,NMD_transcript_variant" => "non_coding",
	"splice_region_variant,intron_variant" => "non_coding",
	"splice_region_variant,intron_variant,nc_transcript_variant" => "non_coding",
	"splice_region_variant,intron_variant,NMD_transcript_variant" => "non_coding",
	"splice_region_variant,non_coding_exon_variant,nc_transcript_variant" => "non_coding",
	"splice_region_variant,synonymous_variant" => "non_coding",
	"splice_region_variant,synonymous_variant,NMD_transcript_variant" => "non_coding",
	"synonymous_variant" => "non_coding",
	"synonymous_variant,NMD_transcript_variant" => "non_coding",
	"TF_binding_site_variant" => "non_coding",
	"upstream_gene_variant" => "non_coding",
	"Complex - deletion inframe" => "aa_sub", # These are used for the WTSI data set
	"Complex - frameshift" => "trunc",
	"Complex - insertion inframe" => "aa_sub",
	"Deletion - Frameshift" => "trunc",
	"Deletion - In frame" => "aa_sub",
	"Essential Splice" => "trunc",
	"Gene Fusion" => "unknown",
	"Genomic Amplification" => "amplification",
	"Homozygous Deletion" => "deletion",
	"Insertion - Frameshift" => "trunc",
	"Insertion - In frame" => "aa_sub",
	"microRNA - Deletion" => "unknown",
	"microRNA - Insertion" => "unknown",
	"microRNA - Substitution" => "unknown",
	"Nonstop extension" => "aa_sub",
	"Substitution - coding silent" => "non_coding",
	"Substitution - Missense" => "aa_sub",
	"Substitution - Nonsense" => "trunc",
	"Unknown" => "unknown",
	"complex sub" => 'aa_sub',				# Added for Osteosarcoma lines
	"downstream" => 'non_coding',
	"ess splice" => 'trunc',
	"frameshift" => 'trunc',
	"inframe" => 'aa_sub',
	"inframe_deletion" => 'aa_sub',
	"inframe_insertion" => 'aa_sub',
	"intergenic" => 'non_coding',
	"intronic" => 'non_coding',
	"missense" => 'aa_sub',
	"nonsense" => 'trunc',
	"Silent" => 'non_coding',
	"start-lost" => 'aa_sub',
	"stop_lost" => 'aa_sub',
	"stop-lost" => 'aa_sub',
	"unknown" => 'non_coding',
	"upstream" => 'non_coding'
	);



# these will store all cell lines and genes seen
# across the various data sets being integrated
my %master_cell_lines_seen;
my %master_genes_seen;


# ===================== #
# Process the CCLE data #
# ===================== #

my %mutations_seen;		# store cell line * genome mutation to prevent counting repeats
my %ccle_truncs;		# store truncating counts for cell * gene
my %ccle_rec_mis;		# store recurrent missense counts for cell * gene
my %ccle_other;			# store other counts for cell * gene

if(defined $ccle_mut_data){
	# Open and read mutations file
	open MUTS, "< $ccle_mut_data" or die "Can't read mutations file $ccle_mut_data: $!\n";
	my $header = <MUTS>;
	while(<MUTS>){
	  my @fields = split(/\t/);
	  my $entrez_gene = $fields[1];			# [1]		Entrez_Gene_Id
	  my $var_class = $fields[8];			# [8]		Variant_Classification	// Missense_Mutation, Intron, Frame_Shift_Del
	  my $sample = $fields[15];				# [15]	Tumor_Sample_Barcode	// JHUEM2_ENDOMETRIUM
	  my $genome_change = $fields[32];		# [32]	Genome_Change
	  my $prot_change = $fields[37];		# [37]	Protein_Change
	
	  # skip processing this variant if the gene is not in the set of output genes
	  next unless exists $entrez_to_output_genes{$entrez_gene};
	  
	  # lookup gene and cell line name in dictionaries
	  my $standard_gene = $entrez_to_output_genes{$entrez_gene};
	  my $standard_cell_line = $sample;
	  if(exists($cell_lines{$sample})){
		$standard_cell_line = $cell_lines{$sample};
	  }
	  else{
		warn "Cell line $sample not found in the cell line dictionary. You may need to add it.\n";
	  }
	
	  # store a hash key for every cell line * genome_change * consequence seen
	  # skip if already counted...
	  my $sample_genome_change_key = $standard_cell_line . "_" . $genome_change . "_" . $var_class;
	  if(exists($mutations_seen{$sample_genome_change_key})){
		next;
	  }
	  else{
		$mutations_seen{$sample_genome_change_key} = 1;
	  }
	  
	  my $ccle_key = "$standard_cell_line\t$standard_gene";
	  
	  # update the master record of all cell lines and genes seen
	  # Note that for the %master_cell_lines_seen hash, we need to
	  # know if the cell line was seen with mutation data or copy
	  # number data. 
	  if(exists($master_cell_lines_seen{$standard_cell_line})){
		$master_cell_lines_seen{$standard_cell_line} .= "\tmut";
	  }
	  else{
		$master_cell_lines_seen{$standard_cell_line} = "\tmut";
	  }
	  $master_genes_seen{$standard_gene} = 1;
	  
	  if($mutation_consequences{$var_class} eq "trunc"){	# if mutation is truncating
		$ccle_truncs{$ccle_key} ++;							# increment cellline+gene truncating count
	  }
	  elsif($mutation_consequences{$var_class} eq "aa_sub"){
		# check if the mutation is recurrent as per Davoli
		my @identifiers = split(/\t/, $standard_gene);
		my $prot_position = $prot_change;
		$prot_position =~ s/^p\.[A-Z]+//;
		$prot_position =~ s/[^\d].*//;
		my $davoli_key = "$identifiers[0]\t$prot_position";
		if(exists($mut_freqs{$davoli_key})){
		  $ccle_rec_mis{$ccle_key} ++;
		}
		else{
		  $ccle_other{$ccle_key} ++;
		}
	  }
	} # finished reading file, close
	close MUTS;
}

# ======================= #
# Process the COSMIC data #
# ======================= #

%mutations_seen = ();	# reset these so we can reuse for COSMIC
my %cos_truncs;			# store truncating counts for cell*gene
my %cos_rec_mis;		# store recurrent missense counts for cell*gene
my %cos_other;			# store other counts for cell*gene

if(defined $cosmic_mut_data){
	# Open and read mutations file
	open COS, "< $cosmic_mut_data" or die "Can't read mutations file $cosmic_mut_data: $!\n";
	my $header = <COS>;
	while(<COS>){
	  my @fields = split(/\t/);
	  my $entrez_gene = $fields[0];  		# this is actually a gene symbol
	  my $var_class = $fields[14];			# [14]  Mutation Description
	  my $sample = $fields[3];				# [3]   Sample name
	  my $genome_change = $fields[18];		# [18]  Mutation GRCh37 genome position
	  my $prot_change = $fields[13];		# [13]  Mutation AA
	  my $mutation_zygosity = $fields[15];	# [15]  Mutation zygosity
	  
	  # sometimes var_class is blank and we can do nothing with it... skip
	  next if $var_class eq '';
	  
	  # lookup gene and cell line name in dictionaries
	  $entrez_gene =~ s/_.+//;	# The COSMIC data often has the gene ID concatenated to the 
								# transcript ID. e.g. CEACAM18_ENST00000451626
								# strip out anything after '_' from entrez_gene
	
	  my $standard_gene = $entrez_gene;
	  if(exists($genes{$entrez_gene})){
		$standard_gene = $genes{$entrez_gene};
	  }
	  $entrez_gene = $standard_gene ; # this may look nuts but we are going to munge standard_gene later and want to still keep the gene symbol handy...
	  
	  my $standard_cell_line = $sample;
	  if(exists($cell_lines{$sample})){
		$standard_cell_line = $cell_lines{$sample};
	  }
	  else{
		warn "Cell line $sample not found in the cell line dictionary. You may need to add it.\n";
	  }
	  
	  # skip processing this variant if the gene symbol is not in the set of output genes
	  next unless exists $symbol_to_output_genes{$standard_gene};
	
	  # change $standard_gene from the symbol to the full trio of IDs and keep a record of
	  # the standardised gene symbol for later use in $entrez_gene;
	  $standard_gene = $symbol_to_output_genes{$standard_gene};
	
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
	  
	  # update the master record of all cell lines and genes seen
	  if(exists($master_cell_lines_seen{$standard_cell_line})){
		$master_cell_lines_seen{$standard_cell_line} .= "\tmut";
	  }
	  else{
		$master_cell_lines_seen{$standard_cell_line} = "\tmut";
	  }
	  $master_genes_seen{$standard_gene} = 1;
	
	  my $prot_position = $prot_change;
	  $prot_position =~ s/^p\.[A-Z]+//;
	  $prot_position =~ s/[^\d].*//;  
	  my $davoli_key = "$entrez_gene\t$prot_position";
	  
	  if($mutation_consequences{$var_class} eq "trunc"){	# if mutation is truncating
		$cos_truncs{$cos_key} ++;							# increment cellline+gene truncating count
	  }
	  elsif($mutation_consequences{$var_class} eq "aa_sub"){
		# check if the mutation is recurrent as per Davoli
		if(exists($mut_freqs{$davoli_key})){
		  $cos_rec_mis{$cos_key} ++;
		}
		else{
		  $cos_other{$cos_key} ++;
		}
	  }
	
	
	  # collect the details for later output:
	  if(exists $details_table{$cos_key}{"cos_mut"}){
		$details_table{$cos_key}{"cos_mut"} .=  "; $prot_change";
		$details_table{$cos_key}{"cos_type"} .=  "; $mutation_consequences{$var_class}";
		if(exists($mut_freqs{$davoli_key})){
		  $details_table{$cos_key}{"cos_davoli_freq"} .= "; $mut_freqs{$davoli_key}";
		}
		else{
		  $details_table{$cos_key}{"cos_davoli_freq"} .= "; NA";
		}
	  }
	  else{
		$details_table{$cos_key}{"cos_mut"} =  $prot_change;
		$details_table{$cos_key}{"cos_type"} .=  $mutation_consequences{$var_class};
		if(exists($mut_freqs{$davoli_key})){
		  $details_table{$cos_key}{"cos_davoli_freq"} = "$mut_freqs{$davoli_key}";
		}
		else{
		  $details_table{$cos_key}{"cos_davoli_freq"} = "NA";
		}
	  }
	  
	  
	  
	} # finished reading file, close
	close COS;
}


# ============================ #
# Process the ICR VCF+VEP data #
# ============================ #

%mutations_seen = ();	# reset these so we can reuse for VCF+VEP data
my %icr_truncs;		# store truncating counts for cell*gene
my %icr_rec_mis;	# store recurrent missense counts for cell*gene
my %icr_other;		# store other counts for cell*gene

if(defined $vcf_mut_data){
	# Open and read mutations file
	open ICR, "< $vcf_mut_data" or die "Can't read mutations file $vcf_mut_data: $!\n";
	
	my $header = <ICR>;
	
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
	  my $entrez_gene = $fields[26]; # this is the Gene symbol... data are old and many ENSG IDs are not in the current dictionary
	  my $var_class = $fields[27];
	  my $sample = $fields[1];
	  my $genome_change = $fields[13];
	  my $prot_change = $fields[28];
	  
	  # sometimes var_class is blank and we can do nothing with it... skip
	  next if $var_class eq '';
	  
	  # lookup gene and cell line name in dictionaries
	  my $standard_gene = $entrez_gene;
	  if(exists($genes{$entrez_gene})){
		$standard_gene = $genes{$entrez_gene};
	  }
	  
	  next unless exists $symbol_to_output_genes{$standard_gene};
	  $standard_gene = $symbol_to_output_genes{$standard_gene};
	
	  my $standard_cell_line = $sample;
	  if(exists($cell_lines{$sample})){
		$standard_cell_line = $cell_lines{$sample};
	  }
	  else{
		warn "Cell line $sample not found in the cell line dictionary. You may need to add it.\n";
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
	  
	 
	  # update the master record of all cell lines and genes seen
	  if(exists($master_cell_lines_seen{$standard_cell_line})){
		$master_cell_lines_seen{$standard_cell_line} .= "\tmut";
	  }
	  else{
		$master_cell_lines_seen{$standard_cell_line} = "\tmut";
	  }
	  $master_genes_seen{$standard_gene} = 1;
	
	
	  my $prot_position = $prot_change;
	  $prot_position =~ s/^p\.[A-Z]+//;
	  $prot_position =~ s/[^\d].*//;
	  my @identifiers = split(/\t/, $standard_gene);
	  my $davoli_key = "$identifiers[0]\t$prot_position";
	
	
	  if($mutation_consequences{$var_class} eq "trunc"){	# if mutation is truncating
		$icr_truncs{$icr_key} ++;							# increment cellline+gene truncating count
	  }
	  elsif($mutation_consequences{$var_class} eq "aa_sub"){
		# check if the mutation is recurrent as per Davoli
		
		if(exists($mut_freqs{$davoli_key})){
		  $icr_rec_mis{$icr_key} ++;
		}
		else{
		  $icr_other{$icr_key} ++;
		}
	  }
	  
	  # collect the details for later output:
	  if(exists $details_table{$icr_key}{"icr_mut"}){
		$details_table{$icr_key}{"icr_mut"} .=  "; $prot_change";
		$details_table{$icr_key}{"icr_type"} .=  "; $mutation_consequences{$var_class}";
		if(exists($mut_freqs{$davoli_key})){
		  $details_table{$icr_key}{"icr_davoli_freq"} .= "; $mut_freqs{$davoli_key}";
		}
		else{
		  $details_table{$icr_key}{"icr_davoli_freq"} .= "; NA";
		}
	  }
	  else{
		$details_table{$icr_key}{"icr_mut"} =  $prot_change;
		$details_table{$icr_key}{"icr_type"} .=  $mutation_consequences{$var_class};
		if(exists($mut_freqs{$davoli_key})){
		  $details_table{$icr_key}{"icr_davoli_freq"} = "$mut_freqs{$davoli_key}";
		}
		else{
		  $details_table{$icr_key}{"icr_davoli_freq"} = "NA";
		}
	  }
	  
	} # finished reading file, close
	close ICR;
}

# ============================ #
# Process the Biankin VEP data #
# ============================ #

%mutations_seen = ();	# reset these so we can reuse for Biankin VEP data
my %biankin_truncs;		# store truncating counts for cell*gene
my %biankin_rec_mis;	# store recurrent missense counts for cell*gene
my %biankin_other;		# store other counts for cell*gene

if(defined $biankin_mut_data){
	# Open and read mutations file
	open BIANKIN, "< $biankin_mut_data" or die "Can't read mutations file $biankin_mut_data: $!\n";
	
	my $header = <BIANKIN>;
	
	while(<BIANKIN>){
	
	  # get gene, cell line, prot. mut., consequence etc.
	  # [0]		sample
	  # [1]		genome_chnage
	  # [3]		ENSG ID
	  # [6]		consequence/variant classification
	  # [9]		protein position
	  # [10]	protein change
	
	  my @fields = split(/\t/);
	#  my $entrez_gene = $fields[3];	# this is the ENSG ID - Note - changed on 141008 to use the symbol because these ENSGs are not in the cancer gene classification file
	  my $entrez_gene = $fields[23];	# this is the symbol
	  my $var_class = $fields[6];
	  my $sample = $fields[0];
	  my $genome_change = $fields[1];
	  
	  my $prot_pos = $fields[9];
	  my @prot_ref_alt_aa = split(/\//, $fields[10]);
	  my $prot_change = undef;
	  if(defined($prot_ref_alt_aa[1])){
		$prot_change = 'p.' . $prot_ref_alt_aa[0] . $prot_pos . $prot_ref_alt_aa[1];
	  }
	  else{
		$prot_change = 'p.' . $prot_ref_alt_aa[0] . $prot_pos;
	  }
	  
	
	  
	  # sometimes var_class is blank and we can do nothing with it... skip
	  next if $var_class eq '';
	  
	
	  # lookup gene and cell line name in dictionaries
	  my $standard_gene = $entrez_gene;
	  if(exists($genes{$entrez_gene})){
		$standard_gene = $genes{$entrez_gene};
	  }
	  next unless exists $symbol_to_output_genes{$standard_gene};
	  $standard_gene = $symbol_to_output_genes{$standard_gene};
	
	  my $standard_cell_line = $sample;
	  if(exists($cell_lines{$sample})){
		$standard_cell_line = $cell_lines{$sample};
	  }
	  else{
		warn "Cell line $sample not found in the cell line dictionary. You may need to add it.\n";
	  }
	
	
	#  print "Prot: $prot_change\t";
	#  print "Stdgn: $standard_gene\t";
	#  print "Entrz: $entrez_gene\t";
	#  print "Consq: $var_class\n";
	  # store a hash key for every cell line * genome_change * consequence seen
	  # skip if already counted...
	  my $sample_genome_change_key = $standard_cell_line . "_" . $genome_change . "_" . $var_class;
	
	
	# This was dropped because if we are processing data with multiple transcripts for each mutation,
	# we want to test all the transcripts...
	#  if(exists($mutations_seen{$sample_genome_change_key})){
	#    print "already seen $sample_genome_change_key - skipping.\n";
	#    next;
	#  }
	#  else{
		$mutations_seen{$sample_genome_change_key} = 1;
	# }
	  
	  my $biankin_key = "$standard_cell_line\t$standard_gene";
	  
	  # update the master record of all cell lines and genes seen
	  if(exists($master_cell_lines_seen{$standard_cell_line})){
		$master_cell_lines_seen{$standard_cell_line} .= "\tmut";
	  }
	  else{
		$master_cell_lines_seen{$standard_cell_line} = "\tmut";
	  }
	  $master_genes_seen{$standard_gene} = 1;
	
	  my @identifiers = split(/\t/, $standard_gene);
	  my $davoli_key = "$identifiers[0]\t$prot_pos";
	
	  if($mutation_consequences{$var_class} eq "trunc"){	# if mutation is truncating
		$biankin_truncs{$biankin_key} ++;					# increment cellline+gene truncating count
	  }
	  elsif($mutation_consequences{$var_class} eq "aa_sub"){
		# check if the mutation is recurrent as per Davoli
	   
		if(exists($mut_freqs{$davoli_key})){
		  $biankin_rec_mis{$biankin_key} ++;
		  print "Found recurrent hit for $davoli_key\n";
		}
		else{
		  $biankin_other{$biankin_key} ++;
		}
	  }
	  
	  # collect the details for later output:
	  if(exists $details_table{$biankin_key}{"biankin_mut"}){
		$details_table{$biankin_key}{"biankin_mut"} .=  "; $prot_change";
		$details_table{$biankin_key}{"biankin_type"} .=  "; $mutation_consequences{$var_class}";
		if(exists($mut_freqs{$davoli_key})){
		  $details_table{$biankin_key}{"biankin_davoli_freq"} .= "; $mut_freqs{$davoli_key}";
		}
		else{
		  $details_table{$biankin_key}{"biankin_davoli_freq"} .= "; NA";
		}
	  }
	  else{
		$details_table{$biankin_key}{"biankin_mut"} =  $prot_change;
		$details_table{$biankin_key}{"biankin_type"} .=  $mutation_consequences{$var_class};
		if(exists($mut_freqs{$davoli_key})){
		  $details_table{$biankin_key}{"biankin_davoli_freq"} = "$mut_freqs{$davoli_key}";
		}
		else{
		  $details_table{$biankin_key}{"biankin_davoli_freq"} = "NA";
		}
	  }
	  
	} # finished reading file, close
	close BIANKIN;
}

# ====================================== #
# Process the WTSI mutation and CNA data #
# ====================================== #

%mutations_seen = ();	# reset these so we can reuse for COSMIC
my %wtsi_truncs;	# store truncating counts for cell*gene
my %wtsi_rec_mis;	# store recurrent missense counts for cell*gene
my %wtsi_other;		# store other counts for cell*gene

my %wtsi_cnas;		# store: none/amp/del/gain/loss

if(defined $wtsi_mut_data){
	# Open and read mutations file
	open WTSI, "< $wtsi_mut_data" or die "Can't read mutations file $wtsi_mut_data: $!\n";
	
	my $header = <WTSI>;
	
	while(<WTSI>){
	
	  # get gene, cell line, prot. mut., consequence etc.
	  # [0] ID_VARIANT
	  # [1] SAMPLE_NAME
	  # [2] CHR
	  # [3] GENOME_START
	  # [4] GENOME_STOP
	  # [5] STRAND
	  # [6] ALGORITHM
	  # [7] GENE_NAME
	  # [8] COSMIC_Name
	  # [9] TRANSCRIPT
	  # [10] WT
	  # [11] MT
	  # [12] CDS_SYNTAX
	  # [13] AA_MUT_SYNTAX
	  # [14] DESCRIPTION
	  # [15] QUALITY
	  # [16] ZYGOSITY
	  # [17] VERIF_STATUS_MINT
	  # [18] FATHMM Score
	  # [19] FATHMM Summary
	  # [20] Comment
	
	  my @fields = split(/\t/);
	  my $entrez_gene = $fields[7];  # this is actually a gene symbol
	  my $var_class = $fields[14];
	  my $sample = $fields[1];
	  my $genome_change = $fields[12];
	  my $prot_change = $fields[13];
	  my $mutation_zygosity = $fields[16];
	  
	  # sometimes var_class is blank and we can do nothing with it... skip
	  next if $var_class eq '' || $var_class eq 'Unknown';
	  
	  # lookup gene and cell line name in dictionaries
	
	  $entrez_gene =~ s/_.+//;	# The COSMIC data often has the gene ID concatenated to the 
								# transcript ID. e.g. CEACAM18_ENST00000451626
								# strip out anything after '_' from entrez_gene
	
	  my $standard_gene = $entrez_gene;
	  
	  if(exists($genes{$entrez_gene})){
		$standard_gene = $genes{$entrez_gene};
	  }
	  $entrez_gene = $standard_gene ; # this may look nuts but we are going to munge standard_gene later and want to still keep the gene symbol handy...
	  
	  my $standard_cell_line = $sample;
	  if(exists($cell_lines{$sample})){
		$standard_cell_line = $cell_lines{$sample};
	  }
	  else{
		warn "Cell line $sample not found in the cell line dictionary. You may need to add it.\n";
	  }
	
	  
	  # skip processing this variant if the gene symbol is not in the set of output genes
	  next unless exists $symbol_to_output_genes{$standard_gene};
	
	  # change $standard_gene from the symbol to the full trio of IDs and keep a record of
	  # the standardised gene symbol for later use in $entrez_gene;
	  $standard_gene = $symbol_to_output_genes{$standard_gene};
	
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
	
	  my $wtsi_key = "$standard_cell_line\t$standard_gene";
	  
	  # update the master record of all cell lines and genes seen
	  if(exists($master_cell_lines_seen{$standard_cell_line})){
		$master_cell_lines_seen{$standard_cell_line} .= "\tmut\tCNA";
	  }
	  else{
		$master_cell_lines_seen{$standard_cell_line} = "\tmut\tCNA";
	  }
	  $master_genes_seen{$standard_gene} = 1;
	
	  my $prot_position = $prot_change;
	  $prot_position =~ s/^p\.[A-Z]+//;
	  $prot_position =~ s/[^\d].*//;
	  my $davoli_key = "$entrez_gene\t$prot_position";
	
	  if($mutation_consequences{$var_class} eq "trunc"){	# if mutation is truncating
		$wtsi_truncs{$wtsi_key} ++;							# increment cellline+gene truncating count
	  }
	  elsif($mutation_consequences{$var_class} eq "aa_sub"){
		# check if the mutation is recurrent as per Davoli
		if(exists($mut_freqs{$davoli_key})){
		  $wtsi_rec_mis{$wtsi_key} ++;
		}
		else{
		  $wtsi_other{$wtsi_key} ++;
		}
	  }
	  elsif($mutation_consequences{$var_class} eq "deletion"){
		$wtsi_cnas{$wtsi_key} = -2; # hom-del
	  }
	  elsif($mutation_consequences{$var_class} eq "amplification"){
		$wtsi_cnas{$wtsi_key} = 2; # amp
		
		### DEBUG
		print "WTSI - found an AMP for $wtsi_key\n";
		
	  }
	  
	  # collect the details for later output:
	  if(exists $details_table{$wtsi_key}{"wtsi_mut"}){
		$details_table{$wtsi_key}{"wtsi_mut"} .=  "; $prot_change";
		$details_table{$wtsi_key}{"wtsi_type"} .=  "; $mutation_consequences{$var_class}";
		if(exists($mut_freqs{$davoli_key})){
		  $details_table{$wtsi_key}{"wtsi_davoli_freq"} .= "; $mut_freqs{$davoli_key}";
		}
		else{
		  $details_table{$wtsi_key}{"wtsi_davoli_freq"} .= "; NA";
		}
	  }
	  else{
		$details_table{$wtsi_key}{"wtsi_mut"} =  $prot_change;
		$details_table{$wtsi_key}{"wtsi_type"} .=  $mutation_consequences{$var_class};
		if(exists($mut_freqs{$davoli_key})){
		  $details_table{$wtsi_key}{"wtsi_davoli_freq"} = "$mut_freqs{$davoli_key}";
		}
		else{
		  $details_table{$wtsi_key}{"wtsi_davoli_freq"} = "NA";
		}
	  }
	  
	} # finished reading file, close
	close WTSI;
}


# ================================= #
# Process the CCLE gistic CNA data #
# ================================= #


# reset these so we can reuse for VCF+VEP data
%mutations_seen = ();

my %ccle_cnas;		# store: none/amp/del/gain/loss

if(defined $ccle_cna_data){
	# Open and read the CCLE CNA file
	open CCLECNA, "< $ccle_cna_data" or die "Can't read mutations file $ccle_cna_data: $!\n";
	
	my $header = <CCLECNA>;
	chomp($header);
	my @ccle_cna_celllines = split(/\t/,$header);
	
	while(<CCLECNA>){
	
	  my @fields = split(/\t/);
	  chomp($fields[$#fields]);
	
	  my $gene_id = shift @fields;
	  my $gene_symbol = shift @fields;
	 # if(exists($genes{$standard_gene})){
	 #   $standard_gene = $genes{$standard_gene};
	 # }
	  
	  unless(exists($entrez_to_output_genes{$gene_id})){
		print "Skipping $gene_id ($gene_symbol)\n";
		next;
	  }
	  
	  my $standard_gene = $entrez_to_output_genes{$gene_id};
	  
	  my $cell_line_counter = 2;	# We start counting from the index of the first cell line name
									# i.e. we skip the first two array elements with geneID and symbol
	  
	  foreach my $cna (@fields){
		
		my $standard_cell_line = $ccle_cna_celllines[$cell_line_counter];
		$cell_line_counter ++;
		
		if(exists($cell_lines{$standard_cell_line})){
		  $standard_cell_line = $cell_lines{$standard_cell_line};
		}
		else{
		  warn "Cell line $standard_cell_line not found in the cell line dictionary. You may need to add it.\n";
		}
		
		
		# update the master record of all cell lines and genes seen
		if(exists($master_cell_lines_seen{$standard_cell_line})){
			$master_cell_lines_seen{$standard_cell_line} .= "\tCNA";
		}
		else{
			$master_cell_lines_seen{$standard_cell_line} = "\tCNA";
		}
		$master_genes_seen{$standard_gene} = 1; 
		
		my $ccle_cna_key = "$standard_cell_line\t$standard_gene";
	
		$ccle_cnas{$ccle_cna_key} = $cna;
		
		# add the CCLE CNA to the details table hash
		if(exists $details_table{$ccle_cna_key}{"ccle_cna"}){
		  $details_table{$ccle_cna_key}{"ccle_cna"} .= "; $cna";
		}
		else{
		  $details_table{$ccle_cna_key}{"ccle_cna"} = "$cna";
		}    
	  }
	
	} # finished reading file, close
	close CCLECNA;
}

# =============================== #
# Process the expression z-scores
# =============================== #

my %wtsi_exprn_z;		# store: -1/0/1 for under/normal/overexpression

if(defined $wtsi_expression_z_data){
	# Open and read the CCLE CNA file
	open EXPRNZ, "< $wtsi_expression_z_data" or die "Can't read expression z-scores file $wtsi_expression_z_data: $!\n";
	
	
	# gene	expression.z	cell.line
	# CHEK2_11200_ENSG00000183765	NA	SW13_ADRENAL_GLAND
	# CHEK2_11200_ENSG00000183765	-2.00063333901991	NH12_AUTONOMIC_GANGLIA
	# CHEK2_11200_ENSG00000183765	-1.10868384234017	CHP212_AUTONOMIC_GANGLIA
	
	my $header = <EXPRNZ>;
	
	while(<EXPRNZ>){
	
	  my @fields = split(/\t/);
	  my $standard_gene = $fields[0];
	  my $exprnz = $fields[1];
	  my $standard_cell_line = $fields[2];
	  
	  chomp $standard_cell_line;
	  $standard_gene =~ s/^([^_]+)_([^_]+)_/$1\t$2\t/; # replace the underscores following the symbol and EntrezGeneID for tabs
	  next if $exprnz =~ /NA/;
	  
	  # lookup gene and cell line name in dictionaries
	
	  # skip processing this variant if the gene symbol is not in the set of output genes
	  if(!exists $output_genes{$standard_gene}){
	#    print "Skipping expression for gene: $standard_gene\n";
		next;
	  }
	
	  my $wtsi_key = "$standard_cell_line\t$standard_gene";
	  
	  # update the master record of all cell lines and genes seen
	  if(exists($master_cell_lines_seen{$standard_cell_line})){
		$master_cell_lines_seen{$standard_cell_line} .= "\texprn";
	  }
	  else{
		$master_cell_lines_seen{$standard_cell_line} = "\texprn";
	  }
	  $master_genes_seen{$standard_gene} = 1;
	
	  if($exprnz >= 2){	# overexpressed
		$wtsi_exprn_z{$wtsi_key} = 1;
	  }
	  elsif($exprnz <= -2){ # underexpressed
		$wtsi_exprn_z{$wtsi_key} = -1;
	  }
	  else{
		$wtsi_exprn_z{$wtsi_key} = 0;
	  }
	  
	  # add expression z-scores to the detail table hash
	  if(exists $details_table{$wtsi_key}{'exprnz'}){
		$details_table{$wtsi_key}{'exprnz'} .= "; $exprnz";
	  }
	  else{
		$details_table{$wtsi_key}{'exprnz'} = $exprnz;
	  }
	  
	  
	} # finished reading file, close
	close EXPRNZ;
}

# ==================================== #
# Process the hashes to format outputs
# ==================================== #


# get lists of all genes and cell lines seen in any data set
my @cell_lines_seen = keys %master_cell_lines_seen;

#my @genes_seen = keys %master_genes_seen;
my @genes_seen = keys %output_genes;		# the list of genes from CGC and C5000


# We need to indicate which data sets were available (mut, cna or expression) as a column
# in the output.

my $output_data = '';				# detailed output with separate columns for truncs, hotspots, others, CNAs and expression 
my $mutation_matrix = '';			# the functionally relevant changes
my $other_matrix = '';				# all changes regardless of relevance
my $mutation_classification = '';	# del=5,amp=4,trunc_hom=3,trunc_het=2,miss=1,overexpr=6,underexpr=7,wt=0

my $output_header = "cell_line\tdatasets";
my $matrix_header = "cell_line";
foreach my $seen_gene (@genes_seen){	# The CGC and C5000s sets
  my $seen_gene_header = $seen_gene;
  $seen_gene_header =~ s/\t/_/g;
  $output_header .= "\t$seen_gene_header (trunc)\t$seen_gene_header (rec_mis)\t$seen_gene_header (other)\t$seen_gene_header (gistic)\t$seen_gene_header (exprnz)";
  $matrix_header .= "\t$seen_gene_header";
}
$output_header .= "\n";
$matrix_header .= "\n";


#foreach my $seen_cell_line (@cell_lines_seen){ # from keys %master_cell_lines_seen
while(my ($seen_cell_line, $cell_line_seen_in_dataset) = each  %master_cell_lines_seen){
  my $datasets = "";
  if($cell_line_seen_in_dataset =~ /mut/){
    $datasets .= 'mut_'
  }
  if($cell_line_seen_in_dataset =~ /CNA/){
    $datasets .= 'CNA_'
  }
  if($cell_line_seen_in_dataset =~ /exprn/){
    $datasets .= 'expression_'
  }
  
  $output_data .= "$seen_cell_line\t$datasets";
  $mutation_matrix .= "$seen_cell_line";
  $other_matrix .= "$seen_cell_line";
  $mutation_classification .= "$seen_cell_line";
  
  foreach my $seen_gene (@genes_seen){
    
    my $hash_key = "$seen_cell_line\t$seen_gene";
    my $hom_del = 0;
    my $loss = 0;
    my $amp = 0;
    my $trunc = 0;
    my $rec_mis = 0;
    my $other = 0;
    my $exprn = 0;
    
    if(exists($ccle_truncs{$hash_key}) || exists($cos_truncs{$hash_key}) || exists($icr_truncs{$hash_key}) || exists($biankin_truncs{$hash_key}) || exists($wtsi_truncs{$hash_key})){
      $output_data .= "\t1";
      $trunc = 1;
    }
    else{
      $output_data .= "\t0";
    }
    
    if(exists($ccle_rec_mis{$hash_key}) || exists($cos_rec_mis{$hash_key}) || exists($icr_rec_mis{$hash_key}) || exists($biankin_rec_mis{$hash_key}) || exists($wtsi_rec_mis{$hash_key})){
      $output_data .= "\t1";
      $rec_mis = 1;
    }
    else{
      $output_data .= "\t0";
    }
    
    if(exists($ccle_other{$hash_key}) || exists($cos_other{$hash_key}) || exists($icr_other{$hash_key}) || exists($biankin_other{$hash_key}) || exists($wtsi_other{$hash_key})){
      $output_data .= "\t1";
      $other = 1;
    }
    else{
      $output_data .= "\t0";
    }
    
    if(exists($ccle_cnas{$hash_key}) || exists($wtsi_cnas{$hash_key})){
      if(exists($wtsi_cnas{$hash_key})){
        $output_data .= "\t$wtsi_cnas{$hash_key}";
      }
      elsif(exists($ccle_cnas{$hash_key})){
        $output_data .= "\t$ccle_cnas{$hash_key}";
      }
      
      
      if(exists($wtsi_cnas{$hash_key})){
        if($wtsi_cnas{$hash_key} == -2){
          $hom_del = 1;
        }
        elsif($wtsi_cnas{$hash_key} == 2){
          $amp = 1;
        }
      }
      if(exists($ccle_cnas{$hash_key})){
        if($ccle_cnas{$hash_key} == -2){
          $hom_del = 1;
        }
        elsif($ccle_cnas{$hash_key} == 2){
          $amp = 1;
        }
        elsif($ccle_cnas{$hash_key} == -1){
          $loss = 1;
        }
      }
    }
    else{
      $output_data .= "\tNA";
    }
   
    if(exists($wtsi_exprn_z{$hash_key})){
      $output_data .= "\t$wtsi_exprn_z{$hash_key}";
      $exprn = $wtsi_exprn_z{$hash_key};
    }
    else{
      $output_data .= "\t0";
    }
    
    # now decide how to classify the cell line / gene...
    
    my $matrix_value = 0;
    
    if($output_genes{$seen_gene} eq 'OG'){
      if($rec_mis == 1 || $amp == 1 || $exprn == 1){
        $matrix_value = 1;
      }
    }
    elsif($output_genes{$seen_gene} eq 'TSG'){
      if($rec_mis == 1 || $hom_del == 1 || $trunc == 1 || $exprn == -1){
        $matrix_value = 1;
      }
    }

    if($matrix_value == 0){
      $mutation_matrix .= "\t0";
      if($other == 1){
        $other_matrix .= "\t1";
      }
      else{
        $other_matrix .= "\t0";
      }
    }
    else{
      $mutation_matrix .= "\t1";
      $other_matrix .= "\t1";
    }
    
    
    # values for each type of mutation
    # 5: hom del in TSG
    # 4: amp in OG
    # 3: hom trunc in TSG
    # 2: het trunc
    # 1: rec missense
    # 7: underexpressed
    # 6: overexpressed
    # 0: none
    
    if($output_genes{$seen_gene} eq 'OG' && $amp == 1){
      $mutation_classification .= "\t4";
    }
    elsif($output_genes{$seen_gene} eq 'TSG' && $hom_del == 1){
      $mutation_classification .= "\t5";
    }
    elsif($output_genes{$seen_gene} eq 'TSG' && $loss == 1 && $trunc == 1){
      $mutation_classification .= "\t3";
    }
    elsif($output_genes{$seen_gene} eq 'TSG' && $trunc == 1){
      $mutation_classification .= "\t2";
    }
    elsif($rec_mis == 1){
      $mutation_classification .= "\t1";
    }
    elsif($output_genes{$seen_gene} eq 'OG' && $exprn == 1){
      $mutation_classification .= "\t6";
    }
    elsif($output_genes{$seen_gene} eq 'TSG' && $exprn == -1){
      $mutation_classification .= "\t7";
    }
    else{
      $mutation_classification .= "\t0";
    }
    
  }
  
  $output_data .= "\n";
  $mutation_matrix .= "\n";
  $other_matrix .= "\n";
  $mutation_classification .= "\n";
  
}
open OUT, "> $output" or die "Can't write to output file $output: $!\n";
print OUT "$output_header$output_data";
close OUT;

open MUTMAT, "> $output.functional_mutations" or die "Can't write to functional output file $output: $!\n";
print MUTMAT "$matrix_header$mutation_matrix";
close MUTMAT;

open OTHERMAT, "> $output.nonfunctional_mutations" or die "Can't write to non-functional output file $output: $!\n";
print OTHERMAT "$matrix_header$other_matrix";
close OTHERMAT;

open MUTCLASS, "> $output.mutation_classifications.txt" or die "Can't write to mutation classification output file $output: $!\n";
print MUTCLASS "$matrix_header$mutation_classification";
close MUTCLASS;


# In addition to the primary key ("$standard_cell_line\t$standard_gene"),
# the secondary keys for %details_table include:
# cos_mut
# cos_type
# cos_davoli_freq
# icr_mut
# icr_type
# icr_davoli_freq
# biankin_mut
# biankin_type
# biankin_davoli_freq
# wtsi_mut
# wtsi_type
# wtsi_davoli_freq
# ccle_cna
# exprnz

open DETAILS, "> $output.details" or die "Can't write details to $output:$! \n";

print DETAILS "cell_line\tgene_symbol\tEntrez_geneID\tEnsembl_GeneID\tCOSMIC_mutation\tCOSMIC_consequence\tCOSMIC_Davoli\tICR_mutation\tICR_consequence\tICR_Davoli\tBIAKN_mutation\tBIAKN_consequence\tBIANK_Davoli\tWTSI_mutation\tWTSI_consequence\tWTSI_Davoli\tCCLE_CNA\tExprnZ\n";

my @details_keys = keys %details_table;
foreach my $details_key (@details_keys){
  print DETAILS "$details_key";

# cos_mut
  if(exists $details_table{$details_key}{"cos_mut"}){
    chomp($details_table{$details_key}{'cos_mut'});
    print DETAILS "\t$details_table{$details_key}{'cos_mut'}";
  }
  else{
    print DETAILS "\tNA";
  }
  
# cos_type
  if(exists $details_table{$details_key}{"cos_type"}){
    chomp($details_table{$details_key}{"cos_type"});
    print DETAILS "\t$details_table{$details_key}{'cos_type'}";
  }
  else{
    print DETAILS "\tNA";
  }

# cos_davoli_freq
  if(exists $details_table{$details_key}{"cos_davoli_freq"}){
    chomp($details_table{$details_key}{"cos_davoli_freq"});
    print DETAILS "\t$details_table{$details_key}{'cos_davoli_freq'}";
  }
  else{
    print DETAILS "\tNA";
  }

# icr_mut
  if(exists $details_table{$details_key}{"icr_mut"}){
    chomp($details_table{$details_key}{"icr_mut"});
    
    $details_table{$details_key}{"icr_mut"} =~ s/[\r\n]+//g;
    
    print DETAILS "\t$details_table{$details_key}{'icr_mut'}";
  }
  else{
    print DETAILS "\tNA";
  }

# icr_type
  if(exists $details_table{$details_key}{"icr_type"}){
    chomp($details_table{$details_key}{"icr_type"});
    print DETAILS "\t$details_table{$details_key}{'icr_type'}";
  }
  else{
    print DETAILS "\tNA";
  }

# icr_davoli_freq
  if(exists $details_table{$details_key}{"icr_davoli_freq"}){
    chomp($details_table{$details_key}{"icr_davoli_freq"});
    print DETAILS "\t$details_table{$details_key}{'icr_davoli_freq'}";
  }
  else{
    print DETAILS "\tNA";
  }

# biankin_mut
  if(exists $details_table{$details_key}{"biankin_mut"}){
    chomp($details_table{$details_key}{"biankin_mut"});
    print DETAILS "\t$details_table{$details_key}{'biankin_mut'}";
  }
  else{
    print DETAILS "\tNA";
  }

# biankin_type
  if(exists $details_table{$details_key}{"biankin_type"}){
    chomp($details_table{$details_key}{'biankin_type'});
    print DETAILS "\t$details_table{$details_key}{'biankin_type'}";
  }
  else{
    print DETAILS "\tNA";
  }

# biankin_davoli_freq
  if(exists $details_table{$details_key}{"biankin_davoli_freq"}){
    chomp($details_table{$details_key}{"biankin_davoli_freq"});
    print DETAILS "\t$details_table{$details_key}{'biankin_davoli_freq'}";
  }
  else{
    print DETAILS "\tNA";
  }

# wtsi_mut
  if(exists $details_table{$details_key}{"wtsi_mut"}){
    chomp($details_table{$details_key}{"wtsi_mut"});
    print DETAILS "\t$details_table{$details_key}{'wtsi_mut'}";
  }
  else{
    print DETAILS "\tNA";
  }

# wtsi_type
  if(exists $details_table{$details_key}{"wtsi_type"}){
    chomp($details_table{$details_key}{"wtsi_type"});
    print DETAILS "\t$details_table{$details_key}{'wtsi_type'}";
  }
  else{
    print DETAILS "\tNA";
  }

# wtsi_davoli_freq
  if(exists $details_table{$details_key}{"wtsi_davoli_freq"}){
    chomp($details_table{$details_key}{"wtsi_davoli_freq"});
    print DETAILS "\t$details_table{$details_key}{'wtsi_davoli_freq'}";
  }
  else{
    print DETAILS "\tNA";
  }

# ccle_cna
  if(exists $details_table{$details_key}{"ccle_cna"}){
    chomp($details_table{$details_key}{"ccle_cna"});
    print DETAILS "\t$details_table{$details_key}{'ccle_cna'}";
  }
  else{
    print DETAILS "\tNA";
  }

# exprnz
  if(exists $details_table{$details_key}{"exprnz"}){
    chomp($details_table{$details_key}{"exprnz"});
    print DETAILS "\t$details_table{$details_key}{'exprnz'}";
  }
  else{
    print DETAILS "\tNA";
  }

  print DETAILS "\n";
}

close DETAILS;





sub usage() {
  my $usage =<<END;
#----------------------------#
#  process_mutation_data.pl  #
#----------------------------#

James Campbell (jamesc\@icr.ac.uk)

Usage:
perl process_mutation_data.pl [options]

Options
  --help					Display this message and quit
  --ccle_muts				Path to the CCLE maf file [required]
  --cosmic_muts				Path to the COSMIC mutations file [required]
  --vcf_muts				Path to the ICR mutation data [required]
  --ccle_cna				Path to the gistic CCLE CN calls from cBioPortal[required]
  --wtsi_muts				Path to the WTSI mutation[required]
  --biankin_muts			Path to the Biankin lab mutation data [required]
  --wtsi_expression_z_data	Path to the processed WTSI expression data [required]
  --mut_freqs				Path to the processed Davoli data [required]
  --genes					Path to the gene name dictionary [required]
  --cell_lines				Path to the cell line name dictionary [required]
  --output					Path to output. Defaults to ccle_data.proc [optional]
  --output_genes_list		Path to the list of genes requested [required]
END

  print $usage;
}








