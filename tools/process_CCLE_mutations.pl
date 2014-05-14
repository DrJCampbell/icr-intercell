#!/usr/bin/perl -w
# ============================================= #
# process_CCLE_mutations.pl
# read data from CCLE, count mutations by type
# for each gene and cell line
# jamesc@icr.ac.uk, 14th May 2014
# ============================================= #

use strict;
use Getopt::Long;

# configs you may need to change in the future:
my $base_url = 'http://www.cbioportal.org/public-portal/webservice.do?';

my ($gene_file, $output, $help) = undef;
GetOptions (
  "gene_file=s" => \$gene_file,
  "output=s" => \$output,
  "help" => \$help,
 );

# print usage message if requested or no args supplied
if(defined($help)) {
  &usage;
  exit(0);
}




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
--cell_lines	Path to the cell line name dictionary [required]
--output		Path to output. Defaults to ccle_data.proc [optional]
END

  print $usage;
}







