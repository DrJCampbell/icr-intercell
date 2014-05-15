#!/usr/bin/perl -w
# ============================================= #
# process_CCLE_mutations.pl
# read data from CCLE, count mutations by type
# for each gene and cell line
# jamesc@icr.ac.uk, 14th May 2014
# ============================================= #

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

# read cell line name dictionary

# read gene name dictionary

# ====================================================== #
# Want a table with one row per cell line and one column
# per (gene*consequence). The consequences could be like
# truncating (fs,*,ess) / recurrent missense / other
# 
# How can we represent more complex data such as position
# of mutations in a protein and z-scores? As a binned 
# image, scaled to protein length.
# ====================================================== #

# define globals to store cellline+gene => count for each of:
# truncating
# recurrent_missense
# other


# Open and read mutations file

# get gene, cell line, prot. mut., consequence etc.

# lookup gene and cell line name in dictionaries

# if mutation is truncating, increment cellline+gene truncating count

# if mutation is missense, check for existence in the davoli hash
# if in the davoli hash then inrement recurrent missense count, 
# else increment the 'other' count


# when finished reading file, close

# get list of cell lines and genes (from hash keys)
# foreach cell lines
#   foreach gene
#     add to table

# write out table

























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







