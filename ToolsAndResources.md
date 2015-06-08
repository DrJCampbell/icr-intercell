# Introduction #

This project contains code and resources used in the analysis of Intercell.


# resources #

resources are stored under /resources.

**cell\_line\_name\_dictionary.txt**
Used to ensure that cell line identifiers from different sources are standardised to those used by the CCLE or follow the CCLE format. The header of the file provides information on how the dictionary was generated. The first column contains cell line name alia. The second column contains the standardised (CCLE) name. The third column contains the source of the alia for reference purposes.

```
# icr-intercell cell line name dictionary v0.1
#
# jamesc@icr.ac.uk, 13th May 2014
#
# This file is used as a dictionary to standardise cell line names in different data sets
# Cell line names found in the first column (cell_line_name) should be converted to those
# given in the second column (intercell_name).
#
# cell line name sources used:
#	COSMIC CLP (http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/download
 -> CosmicCLP_CompleteExport_v68.tsv.gz)
#	CCLE (http://www.broadinstitute.org/ccle/data/browseData
 -> registration required -> CCLE_sample_info_file_2012-10-18.txt)
#	Intercell (http://www.ncbi.nlm.nih.gov/pubmed/21984977)
# The cell line names used in the Intercell analysis attempt to follow the CCLE naming convention.
# the convention is the cell line name (containing only numbers and upper case letters)
# followed by an underscore, followed by the tissue in upper case.
# Where a cell line was not present in the CCLE data, a best effort was made to use
# the tissue type as given by the CCLE.
# This version contains 1537 distinct cell lines names (intercell_name)
#
#cell_line_name	intercell_name	source
22RV1	22RV1_PROSTATE	cosmic_clp
22RV1	22RV1_PROSTATE	gdsc_en_input
22RV1_PROSTATE	22RV1_PROSTATE	CCLE
23132-87	2313287_STOMACH	cosmic_clp
23132-87	2313287_STOMACH	gdsc_en_input
2313287_STOMACH	2313287_STOMACH	CCLE
```

**gene\_name\_dictionary.txt**
Used to standardise gene names used in different data sources. The official gene names and alia are provided by NCBI.

```
# icr-intercell gene name dictionary v0.1
#
# jamesc@icr.ac.uk, 13th May 2014
#
# This file is used as a dictionary to standardise gene names in different data sets
# Cell line names found in the first column (alias) should be converted to those
# given in the second column (gene_name).
#
# Gene name information was extracted from columns 4 and 7 in the file:
# 	ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
#
# This version contains 102615 alia for 47634 distinct gene names
#
#alias	gene_name
DDC8	TIMP2
3-OST-5	HS3ST5
M7V1	SLC1A5
RPS11	RPS11
ECGF-beta	FGF1
RPL10P14	RPL10P14
PNMA1	PNMA1
ZFP318	ZNF318
AT	GATM
```


# tools #

**process\_CCLE\_mutations.pl**
```
usage
perl process_CCLE_mutations.pl [options]

Options
--help                  Display this message and quit
--ccle_data             Path to the CCLE maf file [required]
--mut_freqs             Path to the processed Davoli data [required]
--genes                 Path to the gene name dictionary [required]
--cell_lines            Path to the cell line name dictionary [required]
--output                Path to output. Defaults to ccle_data.proc [optional]
```