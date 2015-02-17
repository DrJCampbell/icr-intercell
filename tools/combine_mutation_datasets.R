#
# combine processed mutation data sets
# jamesc@icr.ac.uk, 6th Jan 2015
#

setwd("/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables")


bone_exome_func_muts_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_bone_exome_data_150129.functional_mutations.txt"
bone_exome_mut_classes_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_bone_exome_data_150129.mutation_classifications.txt"
bone_exome_all_muts_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_bone_exome_data_150129.nonfunctional_mutations.txt"

ovarian_exome_func_muts_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_ovarian_exome_data_150129.functional_mutations.txt"
ovarian_exome_mut_classes_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_ovarian_exome_data_150129.mutation_classifications.txt"
ovarian_exome_all_muts_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_ovarian_exome_data_150129.nonfunctional_mutations.txt"

wtsi_exome_func_muts_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_WTSI_exome_data_150129.functional_mutations.txt"
wtsi_exome_mut_classes_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_WTSI_exome_data_150129.mutation_classifications.txt"
wtsi_exome_all_muts_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_WTSI_exome_data_150129.nonfunctional_mutations.txt"

cosmic_exome_func_muts_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_cosmic_exome_150129.functional_mutations.txt"
cosmic_exome_mut_classes_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_cosmic_exome_150129.mutation_classifications.txt"
cosmic_exome_all_muts_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_cosmic_exome_150129.nonfunctional_mutations.txt"

wtsi_cnv_func_muts_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_WTSI_cnv_data_150202.functional_mutations.txt"
wtsi_cnv_mut_classes_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_WTSI_cnv_data_150202.mutation_classifications.txt"
wtsi_cnv_all_muts_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_WTSI_cnv_data_150202.nonfunctional_mutations.txt"

ccle_cnv_func_muts_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_CCLE_cnv_data_150202.functional_mutations.txt"
ccle_cnv_mut_classes_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_CCLE_cnv_data_150202.mutation_classifications.txt"
ccle_cnv_all_muts_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_CCLE_cnv_data_150202.nonfunctional_mutations.txt"


wtsi_exprn_zscores_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_WTSI_exprn_data_150129.expression_zscores.txt"
wtsi_exprn_func_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_WTSI_exprn_data_150129.functional_mutations.txt"
wtsi_exprn_classes_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_WTSI_exprn_data_150129.mutation_classifications.txt"
wtsi_exprn_all_file <- "/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/processed_WTSI_exprn_data_150129.nonfunctional_mutations.txt"

combined_exome_func_muts_file <- "combined_exome_func_muts_150203.txt"
combined_exome_all_muts_file <- "combined_exome_all_muts_150203.txt"
combined_exome_mut_classes_file <- "combined_exome_mut_classes_150203.txt"
combined_cnv_func_muts_file <- "combined_cnv_func_muts_150203.txt"
combined_cnv_all_muts_file <- "combined_cnv_all_muts_150203.txt"
combined_cnv_mut_classes_file <- "combined_cnv_mut_classes_150203.txt"
combined_exome_cnv_func_muts_file <- "combined_exome_cnv_func_muts_150203.txt"
combined_exome_cnv_all_muts_file <- "combined_exome_cnv_all_muts_150203.txt"
combined_exome_cnv_mut_classes_file <- "combined_exome_cnv_mut_classes_150203.txt"
combined_exome_cnv_exprn_func_muts_file <- "combined_exome_cnv_exprn_func_muts_150203.txt"
combined_exome_cnv_exprn_all_muts_file <- "combined_exome_cnv_exprn_all_muts_150203.txt"
combined_exome_cnv_exprn_mut_classes_file <- "combined_exome_cnv_exprn_mut_classes_150203.txt"


bone_exome_func_muts <- read.table(
	bone_exome_func_muts_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

bone_exome_mut_classes <- read.table(
	bone_exome_mut_classes_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

bone_exome_all_muts <- read.table(
	bone_exome_all_muts_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

ovarian_exome_func_muts <- read.table(
	ovarian_exome_func_muts_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

ovarian_exome_mut_classes <- read.table(
	ovarian_exome_mut_classes_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

ovarian_exome_all_muts <- read.table(
	ovarian_exome_all_muts_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

wtsi_exome_func_muts <- read.table(
	wtsi_exome_func_muts_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

wtsi_exome_mut_classes <- read.table(
	wtsi_exome_mut_classes_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

wtsi_exome_all_muts <- read.table(
	wtsi_exome_all_muts_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

cosmic_exome_func_muts <- read.table(
	cosmic_exome_func_muts_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

cosmic_exome_mut_classes <- read.table(
	cosmic_exome_mut_classes_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

cosmic_exome_all_muts <- read.table(
	cosmic_exome_all_muts_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

wtsi_cnv_func_muts <- read.table(
	wtsi_cnv_func_muts_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

wtsi_cnv_mut_classes <- read.table(
	wtsi_cnv_mut_classes_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

wtsi_cnv_all_muts <- read.table(
	wtsi_cnv_all_muts_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

# Need to fix the fact that some cell lines in the WTSI data set will have
# have CNV data but contain no alterations in the set of 572 cancer genes.
# These will be eroneously excluded if we do not set the CNV status to 0.

wtsi_cell_lines <- rownames(wtsi_exome_all_muts)
wtsi_cell_lines_needing_cnvs <- setdiff(wtsi_cell_lines, rownames(wtsi_cnv_all_muts))
wtsi_cnv_block_to_add <- matrix(data=0, nrow=length(wtsi_cell_lines_needing_cnvs), ncol=ncol(wtsi_cnv_all_muts))
colnames(wtsi_cnv_block_to_add) <- colnames(wtsi_cnv_all_muts)
rownames(wtsi_cnv_block_to_add) <- wtsi_cell_lines_needing_cnvs

wtsi_cnv_func_muts <- rbind(
	wtsi_cnv_func_muts,
	wtsi_cnv_block_to_add
	)
wtsi_cnv_mut_classes <- rbind(
	wtsi_cnv_mut_classes,
	wtsi_cnv_block_to_add
	)
wtsi_cnv_all_muts <- rbind(
	wtsi_cnv_all_muts,
	wtsi_cnv_block_to_add
	)



ccle_cnv_func_muts <- read.table(
	ccle_cnv_func_muts_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

ccle_cnv_mut_classes <- read.table(
	ccle_cnv_mut_classes_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

ccle_cnv_all_muts <- read.table(
	ccle_cnv_all_muts_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)


wtsi_exprn_func_muts <- read.table(
	wtsi_exprn_func_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

wtsi_exprn_all_muts <- read.table(
	wtsi_exprn_all_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)

wtsi_exprn_mut_classes <- read.table(
	wtsi_exprn_classes_file,
	sep="\t",
	header=TRUE,
	row.names=1
	)


# ======================================= #
# Join all exome data sets
# 1. find union of cell line names
# 2. for each cell line, OR together the
#    mutation status values for each gene
# ======================================= #

common_exome_cell_lines <- 
	union(
		rownames(bone_exome_func_muts),
		union(
			rownames(ovarian_exome_func_muts),
			union(
				rownames(wtsi_exome_func_muts),
				rownames(cosmic_exome_func_muts)
				)
			)
		)
		# note that most cell lines are covered by the COSMIC data
		# all 97 of the WTSI data set (exomes) plus five of the 11
		# ovarian cell lines and nine of the 17 bone lines we have
		# this leaves 1023 individual cell lines

common_exome_genes <- 
	union(
		colnames(bone_exome_func_muts),
		union(
			colnames(ovarian_exome_func_muts),
			union(
				colnames(wtsi_exome_func_muts),
				colnames(cosmic_exome_func_muts)
				)
			)
		)
		# all 572 genes should be included

# set up matrices to hold the joined data sets
combined_exome_func_muts <- matrix(
	data=NA,
	nrow=length(common_exome_cell_lines),
	ncol=length(common_exome_genes),
	dimnames=list(
		common_exome_cell_lines,
		common_exome_genes
		)
	)
combined_exome_all_muts <- matrix(
	data=NA,
	nrow=length(common_exome_cell_lines),
	ncol=length(common_exome_genes),
	dimnames=list(
		common_exome_cell_lines,
		common_exome_genes
		)
	)
combined_exome_mut_classes <- matrix(
	data=NA,
	nrow=length(common_exome_cell_lines),
	ncol=length(common_exome_genes),
	dimnames=list(
		common_exome_cell_lines,
		common_exome_genes
		)
	)

# loop through each set of matrices and OR the values (take the max)
cell_line <- NULL
for(cell_line in common_exome_cell_lines){
	gene <- NULL
	for(gene in common_exome_genes){
		
		if(!is.na(cosmic_exome_func_muts[cell_line,gene])){
			combined_exome_func_muts[cell_line,gene] <- cosmic_exome_func_muts[cell_line,gene]
		}else{
			combined_exome_func_muts[cell_line,gene] = 0
		}
		if(!is.na(wtsi_exome_func_muts[cell_line,gene])){
			combined_exome_func_muts[cell_line,gene] <- max(
				combined_exome_func_muts[cell_line,gene],
				wtsi_exome_func_muts[cell_line,gene]
				)
		}
		if(!is.na(bone_exome_func_muts[cell_line,gene])){
			combined_exome_func_muts[cell_line,gene] <- max(
				combined_exome_func_muts[cell_line,gene],
				bone_exome_func_muts[cell_line,gene]
				)
		}
		if(!is.na(ovarian_exome_func_muts[cell_line,gene])){
			combined_exome_func_muts[cell_line,gene] <- max(
				combined_exome_func_muts[cell_line,gene],
				ovarian_exome_func_muts[cell_line,gene]
				)
		}

		if(!is.na(cosmic_exome_all_muts[cell_line,gene])){
			combined_exome_all_muts[cell_line,gene] <- cosmic_exome_all_muts[cell_line,gene]
		}else{
			combined_exome_all_muts[cell_line,gene] = 0
		}
		
		if(!is.na(wtsi_exome_all_muts[cell_line,gene])){
			combined_exome_all_muts[cell_line,gene] <- max(
				combined_exome_all_muts[cell_line,gene],
				wtsi_exome_all_muts[cell_line,gene]
				)
		}
		if(!is.na(bone_exome_all_muts[cell_line,gene])){
			combined_exome_all_muts[cell_line,gene] <- max(
				combined_exome_all_muts[cell_line,gene],
				bone_exome_all_muts[cell_line,gene]
				)
		}
		if(!is.na(ovarian_exome_all_muts[cell_line,gene])){
			combined_exome_all_muts[cell_line,gene] <- max(
				combined_exome_all_muts[cell_line,gene],
				ovarian_exome_all_muts[cell_line,gene]
				)
		}


		# do mut classes - for exome, take the max integer value class observed (i.e. truncation tops missense...)
		if(!is.na(cosmic_exome_mut_classes[cell_line,gene])){
			combined_exome_mut_classes[cell_line,gene] <- cosmic_exome_mut_classes[cell_line,gene]
		}else{
			combined_exome_mut_classes[cell_line,gene] = 0
		}
		
		if(!is.na(wtsi_exome_mut_classes[cell_line,gene])){
			combined_exome_mut_classes[cell_line,gene] <- max(combined_exome_mut_classes[cell_line,gene], wtsi_exome_mut_classes[cell_line,gene])
		}
		
		if(!is.na(bone_exome_mut_classes[cell_line,gene])){
			combined_exome_mut_classes[cell_line,gene] <- max(combined_exome_mut_classes[cell_line,gene], bone_exome_mut_classes[cell_line,gene])
		}
		
		if(!is.na(ovarian_exome_mut_classes[cell_line,gene])){
			combined_exome_mut_classes[cell_line,gene] <- max(combined_exome_mut_classes[cell_line,gene], ovarian_exome_mut_classes[cell_line,gene])
		}
		
	}
}



# write out the combined data sets
write.table(
	combined_exome_func_muts,
	file=combined_exome_func_muts_file,
	sep="\t",
	col.names=TRUE,
	row.names=TRUE,
	quote=FALSE
	)

write.table(
	combined_exome_all_muts,
	file=combined_exome_all_muts_file,
	sep="\t",
	col.names=TRUE,
	row.names=TRUE,
	quote=FALSE
	)


write.table(
	combined_exome_mut_classes,
	file=combined_exome_mut_classes_file,
	sep="\t",
	col.names=TRUE,
	row.names=TRUE,
	quote=FALSE
	)

#temp <- apply(combined_exome_func_muts,2,sum)
#names(temp) <- colnames(combined_exome_func_muts)
#head(sort(temp, decreasing=TRUE),50)



# ======================================= #
# Join all CNV data sets (CCLE and WTSI)
# 1. find union of cell line names
# 2. for each cell line, OR together the
#    mutation status values for each gene
# ======================================= #

common_cnv_cell_lines <- 
	union(
		rownames(wtsi_cnv_func_muts),
		rownames(ccle_cnv_func_muts)
		)

common_cnv_genes <- 
	union(
		colnames(wtsi_cnv_func_muts),
		colnames(ccle_cnv_func_muts)
		)
		# all 572 genes should be included

# set up matrices to hold the joined data sets
combined_cnv_func_muts <- matrix(
	data=NA,
	nrow=length(common_cnv_cell_lines),
	ncol=length(common_cnv_genes),
	dimnames=list(
		common_cnv_cell_lines,
		common_cnv_genes
		)
	)
combined_cnv_all_muts <- matrix(
	data=NA,
	nrow=length(common_cnv_cell_lines),
	ncol=length(common_cnv_genes),
	dimnames=list(
		common_cnv_cell_lines,
		common_cnv_genes
		)
	)
combined_cnv_mut_classes <- matrix(
	data=NA,
	nrow=length(common_cnv_cell_lines),
	ncol=length(common_cnv_genes),
	dimnames=list(
		common_cnv_cell_lines,
		common_cnv_genes
		)
	)

# loop through each set of matrices and OR the values (take the max)
cell_line <- NULL
for(cell_line in common_cnv_cell_lines){
	gene <- NULL
	for(gene in common_cnv_genes){
		
		if(!is.na(ccle_cnv_func_muts[cell_line,gene])){
			combined_cnv_func_muts[cell_line,gene] <- ccle_cnv_func_muts[cell_line,gene]
		}else{
			combined_cnv_func_muts[cell_line,gene] = 0
		}
		if(!is.na(wtsi_cnv_func_muts[cell_line,gene])){
			combined_cnv_func_muts[cell_line,gene] <- max(
				combined_cnv_func_muts[cell_line,gene],
				wtsi_cnv_func_muts[cell_line,gene]
				)
		}

		if(!is.na(ccle_cnv_all_muts[cell_line,gene])){
			combined_cnv_all_muts[cell_line,gene] <- ccle_cnv_all_muts[cell_line,gene]
		}else{
			combined_cnv_all_muts[cell_line,gene] = 0
		}
		
		if(!is.na(wtsi_cnv_all_muts[cell_line,gene])){
			combined_cnv_all_muts[cell_line,gene] <- max(
				combined_cnv_all_muts[cell_line,gene],
				wtsi_cnv_all_muts[cell_line,gene]
				)
		}


		# do mut classes - for cnv, take the max integer value class observed (i.e. truncation tops missense...)
		if(!is.na(ccle_cnv_mut_classes[cell_line,gene])){
			combined_cnv_mut_classes[cell_line,gene] <- ccle_cnv_mut_classes[cell_line,gene]
		}else{
			combined_cnv_mut_classes[cell_line,gene] = 0
		}
		
		if(!is.na(wtsi_cnv_mut_classes[cell_line,gene])){
			if(wtsi_cnv_mut_classes[cell_line,gene] == combined_cnv_mut_classes[cell_line,gene]){
				next
			}else{
				combined_cnv_mut_classes[cell_line,gene] <- max(combined_cnv_mut_classes[cell_line,gene], wtsi_cnv_mut_classes[cell_line,gene])
			}
		}
		
	}
}

# write out the combined data sets
write.table(
	combined_cnv_func_muts,
	file=combined_cnv_func_muts_file,
	sep="\t",
	col.names=TRUE,
	row.names=TRUE,
	quote=FALSE
	)

write.table(
	combined_cnv_all_muts,
	file=combined_cnv_all_muts_file,
	sep="\t",
	col.names=TRUE,
	row.names=TRUE,
	quote=FALSE
	)

write.table(
	combined_cnv_mut_classes,
	file=combined_cnv_mut_classes_file,
	sep="\t",
	col.names=TRUE,
	row.names=TRUE,
	quote=FALSE
	)




# =================================== #
# combine the exome and CNV data sets
# =================================== #

common_genes <- intersect(
	colnames(combined_exome_func_muts),
	colnames(combined_cnv_func_muts)
	)

common_cell_lines <- intersect(
	rownames(combined_exome_func_muts),
	rownames(combined_cnv_func_muts)
	)

# loop through the cell lines (rows), then through the genes (cols)
# for the common sets. OR the values together into combined matrices

combined_exome_cnv_func_muts <- matrix(
	data=NA,
	nrow=length(common_cell_lines),
	ncol=length(common_genes),
	dimnames=list(common_cell_lines,common_genes)
	)
combined_exome_cnv_all_muts <- matrix(
	data=NA,
	nrow=length(common_cell_lines),
	ncol=length(common_genes),
	dimnames=list(common_cell_lines,common_genes)
	)
combined_exome_cnv_mut_classes <- matrix(
	data=NA,
	nrow=length(common_cell_lines),
	ncol=length(common_genes),
	dimnames=list(common_cell_lines,common_genes)
	)

cell_line <- NULL
for(cell_line in common_cell_lines){
	gene <- NULL
	for(gene in common_genes){
		
		if(
			combined_exome_func_muts[cell_line,gene] == 1 |
		 	combined_cnv_func_muts[cell_line,gene] == 1){
			combined_exome_cnv_func_muts[cell_line,gene] = 1
		}else{
			combined_exome_cnv_func_muts[cell_line,gene] = 0
		}
		
		if(
			combined_exome_all_muts[cell_line,gene] == 1 |
		 	combined_cnv_all_muts[cell_line,gene] == 1){
			combined_exome_cnv_all_muts[cell_line,gene] = 1
		}else{
			combined_exome_cnv_all_muts[cell_line,gene] = 0
		}

		if(
			combined_exome_mut_classes[cell_line,gene] > 0 |
		 	combined_cnv_mut_classes[cell_line,gene] > 0){
			if(combined_exome_mut_classes[cell_line,gene] > 0){
				combined_exome_cnv_mut_classes[cell_line,gene] = combined_exome_mut_classes[cell_line,gene]
			}else{
				combined_exome_cnv_mut_classes[cell_line,gene] = combined_cnv_mut_classes[cell_line,gene]
			}
		}else{
			combined_exome_cnv_mut_classes[cell_line,gene] = 0
		}
	
	}
}

write.table(
	combined_exome_cnv_func_muts,
	file=combined_exome_cnv_func_muts_file,
	sep="\t",
	col.names=TRUE,
	row.names=TRUE,
	quote=FALSE
	)

write.table(
	combined_exome_cnv_all_muts,
	file=combined_exome_cnv_all_muts_file,
	sep="\t",
	col.names=TRUE,
	row.names=TRUE,
	quote=FALSE
	)

write.table(
	combined_exome_cnv_mut_classes,
	file=combined_exome_cnv_mut_classes_file,
	sep="\t",
	col.names=TRUE,
	row.names=TRUE,
	quote=FALSE
	)


# ==================================== #
# combine the expression data set with
# the exome and copy number data
# ==================================== #

common_genes_exprn <- intersect(
	common_genes,
	colnames(wtsi_exprn_func_muts)
	)

common_cell_lines_exprn <- intersect(
	common_cell_lines,
	rownames(wtsi_exprn_func_muts)
	)

# loop through the cell lines (rows), then through the genes (cols)
# for the common sets. OR the values together into combined matrices

combined_exome_cnv_exprn_func_muts <- matrix(
	data=NA,
	nrow=length(common_cell_lines_exprn),
	ncol=length(common_genes_exprn),
	dimnames=list(common_cell_lines_exprn,common_genes_exprn)
	)
combined_exome_cnv_exprn_all_muts <- matrix(
	data=NA,
	nrow=length(common_cell_lines_exprn),
	ncol=length(common_genes_exprn),
	dimnames=list(common_cell_lines_exprn,common_genes_exprn)
	)
combined_exome_cnv_exprn_mut_classes <- matrix(
	data=NA,
	nrow=length(common_cell_lines_exprn),
	ncol=length(common_genes_exprn),
	dimnames=list(common_cell_lines_exprn,common_genes_exprn)
	)

cell_line <- NULL
for(cell_line in common_cell_lines_exprn){
	gene <- NULL
	for(gene in common_genes_exprn){
		
		if(
			combined_exome_cnv_func_muts[cell_line,gene] == 1 |
		 	wtsi_exprn_func_muts[cell_line,gene] == 1){
			combined_exome_cnv_exprn_func_muts[cell_line,gene] = 1
		}else{
			combined_exome_cnv_exprn_func_muts[cell_line,gene] = 0
		}
		
		if(
			combined_exome_cnv_all_muts[cell_line,gene] == 1 |
		 	wtsi_exprn_all_muts[cell_line,gene] == 1){
			combined_exome_cnv_exprn_all_muts[cell_line,gene] = 1
		}else{
			combined_exome_cnv_exprn_all_muts[cell_line,gene] = 0
		}

		if(
			combined_exome_cnv_mut_classes[cell_line,gene] > 0 |
		 	wtsi_exprn_mut_classes[cell_line,gene] > 0){
			if(combined_exome_cnv_mut_classes[cell_line,gene] > 0){
				combined_exome_cnv_exprn_mut_classes[cell_line,gene] = combined_exome_cnv_mut_classes[cell_line,gene]
			}else{
				combined_exome_cnv_exprn_mut_classes[cell_line,gene] = wtsi_exprn_mut_classes[cell_line,gene]
			}
		}else{
			combined_exome_cnv_exprn_mut_classes[cell_line,gene] = 0
		}
	
	}
}



write.table(
	combined_exome_cnv_exprn_func_muts,
	file=combined_exome_cnv_exprn_func_muts_file,
	sep="\t",
	col.names=TRUE,
	row.names=TRUE,
	quote=FALSE
	)

write.table(
	combined_exome_cnv_exprn_all_muts,
	file=combined_exome_cnv_exprn_all_muts_file,
	sep="\t",
	col.names=TRUE,
	row.names=TRUE,
	quote=FALSE
	)

write.table(
	combined_exome_cnv_exprn_mut_classes,
	file=combined_exome_cnv_exprn_mut_classes_file,
	sep="\t",
	col.names=TRUE,
	row.names=TRUE,
	quote=FALSE
	)


