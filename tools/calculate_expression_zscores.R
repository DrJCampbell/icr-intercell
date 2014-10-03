#
#
#

require(mixtools)

#setwd("/Users/jamesc/Documents/113_data_sources/Sanger_cell_line_genetic_and_expression_data")

setwd("/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/resources/expression_copy")

# ================================ #
# read in and prepare the data set
# ================================ #


exprn <- read.table(
	file="en_input_w5.cancer_gene_subset.tsv",
	sep="\t",
	header=TRUE,
	row.names=1
	)

tissues <- NULL
i <- NULL
for(i in 1:ncol(exprn)){
	tissues[i] <- sub(
		"^[^_]+_",
		"",
		colnames(exprn)[i],
		perl=TRUE
		)	
}
names(tissues) <- colnames(exprn)
tissue_types <- levels(as.factor(tissues))
tissue_cols <- sample(rainbow(30))


mut_classes <- read.table(
	file="/Users/jamesc/Dropbox/Intercell_II_analysis/icr-intercell/data_tables/combined_mutation_data_v7_excluding_CCLE_hybcap.txt.mutation_classifications.txt",
	sep="\t",
	header=TRUE,
	row.names=1
	)
mut_classes.t <- t(mut_classes)

common.celllines <- intersect(
	colnames(exprn),
	colnames(mut_classes.t)
	)

# get the common cell lines for mut_classes, exprn and tissues
mut_classes.cmn <- NULL
exprn.cmn <- NULL
#exprn.z.cmn <- NULL
tissues.cmn <- NULL
cell_line <- NULL
for(cell_line in common.celllines){
	mut_classes.cmn <- cbind(
		mut_classes.cmn,
		mut_classes.t[,cell_line]
		)
	exprn.cmn <- cbind(
		exprn.cmn,
		exprn[,cell_line]
		)
	tissues.cmn <- c(
		tissues.cmn,
		tissues[cell_line]
		)
#	exprn.z.cmn <- cbind(
#		exprn.z.cmn,
#		exprn.z[,cell_line]
#		)
}
colnames(exprn.cmn) <- common.celllines
#colnames(exprn.z.cmn) <- common.celllines
colnames(mut_classes.cmn) <- common.celllines
rownames(exprn.cmn) <- rownames(exprn)
#rownames(exprn.z.cmn) <- rownames(exprn)
rownames(mut_classes.cmn) <- rownames(mut_classes.t)



# In the input file, rows 2:13,322 have expression data
# 

#exprn.z <- NULL
#for(i in 1:nrow(exprn)){
#	this.mad <- mad(log2(as.numeric(exprn[i,])))
#	this.med <- median(log2(as.numeric(exprn[i,])))
#	this.z <- (log2(as.numeric(exprn[i,])) - this.med) / this.mad
#	exprn.z <- rbind(
#		exprn.z,
#		this.z
#		)
#}
#colnames(exprn.z) <- colnames(exprn)
#rownames(exprn.z) <- rownames(exprn)

#write.table(
#	exprn.z,
#	file="en_input_w5.cancer_gene_subset.log2.zscores_140930.tsv",
#	sep="\t",
#	col.names=TRUE,
#	row.names=TRUE,
#	quote=FALSE
#	)

#exprn.z <- read.table(
#	file="en_input_w5.cancer_gene_subset.log2.zscores_140930.tsv",
#	sep="\t",
#	header=TRUE,
#	row.names=1
#	)

#pdf(file="expression_zscores_140930.pdf")
#boxplot_exprn_z("ERBB2_2064_ENSG00000141736", exprn.z)
#boxplot_exprn_z("PIK3CA_5290_ENSG00000121879", exprn.z)
#boxplot_exprn_z("CDH1_999_ENSG00000039068", exprn.z)
#boxplot_exprn_z("RB1_5925_ENSG00000139687", exprn.z, upper=5, lower=-10)
#boxplot_exprn_z("ARID1A_8289_ENSG00000117713", exprn.z, lower=-8)
#boxplot_exprn_z("CDKN2A_1029_ENSG00000147889", exprn.z)
#boxplot_exprn_z("ZMYM3_9203_ENSG00000147130", exprn.z) # 7.197712e-01
#boxplot_exprn_z("ZNF318_24149_ENSG00000171467", exprn.z) # 6.150296e-12
#boxplot_exprn_z("VHL_7428_ENSG00000134086", exprn.z) # 1.523762e-16
#dev.off()

# foreach gene, test if the distribution is approximately normal
#shapiro_results <- NULL
#i <- NULL
#for(i in 1:nrow(exprn.z)){
#	this.test <- shapiro.test(
#		exprn.z[i,]
#		)
#	shapiro_results <- c(shapiro_results ,this.test$p.value)
#	names(shapiro_results)[i] <- rownames(exprn.z)[i]
#}

#pdf("normalmixEM_plots.pdf")
#mix_results <- NULL
#i <- NULL
#for(i in 1:nrow(exprn.z)){
#	try(
#		this.mix <- normalmixEM(
#			exprn.z[i,], k=2
#			)
#		)
#	mix_results <- rbind(
#		mix_results,
#		c(this.mix$mu[1], this.mix$sigma[1],this.mix$mu[2], this.mix$sigma[2],this.mix$loglik)
#		)
#	rownames(mix_results)[i] <- rownames(exprn.z)[i]
#	try(
#		plot(this.mix,2)
#		)
#}
#colnames(mix_results) <- c("mu1", "sigma1", "mu2", "sigma2", "loglik")
#dev.off()

#cdh1.mix <- normalmixEM(
#			exprn.z[86,], k=2
#			)
#
#str(cdh1.mix)
#plot(cdh1.mix,2, breaks=100)


# ========================================= #
# plot expression values across the tissues #
# old function using plot, not stripchart   #
# ========================================= #

# define colours to plot each tissue
tissue_cols <- sample(rainbow(30))
names(tissue_cols) <- tissue_types

pdf("plot_expression_by_tissue.pdf", width=10, height=5)

plot_exprn_by_tissue(
	exprn,
	gene="CDH1_999_ENSG00000039068",
	tissues,
	tissue_cols=tissue_cols
	)

plot_exprn_by_tissue(
	exprn,
	gene="ERBB2_2064_ENSG00000141736",
	tissues,
	tissue_cols=tissue_cols
	)

plot_exprn_by_tissue(
	exprn,
	gene="RB1_5925_ENSG00000139687",
	tissues,
	tissue_cols=tissue_cols
	)

plot_exprn_by_tissue(
	exprn,
	gene="ARID1A_8289_ENSG00000117713",
	tissues,
	tissue_cols=tissue_cols
	)

dev.off()



# ============================================== #
# plot expression values across the tissues - v2 #
# ============================================== #

# define colours to plot each tissue
tissue_cols <- sample(rainbow(30))
names(tissue_cols) <- tissue_types

pdf("plot_expression_by_tissue_v2.pdf", width=10, height=5)

plot_exprn_by_tissue2(
	exprn,
	gene="CDH1_999_ENSG00000039068",
	tissues,
	tissue_cols=tissue_cols
	)

plot_exprn_by_tissue2(
	exprn,
	gene="ERBB2_2064_ENSG00000141736",
	tissues,
	tissue_cols=tissue_cols
	)

plot_exprn_by_tissue2(
	exprn,
	gene="RB1_5925_ENSG00000139687",
	tissues,
	tissue_cols=tissue_cols
	)

plot_exprn_by_tissue2(
	exprn,
	gene="ARID1A_8289_ENSG00000117713",
	tissues,
	tissue_cols=tissue_cols
	)

dev.off()




# ==================================================== #
# focus on CDH1 in breast. Divide cell line expression
# according to known mutations and wt groups
# ==================================================== #

boxplot(
	log2(exprn.cmn["CDH1_999_ENSG00000039068",which(tissues.cmn == "BREAST")]) ~ mut_classes.cmn["CDH1_999_ENSG00000039068",which(tissues.cmn == "BREAST")],
	pch="",
	names=c(
		"wt",
		"truncation (het)",
		"deletion"
		),
	ylab="log2 normalised CDH1 expression",
	main="WTSI expression data\nfor breast cell lines"
	)

stripchart(
	log2(exprn.cmn["CDH1_999_ENSG00000039068",which(tissues.cmn == "BREAST")])~
	mut_classes.cmn["CDH1_999_ENSG00000039068",which(tissues.cmn == "BREAST")],
	pch=19,
	col=rgb(0,0,0,0.3),
	vertical=TRUE,
	method="jitter",
	jitter=0.2,
	add=TRUE
	)

# print out expression values for:
# wt lines where log2(CHD1) ≤ 8
# trunc lines
# del line

# get wt below 8
log2(exprn.cmn["CDH1_999_ENSG00000039068",which(
	tissues.cmn == "BREAST" &
	mut_classes.cmn["CDH1_999_ENSG00000039068",] == 0 &
	log2(exprn.cmn["CDH1_999_ENSG00000039068",]) <= 8
	)])

# get truncs
log2(exprn.cmn["CDH1_999_ENSG00000039068",which(
	tissues.cmn == "BREAST" &
	mut_classes.cmn["CDH1_999_ENSG00000039068",] == 2
	)])

# get dels
log2(exprn.cmn["CDH1_999_ENSG00000039068",which(
	tissues.cmn == "BREAST" &
	mut_classes.cmn["CDH1_999_ENSG00000039068",] == 5
	)])
colnames(exprn.cmn)[which( # needed because we don't get the cell line name for dels!?!
	tissues.cmn == "BREAST" &
	mut_classes.cmn["CDH1_999_ENSG00000039068",] == 5
	)]





# ------------------------------------------------------------------- #
# Calculate within-tissue expression z-scores for all genes * tissues
# ------------------------------------------------------------------- #

exprn.z.bytissue.cmn <- NULL
for(i in 1:nrow(exprn.cmn)){
	
	this.tissue <- NULL
	
	z.row <- NULL
	
	for(this.tissue in tissue_types){

		this.mad <- mad(log2(as.numeric(exprn.cmn[i,which(tissues.cmn == this.tissue)])))
		
		this.med <- median(log2(as.numeric(exprn.cmn[i,which(tissues.cmn == this.tissue)])))
		
#		print(this.mad)
#		print(this.med)
		
		# if mad is zero, need to set z to NA
		if(this.mad > 0 & !is.na(this.mad) & !is.na(this.med)){
			this.z <- (
				log2(
					as.numeric(
						exprn.cmn[i,which(tissues.cmn == this.tissue)]
						)
					) - this.med
				) / this.mad			
		}else{
			this.z <- rep(NA, length(exprn.cmn[i,which(tissues.cmn == this.tissue)]))
		}
		
		z.row <- c(
			z.row,
			this.z
			)
		
	}
	exprn.z.bytissue.cmn <- rbind(
		exprn.z.bytissue.cmn,
		z.row
		)

}

colnames(exprn.z.bytissue.cmn) <- (common.celllines)
rownames(exprn.z.bytissue.cmn) <- rownames(exprn)

# ================================================================================= #
# reshape the data to long format and write out for use by process_mutation_data.pl
# ================================================================================= #

exprn.z.bytissue.cmn.wide <- as.data.frame(t(exprn.z.bytissue.cmn))

#exprn.z.bytissue.cmn.wide <- cbind(
#	rownames(exprn.z.bytissue.cmn),
#	exprn.z.bytissue.cmn
#	)
#colnames(exprn.z.bytissue.cmn.wide) <- c(
#	"cell.line",
#	colnames(exprn.z.bytissue.cmn)
#	)

exprn.z.bytissue.cmn.long <- reshape(
	exprn.z.bytissue.cmn.wide,
	idvar="cell.line",
	ids=rownames(exprn.z.bytissue.cmn.wide),
	times=names(exprn.z.bytissue.cmn.wide),
	timevar="gene",
	varying=list(names(exprn.z.bytissue.cmn.wide)),
	direction = "long"
	)
colnames(exprn.z.bytissue.cmn.long) <- c(
	"gene",
	"expression.z",
	"cell.line"
	)


write.table(
	exprn.z.bytissue.cmn.long,
	file="expression_zscores_within_tissues_141001.txt",
	sep="\t",
	col.names=TRUE,
	row.names=FALSE,
	quote=FALSE
	)



# ====================================================== #
# box plot expression z within tissue with mutation info
# ====================================================== #

pdf("boxplots_within_tissue_expression_zscores_all_genes_and_tissues.pdf", width=4.5, height=5)

# == ERBB2 == #
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="ERBB2_2064_ENSG00000141736",
	tissue="BREAST",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="ERBB2_2064_ENSG00000141736",
	tissue="OESOPHAGUS",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="ERBB2_2064_ENSG00000141736",
	tissue="LUNG",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="ERBB2_2064_ENSG00000141736",
	tissue="STOMACH",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)

# == RB1 == #
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="RB1_5925_ENSG00000139687",
	tissue="BONE",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="RB1_5925_ENSG00000139687",
	tissue="BREAST",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="RB1_5925_ENSG00000139687",
	tissue="HAEMATOPOIETIC_AND_LYMPHOID_TISSUE",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="RB1_5925_ENSG00000139687",
	tissue="LUNG",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="RB1_5925_ENSG00000139687",
	tissue="OVARY",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)

# == ARID1A == #
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="ARID1A_8289_ENSG00000117713",
	tissue="OVARY",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="ARID1A_8289_ENSG00000117713",
	tissue="KIDNEY",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="ARID1A_8289_ENSG00000117713",
	tissue="STOMACH",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)

# == CDH1 == #
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="CDH1_999_ENSG00000039068",
	tissue="BREAST",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="CDH1_999_ENSG00000039068",
	tissue="LARGE_INTESTINE",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="CDH1_999_ENSG00000039068",
	tissue="LUNG",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="CDH1_999_ENSG00000039068",
	tissue="OESOPHAGUS",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="CDH1_999_ENSG00000039068",
	tissue="OVARY",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="CDH1_999_ENSG00000039068",
	tissue="STOMACH",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)
boxplot_exprn_z_by_tissue_and_mutation(
	exprn.z.bytissue.cmn=exprn.z.bytissue.cmn,
	gene="CDH1_999_ENSG00000039068",
	tissue="URINARY_TRACT",
	tissues.cmn=tissues.cmn,
	mut_classes.cmn=mut_classes.cmn
	)

dev.off()












# ==================== #
#      Functions
# ==================== #



boxplot_exprn_z_by_tissue_and_mutation <- function(
	exprn.z.bytissue.cmn,
	gene,
	tissue,
	tissues.cmn,
	mut_classes.cmn
	){
	
	gene_symbol <- sub("_.+", "", gene, perl=TRUE)
	
	my.names <- NULL
	my.mutclasses <- mut_classes.cmn[gene,which(tissues.cmn == tissue)]
	if(length(
		which(
			my.mutclasses == 0
			)
		) > 0){
		my.names <- c(my.names, "wt")
	}
	if(length(
		which(
			my.mutclasses == 1
			)
		) > 0){
		my.names <- c(my.names, "missense")
	}
	if(length(
		which(
			my.mutclasses == 2
			)
		) > 0){
		my.names <- c(my.names, "trunc (het)")
	}
	if(length(
		which(
			my.mutclasses == 3
			)
		) > 0){
		my.names <- c(my.names, "trunc (hom)")
	}
	if(length(
		which(
			my.mutclasses == 4
			)
		) > 0){
		my.names <- c(my.names, "amp")
	}
	if(length(
		which(
			my.mutclasses == 5
			)
		) > 0){
		my.names <- c(my.names, "del")
	}
	
	boxplot(
		exprn.z.bytissue.cmn[gene,which(tissues.cmn == tissue)] ~ mut_classes.cmn[gene,which(tissues.cmn == tissue)],
		pch="",
		names=my.names,
		ylab=paste(gene_symbol, " expression z-score", sep=""),
		main=paste("WTSI expression z-scores\n", tissue, sep="")
		)
	
	stripchart(
		exprn.z.bytissue.cmn[gene,which(tissues.cmn == tissue)]~
		mut_classes.cmn[gene,which(tissues.cmn == tissue)],
		pch=19,
		col=rgb(0,0,0,0.3),
		vertical=TRUE,
		method="jitter",
		jitter=0.2,
		add=TRUE
		)
	
	abline(
		-2,0,lty=2,col="red"
		)
	abline(
		2,0,lty=2,col="red"
		)
	
}




boxplot_exprn_z <- function(
	gene,
	zscores,
	upper = 5,
	lower=-5
	){

	xvals <- jitter(rep(1, length(zscores[gene,])), amount=0.15)
	overexprn <- which(as.numeric(zscores[gene,]) > 2)
	underexprn <- which(as.numeric(zscores[gene,]) < -2)
	
	boxplot(
		as.numeric(zscores[gene,]),
		pch="",
		ylab="log2 expression z-score",
		main=paste(gene, "expression"),
		ylim=c(lower,upper)
		)
	points(
		xvals,
		as.numeric(zscores[gene,]),
		pch=19,
		col=rgb(0,0,0,0.25)
		)
	abline(
		-2,
		0,
		col="red",
		lty=2
		)
	abline(
		2,
		0,
		col="red",
		lty=2
		)
	
	ytop <- upper - 0.5
	i <- NULL
	for(i in overexprn){
		text(
			1.25,
			ytop,
			colnames(zscores)[i],
			pos=4
			)
		lines(
			c(1.25,1.25),
			c(2,max(as.numeric(zscores[gene,]))),
			col=rgb(0,0,0,0.25)
			)
		lines(
			c(1.2,1.25),
			c(2,2),
			col=rgb(0,0,0,0.25)
			)
		lines(
			c(1.2,1.25),
			c(max(as.numeric(zscores[gene,])),max(as.numeric(zscores[gene,]))),
			col=rgb(0,0,0,0.25)
			)

		ytop = ytop - 0.3
	}

	ybot <- lower + 0.5
	i <- NULL
	for(i in underexprn){
		text(
			0.75,
			ybot,
			colnames(zscores)[i],
			pos=2
			)
		lines(
			c(0.75,0.75),
			c(-2,min(as.numeric(zscores[gene,]))),
			col=rgb(0,0,0,0.25)
			)
		lines(
			c(0.8,0.75),
			c(-2,-2),
			col=rgb(0,0,0,0.25)
			)
		lines(
			c(0.8,0.75),
			c(min(as.numeric(zscores[gene,])),min(as.numeric(zscores[gene,]))),
			col=rgb(0,0,0,0.25)
			)

		ybot = ybot + 0.3
	}
}

# ============================== #

plot_exprn_by_tissue <- function(
	exprn,
	gene,
	tissues,
	tissue_cols
	){
	tissues <- NULL
	i <- NULL
	for(i in 1:ncol(exprn)){
		tissues[i] <- sub(
			"^[^_]+_",
			"",
			colnames(exprn)[i],
			perl=TRUE
			)	
	}
	
	tissue_types <- levels(as.factor(tissues))
	
#	tissue_cols <- sample(rainbow(30))
#	names(tissue_cols) <- tissue_types
	
	x_positions <- 1:ncol(exprn)
	
	par(mai=c(2,1,0.5,0.5))
	plot(
		NULL,
		NULL,
		xlim=c(1,ncol(exprn)),
		ylim=c(
			min(
				log2(exprn[gene,])
				),
			max(
				log2(exprn[gene,])
				)
			),
		xlab="",
		ylab="log2 normalised expression",
		main=gene,
		xaxt="n",
		bty="n"
		)
	tissue_type <- NULL
	for(tissue_type in tissue_types){
		points(
			x_positions[which(tissues == tissue_type)],
			log2(exprn[gene,which(tissues == tissue_type)]),
			col=tissue_cols[tissue_type]
			)
		axis(
			side=1,
			at=c(min(which(tissues == tissue_type)),max(which(tissues == tissue_type))),
			labels=FALSE,
			lwd=0.75,
			tcl=-0.2
			)
		label_x_pos <- min(
			which(
				tissues == tissue_type
				)
			)+(
			max(
				which(
					tissues == tissue_type
					)
				)-min(
					which(
						tissues == tissue_type
						)
					)
				)/2 # end label_x_pos calculation
		mtext(
			tissue_type,
			side=1,
			at=label_x_pos,
			las=2,
			cex=0.75,
			line=0.5
			)
	}
	
}


# ============================== #

plot_exprn_by_tissue2 <- function(
	exprn,
	gene,
	tissues,
	tissue_cols
	){
	tissues <- NULL
	i <- NULL
	for(i in 1:ncol(exprn)){
		tissues[i] <- sub(
			"^[^_]+_",
			"",
			colnames(exprn)[i],
			perl=TRUE
			)	
	}
	
	tissue_types <- levels(as.factor(tissues))
	
#	tissue_cols <- sample(rainbow(30))
#	names(tissue_cols) <- tissue_types
	
	x_positions <- 1:ncol(exprn)
	
	par(mai=c(2,1,0.5,0.5))
	
	stripchart(
		log2(exprn.cmn[gene,]) ~ tissues.cmn,
		method="jitter",
		vertical=TRUE,
		pch=19,
		col=rgb(0,0,0,0.25),
		jitter=0.2,
		las=2,
		cex.axis=0.75,
		ylab="log2 normalised expression",
		main=gene
		)
	
}
