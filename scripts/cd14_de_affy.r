# test
library(affy)
library(limma)
library(BiocInstaller)
library(annotate)
# biocLite("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)

setwd("/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD14")

phen <- read.csv("CD14.Caucasian.txt",header=TRUE) # Import Caucasian data, file list made with greping and appending to header

if (file.exists("cd14_cau_raw.RData")) {
	load("cd14_cau_raw.RData")
} else {
	# Memory-intensive step, so better to save after reading in
	cau_data <- ReadAffy(filenames=as.character(phen$FileName),phenoData = data.frame(phen))
	save(list=c("cau_data"),file="cd14_cau_raw.RData")
}

if (file.exists("cd14_cau_expr.RData")) {
	load("cd14_cau_expr.RData")
} else {
	# Save normalized, corrected expression data
	eset <- expresso(cau_data,bgcorrect.method="rma",normalize.method="quantiles",pmcorrect.method="pmonly",summary.method="medianpolish")
	save(list=c("eset"),file="cd14_cau_expr.RData")
}

run_limma <- function(affyset,expset) {
# Takes in affy objects, runs sex-based DE analysis

	samples <- affyset$Sex
	samples <- as.factor(samples)
	design <- model.matrix(~samples)
	# design <- model.matrix(~0)
	colnames(design) <- c("Male","Female")

	fit <- lmFit(exprs(expset),design)
	# fit <- lmFit(exprs(expset))
	contrast.matrix <- makeContrasts(c("Female-Male"),levels=design)
	contrast_fit <- contrasts.fit(fit,contrast.matrix)

	# return(eBayes(fit,robust=TRUE))
	return(eBayes(contrast_fit,robust=TRUE))
}

if (file.exists("cd14_all_limma.RData")) {
	load("cd14_all_limma.RData")

} else {
	ebFit <- run_limma(cau_data,eset)
	save(list=c("ebFit"),file="cd14_all_limma.RData")
}

write.exprs(eset,file="cd14_expr.txt") # For use in PEER, other applications
write.table(pData(eset),file="cd14_pdata.txt")

source("cd14_peer.r")

if (file.exists("cd14_res_limma.RData")) {
	load("cd14_res_limma.RData")

} else {
	resid_exp <- as.matrix(read.table("cd14_cau_residuals.txt",header=TRUE,row.names=1,as.is=TRUE))
	res_eset <- ExpressionSet(assayData=resid_exp)
	
	eb_res_Fit <- run_limma(cau_data,res_eset)

	save(list=c("eb_res_Fit"),file="cd14_res_limma.RData")
}

resid_row_names <- read.table('expr_to_peer.txt',header=TRUE,as.is=TRUE)
sample_names <- eset[[2]]

top_hits <- topTable(ebFit,number=10000)
affy_ids <- rownames(top_hits)

get_gene <- function(id) {hugene10sttranscriptclusterGENENAME[[as.character(id)]]}
genes <- apply(as.matrix(affy_ids),1,get_gene)

res_top_hits <- topTable(eb_res_Fit,number=10000) # Save top (by pvalue) hits
res_affy_ids <- rownames(res_top_hits)

get_res_gene <- function(id) {conv <- resid_row_names;hugene10sttranscriptclusterGENENAME[[as.character(conv[id,1])]]}
res_genes <- apply(as.matrix(res_affy_ids),1,get_res_gene)

write.table(cbind(top_hits,genes),file="cd14_all_results.txt",sep="\t")
write.table(cbind(res_top_hits,res_genes),file="cd14_residual_results.txt",sep="\t")

if (TRUE) {
# Interesting plots
# pdf('cd14_all_pval.pdf')
# hist(ebFit$p.value)
# dev.off()

pdf('cd14_all_volcano.pdf')
volcanoplot(ebFit)
dev.off()

pdf('cd14_res_volcano.pdf')
volcanoplot(eb_res_Fit)
dev.off()

pdf('cd14_residual_pval.pdf')
hist(res_top_hits$P.Value)
dev.off()

pdf('cd14_top_pval.pdf')
hist(top_hits$P.Value)
dev.off()
}
