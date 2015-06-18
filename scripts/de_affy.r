library(affy)
library(limma)
library(BiocInstaller)
library(annotate)
biocLite("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)

setwd("/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD4")

phen <- read.csv("CD4.Samples.Caucasian.txt",header=TRUE) # Import Caucasian data, file list made with greping and appending to header
malephen <- which(phen$Sex == "Male")
femalephen <- which(phen$Sex == "Female")

peer_factors <- as.matrix(read.table("cd4_cau_factors.txt",header=TRUE,row.names=1,as.is=TRUE))
peer_residuals <- as.matrix(read.table("cd4_cau_residuals.txt",header=TRUE,row.names=1,as.is=TRUE))
res_eset <- ExpressionSet(assayData=peer_residuals)

if (file.exists("r_affy_data.RData")) {
	load("r_affy_data.RData")
} else {
	# Memory-intensive step, so better to save after reading in
	cau_data <- ReadAffy(filenames=as.character(phen$File.Name),phenoData = data.frame(phen))
	save(list=c("cau_data"),file="r_affy_data.RData")
}

if (file.exists("r_affy_expr.RData")) {
	load("r_affy_expr.RData")
} else {
	# Save normalized, corrected expression data
	eset <- expresso(cau_data,bgcorrect.method="rma",normalize.method="quantiles",pmcorrect.method="pmonly",summary.method="medianpolish")
	save(list=c("eset"),file="r_affy_expr.RData")
}

run_limma <- function(affyset,expset) {
# Takes in affy objects, runs sex-based DE analysis

	samples <- affyset$Sex
	samples <- as.factor(samples)
	# design <- model.matrix(~0 + samples)
	# design <- model.matrix(~samples + resexpset)
	design <- model.matrix(~samples)
	colnames(design) <- c("Male","Female")

	fit <- lmFit(exprs(expset),design)
	# fit <- lmFit(exprs(expset))
	contrast.matrix <- makeContrasts(c("Female-Male"),levels=design)
	contrast_fit <- contrasts.fit(fit,contrast.matrix)

	return(eBayes(contrast_fit,robust=TRUE))
	# return(eBayes(fit,robust=TRUE))
}

if (file.exists("cd4_all_limma.RData")) {
	load("cd4_all_limma.RData")

} else {
	 ebFit <- run_limma(cau_data,eset)
	 save(list=c("ebFit"),file="cd4_all_limma.RData")
}

# write.exprs(eset,file="cd4_expr.txt") # For use in PEER, other applications
# write.table(pData(eset),file="cd4_pdata.txt")

# if (file.exists("cd4_res_limma.RData")) {
	# load("cd4_res_limma.RData")

# } else {
	resid_exp <- as.matrix(read.table("cd4_cau_residuals.txt",header=TRUE,row.names=1,as.is=TRUE))
	res_eset <- ExpressionSet(assayData=resid_exp)
	
	alt_fit <- lm(exprs(eset) ~ peer_factors)
	resid <- residuals(alt_fit)

	samples <- cau_data$Sex
	samples <- as.factor(samples)
	# One column, sex binary
	design <- model.matrix(~0 + samples)
	# design <- model.matrix(~samples)
	colnames(design) <- c("Male","Female")
	head(design)	
	# exp_dif <- exprs(eset) - exprs(res_eset)
	fit <- lmFit(resid,design)
	# fit <- lmFit(exprs(expset))
	contrast.matrix <- makeContrasts(mf = Male - Female,levels=design)
	head(contrast.matrix)
	contrast_fit <- contrasts.fit(fit,contrast.matrix)
	# eb_res_Fit <- eBayes(contrast_fit,robust=FALSE)
	eb_res_Fit <- eBayes(contrast_fit)
	# eb_res_Fit <- run_limma(cau_data,eset,res_eset)

	save(list=c("eb_res_Fit"),file="cd4_res_limma.RData")
# }

resid_row_names <- read.table('expr_to_peer.txt',header=TRUE,as.is=TRUE)
sample_names <- eset[[2]]

top_hits <- topTable(ebFit,number=10000)
affy_ids <- rownames(top_hits)

get_gene <- function(id) {hugene10sttranscriptclusterGENENAME[[as.character(id)]]}
genes <- apply(as.matrix(affy_ids),1,get_gene)

res_top_hits <- topTable(eb_res_Fit,number=10000) # Save top (by pvalue) hits
head(res_top_hits)
res_affy_ids <- rownames(res_top_hits)
# write.table(res_top_hits,file="cd4_residual_results.txt",sep="\t")

x <- hugene10sttranscriptclusterGENENAME
mapped_probes <- mappedkeys(x)
xx <<- as.list(x[mapped_probes])

# get_res_gene <- function(id) {conv <- resid_row_names;xx[[(conv[id,1])]]}
get_res_gene <- function(id) {conv <- resid_row_names;hugene10sttranscriptclusterGENENAME[[as.character(conv[id,1])]]}
# get_gene_id <- function(id) {conv <- resid_row_names;hugene10sttranscriptclusterENTREZID[[as.character(conv[id,1])]]}
# head(xx[as.vector(res_affy_ids,mode='integer')])
# res_genes <- apply(as.matrix(res_affy_ids),1,get_res_gene)
# res_gene_id <- apply(as.matrix(res_affy_ids),1,get_gene_id)
# res_genes 
# write.table(res_gene_id,file='top_res_gene_ids.txt')

# write.table(cbind(top_hits,genes),file="cd4_all_results.txt",sep="\t")
# write.table(cbind(res_top_hits,res_genes),file="cd4_residual_results.txt",sep="\t")
write.table(res_top_hits,file="cd4_residual_results.txt",sep="\t")

# exp_diff <- mean(resid_exp[,malephen]) - mean(resid_exp[,femalephen])
# pdf('sex_diff.pdf')
# hist(exp_diff)
# dev.off()

if (TRUE) {
# Interesting plots
# pdf('pval_dist.pdf')
# hist(ebFit$p.value)
# dev.off()

# pdf('cd4_all_volcano.pdf')
# volcanoplot(ebFit)
# dev.off()

pdf('cd4_res_volcano.pdf')
volcanoplot(eb_res_Fit)
dev.off()

pdf('residual_pval.pdf')
hist(res_top_hits$P.Value)
dev.off()

# pdf('top_pval.pdf')
# hist(top_hits$P.Value)
# dev.off()
