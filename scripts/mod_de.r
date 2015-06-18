# test
library(affy)
library(ggplot2)
library(limma)
#library(BiocInstaller)
#biocLite("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)

wd = "/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD4"
setwd(wd)

# Rdata file names
#In
raw_exp_rfile <- "r_affy_data.RData"
norm_exp_rfile <- "r_affy_expr.RData"
phenotypes <- "CD4.Samples.Caucasian.txt"
peer_factors_file <- "cd4_cau_factors.txt"
covarites <- "cd4.cau.cov"
lm_peer_residuals_file <- "cd4_factor_residuals.txt"

phen <- read.csv(file=phenotypes,header=TRUE)
covs <- read.csv('cd4.cau.cov',header=TRUE,row.names=1)

if (file.exists(raw_exp_rfile)) load(raw_exp_rfile)
else {
	cau_data <- ReadAffy(filenames=as.character(phen$FileName),phenoData = data.frame(phen))
	save(list=c("cau_data"),file=raw_exp_rfile)
}

if (file.exists(norm_exp_rfile)) load(norm_exp_rfile)
else {
	eset <- expresso(cau_data,bgcorrect.method="rma",
		normalize.method="quantiles",
		pmcorrect.method="pmonly",
		summary.method="medianpolish")
	save(list=c("eset"),file=norm_exp_rfile)
}

sex <<- unlist(as.matrix(covs[,1],nrow=1))
exp <- exprs(eset)
make_residuals <- function(in_row) {
	residuals(lm(in_row ~ 0 + peer_factors[,1] + peer_factors[,2] + peer_factors[,3] + peer_factors[,4] + peer_factors[,5] + peer_factors[,6] + peer_factors[,7] + peer_factors[,8] + peer_factors[,9] + peer_factors[,10] + peer_factors[,11] + peer_factors[,12] + peer_factors[,13] + peer_factors[,14] + peer_factors[,15] + peer_factors[,16] + peer_factors[,17] + peer_factors[,18] + peer_factors[,19] + peer_factors[,20]))
}

if (file.exists(file=lm_peer_residuals_file)) {
	res <- as.matrix(read.table(file=lm_peer_residuals_file,header=TRUE,row.names=1))
	} else {
	peer_factors <<- as.matrix(read.table(file=peer_factors_file,header=TRUE,row.names=1,as.is=TRUE))
	res <- apply(t(exp),1,make_residuals)
	write.table(res,file=lm_peer_residuals_file)
}

make_ebayes_fit <- function(exp,samples) {

	samples <- as.factor(samples)
	design <- model.matrix(~0 + samples)
	colnames(design) <- c("Female","Male")
	fit <- lmFit(exp,design)
	contrast.matrix <- makeContrasts(mf = Male - Female,levels=design)
	contrast_fit <- contrasts.fit(fit,contrast.matrix)
	eb_Fit <- eBayes(contrast_fit,robust=TRUE)
	return(eb_Fit)
}

samples <- cau_data$Sex
n_eb <- make_ebayes_fit(exprs(eset),samples)
p_eb <- make_ebayes_fit(t(res),samples)

########################################## Analysis of results ##################################
sex_expr <- NULL
for (id in row.names(p_eb)) {
	sex_expr <- rbind(sex_expr,data.frame(id=id,male=mean(res[id,!!samples]),female=mean(res[id,!samples])))
}

n_top_hits <<- topTable(n_eb,number=1000,coef=1)
p_top_hits <<- topTable(p_eb,number=1000,coef=1)
write.table(row.names(p_top_hits),file="cd4_peer_topprobes.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

comp <- function(id) {
	a <- n_top_hits[id,]
	b <- p_top_hits[id,]
	
	a_p <- a$adj.P.Value 
	b_p <- b$adj.P.Value

	if (is.null(a_p) || is.null(b_p)) return(1) 
	else return(as.numeric(a_p) - as.numeric(b_p)) }

diffs <- apply(as.matrix(row.names(n_top_hits)),1,comp)

get_gene <- function(id) {
	gene <- hugene10sttranscriptclusterGENENAME[[as.character(id)]]
	chr <- hugene10sttranscriptclusterCHR[[as.character(id)]]
	return(matrix(c(id,gene,chr),nrow=1))
}

n_ids <- as.matrix(row.names(n_top_hits),ncol=1)
p_ids <- as.matrix(row.names(p_top_hits),ncol=1)

n_probeinfo <- apply(n_ids,1,get_gene)
p_probeinfo <- apply(p_ids,1,get_gene)
# Gx1 list. Index each row with [[i]]. Returns 1x3 character array with id, gene_name, and chr as [,1],[,2],and [,3] indicies

if (FALSE) {
n_genes <- NULL
for (i in 1:length(n_probeinfo)) {
	n_genes[i] <- n_probeinfo[[i]][,2]
} 
p_genes <- NULL
for (i in 1:length(p_probeinfo)) {
	p_genes[i] <- p_probeinfo[[i]][,2]
} 
n_chrs <- NULL
for (i in 1:length(n_probeinfo)) {
	n_chrs[i] <- n_probeinfo[[i]][,3]
} 
p_chrs <- NULL
for (i in 1:length(p_probeinfo)) {
	p_chrs[i] <- p_probeinfo[[i]][,3]
} 
}
chrs <- NULL
i <- 1
for (id in row.names(p_eb)) {
	chrs[i] <- hugene10sttranscriptclusterCHR[[as.character(id)]]
	i <- i+1
}

# write.table(cbind(n_top_hits,n_genes,n_chrs),file="cd14_n_results.txt",sep="\t")
# write.table(cbind(p_top_hits,p_genes,p_chrs),file="cd14_p_results.txt",sep="\t")

chrs[is.na(chrs)] <- 0

if (FALSE) {
pdf('volcanoa.pdf')
volcanoplot(n_eb)
dev.off()

pdf('volcanob.pdf')
volcanoplot(p_eb)
dev.off()

pdf('hist_n.pdf')
hist(as.numeric(n_top_hits$adj.P.Val))
dev.off()

pdf('hist_p.pdf')
hist(as.numeric(p_top_hits$adj.P.Val))
dev.off()

pdf('diffs.pdf')
hist(diffs)
dev.off()
}
g = ggplot(data=data.frame(p_eb),aes(x=coefficients,y=lods,color=as.factor(chrs),label=chrs)) 

pdf('volcano.pdf')
g + geom_text() + xlab("fold change") + ylab("log odds")
dev.off() 
