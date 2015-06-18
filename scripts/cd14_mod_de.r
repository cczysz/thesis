library(affy)
library(ggplot2)
require(ggplot2)
library(limma)
#library(BiocInstaller)
#biocLite("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)

wd = "/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD14"
setwd(wd)

# Rdata file names
raw_exp_rfile = "cd14_cau_raw.RData"
norm_exp_rfile = "cd14_cau_expr.RData"
probeinfofile = "/home/t.cczysz/HuGeneProbeInfo.csv"
peer_factors_file = "cd14_cau_factors.txt"
lm_peer_residuals_file = "cd14_factor_residuals.txt"

phen <- read.csv("CD14.Caucasian.txt",header=TRUE)

if (file.exists(raw_exp_rfile)) {
	load(raw_exp_rfile)
} else {
	cau_data <- ReadAffy(filenames=as.character(phen$FileName),phenoData = data.frame(phen))
	save(list=c("cau_data"),file=raw_exp_rfile)
}

if (file.exists(norm_exp_rfile)) {
	load(norm_exp_rfile)
	# write.table(exprs(eset),file='cd14_expr.txt')
} else {
	# Save normalized, corrected expression data
	eset <- expresso(cau_data,bgcorrect.method="rma",normalize.method="quantiles",pmcorrect.method="pmonly",summary.method="medianpolish")
	save(list=c("eset"),file=norm_exp_rfile)
}

covs <- read.csv('cd14.cau.cov',header=TRUE,row.names=1)

sex <<- unlist(as.matrix(covs[,1],nrow=1))
# exp <- exprs(eset)
make_residuals <- function(in_row) {
	residuals(lm(in_row ~0+ peer_factors[,1] + peer_factors[,2] + peer_factors[,3] + peer_factors[,4] + peer_factors[,5] + peer_factors[,6] + peer_factors[,7] + peer_factors[,8] + peer_factors[,9] + peer_factors[,10] + peer_factors[,11] + peer_factors[,12] + peer_factors[,13] + peer_factors[,14] + peer_factors[,15] + peer_factors[,16] + peer_factors[,17] + peer_factors[,18] + peer_factors[,19] + peer_factors[,20]))
}

if (file.exists(file=lm_peer_residuals_file)) {
	res <- as.matrix(read.table(file=lm_peer_residuals_file,header=TRUE,row.names=1))
	} else {
	peer_factors <- as.matrix(read.table(file=peer_factors_file,header=TRUE,row.names=1,as.is=TRUE))

	res <- apply(exp,1,make_residuals)
	res <- t(res)
	write.table(res,file=lm_peer_residuals_file)
}

make_ebayes_fit <- function(exp,samples) {

	samples <- as.factor(samples)
	design <- model.matrix(~0 + samples)
	colnames(design) <- c("Male","Female")
	fit <- lmFit(exp,design)
	contrast.matrix <- makeContrasts(Female-Male,levels=design)
	contrast_fit <- contrasts.fit(fit,contrast.matrix)
	eb_Fit <- eBayes(contrast_fit,robust=TRUE)
	return(eb_Fit)
}

samples <- cau_data$Sex
n_eb <- make_ebayes_fit(exprs(eset),samples)
p_eb <- make_ebayes_fit(res,samples)

mean_sex_exp <- function(id,res) {
	# data.frame(id=id,male=mean(res[id,samples=="Male"],female=mean(res[id,samples=="Female"])))
	male <- res[id,!!sex]
	female <- res[id,!sex]
	result <- cbind(id,mean(male),mean(female))
}
sex_expr <- apply(as.matrix(row.names(p_eb),ncol=1),1,mean_sex_exp,res=res)
# head(sex_expr)
sex_expr <- t(sex_expr)
sex_expr <- data.frame(id=sex_expr[,1],male=sex_expr[,2],female=sex_expr[,3])
chrs <- NULL
i <- 1
for (id in row.names(data.frame(p_eb))) {
	chrs[i] <- hugene10sttranscriptclusterCHR[[as.character(id)]]
	i <- i+1
}

chrs[is.na(chrs)] <- 0

s <- ggplot(sex_expr,aes(x=male,y=female,color=as.factor(chrs),label=chrs))
pdf('malevsfemale.pdf')
s + geom_label() + theme(axis.ticks=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())
dev.off()

n_top_hits <<- topTable(n_eb,number=1000,coef=1)
p_top_hits <<- topTable(p_eb,number=1000,coef=1)

comp <- function(id) {
	a <- n_top_hits[id,]
	b <- p_top_hits[id,]
	
	a_p <- a$adj.P.Value 
	b_p <- b$adj.P.Value

	if (is.null(a_p) || is.null(b_p)){
		return(1)
	} else {	
		return(as.numeric(a_p) - as.numeric(b_p))
	}
}

diffs <- apply(as.matrix(row.names(n_top_hits)),1,comp)
# Returns norm - peer p_values

get_gene <- function(id) {
	gene <- hugene10sttranscriptclusterGENENAME[[as.character(id)]]
	chr <- hugene10sttranscriptclusterCHR[[as.character(id)]]
	return(as.matrix(c(id,gene,chr),nrow=1))
}

n_ids <- as.matrix(row.names(n_top_hits),ncol=1)
p_ids <- as.matrix(row.names(p_top_hits),ncol=1)
write.table(p_ids,file="cd14_peer_topprobes.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

probe_info <- read.csv(probeinfofile,header=TRUE)
probe_info <- data.frame(probe_info,row.names=2)
# n_probeinfo <- apply(n_ids,1,get_gene)
# p_probeinfo <- apply(p_ids,1,get_gene)
# Gx1 list. Index each row with [[i]]. Returns 1x3 character array with id, gene_name, and chr as [,1],[,2],and [,3] indicies

get_info <- function(id) {unlist(probe_info[id,])}
# top_probe_info <- NULL
top_probe_info <- lapply(p_ids,get_info)
# top_probe_info <- rbind(top_probe_info,lapply(p_ids,get_info))

write.csv(top_probe_info,file="top_probes.csv")
write.table(top_probe_info,file="top_probes.txt")

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


chrs <- NULL
i <- 1
for (id in row.names(data.frame(p_eb))) {
	chrs[i] <- hugene10sttranscriptclusterCHR[[as.character(id)]]
	i <- i+1
}

go_ids <- NULL
i <- 1
for (id in p_ids) {
	res <- hugene10sttranscriptclusterGO[[as.character(id)]]
	if (!is.na(res)){
		go_ids[i] <- res$GOID
		i <- i+1
	}
}

chrs <- NULL
i <- 1
for (id in row.names(data.frame(p_eb))) {
	chrs[i] <- hugene10sttranscriptclusterCHR[[as.character(id)]]
	i <- i+1
}
}
# chrs[is.na(chrs)] <- 0

# g = ggplot(data=data.frame(p_eb),aes(x=coefficients,y=lods,color=as.factor(chrs),label=chrs)) 

# pdf('volcano.pdf')
# g + geom_text() + xlab("fold change") + ylab("log odds") + ylim(0,25)
#dev.off() 
