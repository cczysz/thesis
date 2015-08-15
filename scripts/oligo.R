# Initially, testing probe removal for Caucasian individuals in the CD14 celltype
## Will be expanded after testing

# This file will source other scripts which each perform probe removal, normalization, PEER analysis, and differential expression

# Sys.setenv(R_THREADS=4)
library(oligo)
library(limma)
library(ggplot2)

load('/group/stranger-lab/moliva/ImmVar/probes_mapping/Robjects/merge_probes_DF.Robj')
#LoadData <- function(population,cell.type) {
files.dir = "/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD14/"
setwd(files.dir)

out.dir = '/scratch/t.cczysz/'

phenotype.file = "CD14.Samples.POSTQC.ImmVarFinal.txt"
probe.info = "/home/t.cri.cczysz/HuGeneProbeInfo.csv"

# Save files
peer.factors.f = "/home/t.cri.cczysz/peer_factors.Robj"
residual.exp.f = "/home/t.cri.cczysz/residuals.Robj"
# raw_exp_rfile = "cd4_cau_raw.RData"
# norm_exp_rfile = "cd4_cau_expr.RData"

samples <- read.csv(phenotype.file, header=T)

data.cau <- samples[samples$Race == 'Caucasian', ]
#data.subset <- data.cau[sample(1:nrow(data.cau), 10), ]
data.subset <- data.cau

data.ids.subset <- data.subset[,1]
data.files.subset <- data.subset[,2]
data.sex.subset <- data.subset$Sex

#raw.data.subset <- read.celfiles(as.character(data.files.subset))

# Remove Probes
	# Write code to remove given probes from raw data, perform normalization
if (F) {
bg <- oligo::backgroundCorrect(raw.data.subset)
normalized2 <- normalize(bg)
}
#}

load('/group/stranger-lab/moliva/ImmVar/Robjects/phen.Robj')
source('/home/t.cri.cczysz/thesis/scripts/generate_exp_object.R')

#raw.data <- load.cel.files("Caucasian","CD14")

# source('/home/t.cri.cczysz/thesis/scripts/removeProbes.R')

# RMA defaults to background subtraction, quantile normalization, and summarization via median polish
# data.rma.subset <- rma(raw.data.subset,target=NULL)
# PEER files
# Required for PEER:
	# Raw Expression as `expression`
	# Vector of sex covariates as `sex`
# Outputs:
	# Factors: #(Genes)x#(Factors) matrix of PEER factors for use in linear model

load('/scratch/t.cczysz/exp_genes.Robj')
#load('/group/stranger-lab/moliva/ImmVar/Robjects/exp_genes.Robj')
expression <- exp_genes

if (file.exists(file=peer.factors.f)) load(file=peer.factors.f) else {
  source('/home/t.cri.cczysz/thesis/scripts/peer.R')
  sex <- as.numeric(phen[phen$Race == 'Caucasian', ]$Sex == 'Male')
  peer.factors <- RunPeer(expression,k=20,sex)
  save(peer.factors, file=peer.factors.f) 
}
save(peer.factors, file='/scratch/t.cczysz/peer_factors.Robj') 

### Calculating Residuals
if (file.exists(file=residual.exp.f)) load(file=residual.exp.f) else {
	MakeResiduals <- function(input.row,peer.factors) {
		residuals(lm(input.row ~ 0 + peer.factors[, -1]))
	}
	expr.residuals <- apply(expression, 1, MakeResiduals, peer.factors=peer.factors)
	expr.residuals <- t(expr.residuals)
}

PerformDEAnalysis <- function(expr,samples) {
	design <- model.matrix(~0 + samples)
	colnames(design) <- c("Female","Male")

	fit <- lmFit(expr, design)

	contrast.matrix <- makeContrasts(mf = Male - Female, levels=design)
	contrast.fit <- contrasts.fit(fit, contrast.matrix)
	eb.fit <- eBayes(contrast.fit, robust=TRUE)
}

eb.fit <- PerformDEAnalysis(expr.residuals,as.factor(sex))

top.de.genes <- row.names(topTable(eb.fit, number=100))

DE_probesets <- unique(merge_probes_DF[merge_probes_DF$gene_ensembl%in%top.de.genes, ]$fsetid)

g = ggplot(data=data.frame(eb.fit),aes(x=coefficients,y=lods)) 
pdf('/home/t.cri.cczysz/volcano.pdf')
g + geom_point() + xlab("fold change") + ylab("log odds")
dev.off() 

pdf(file='/home/t.cri.cczysz/de_probesets.pdf')
for (set in top.probesets) {
	plot(density(expr.residuals[set, !!sex]),col='blue',xlim=c(-1,1),ylim=c(0,1))
	lines(density(expr.residuals[set, !sex]),col='red')
}
dev.off()
