# Initially, testing probe removal for Caucasian individuals in the CD14 celltype
## Will be expanded after testing

# This file will source other scripts which each perform probe removal, normalization, PEER analysis, and differential expression

# Sys.setenv(R_THREADS=4)
library(oligo)
library(limma)
library(ggplot2)

files.dir = "/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD14/"
# out.dir = paste(files.dir,'out',sep='')
out.dir = '/scratch/t.cczysz/'
setwd(files.dir)

# Normalization files
phenotype.file = "CD14.Samples.POSTQC.ImmVarFinal.txt"
peer.factors.f = "/scratch/t.cczysz/peer_factors.Robj"
# raw_exp_rfile = "cd4_cau_raw.RData"
# norm_exp_rfile = "cd4_cau_expr.RData"

probe.info = "/home/t.cczysz/HuGeneProbeInfo.csv"

samples <- read.csv(phenotype.file, header=T)

data.cau <- samples[samples$Race == 'Caucasian', ]
# data.subset <- data.cau[sample(1:nrow(data.cau), 20), ]
data.subset <- data.cau

data.ids.subset <- data.subset[,1]
data.files.subset <- data.subset[,2]
data.sex.subset <- data.subset$Sex

raw.data.subset <- read.celfiles(as.character(data.files.subset))

#row.names(phenoData(rawData)) <- data_subset_ids
# fit1 <- fitProbeLevelModel(rawData,normalize=F,background=F)
# RMA defaults to background subtraction, quantile normalization, and summarization via median polish
data.rma.subset <- rma(raw.data.subset)

# Remove Probes
	# Write code to remove given probes from raw data, perform normalization
# PEER files
# Required for PEER:
	# Raw Expression as `expression`
	# Vector of sex covariates as `sex`
# Outputs:
	# Factors: #(Genes)x#(Factors) matrix of PEER factors for use in linear model
sex <- as.numeric(data.sex.subset == 'Male')
expression <- exprs(data.rma.subset)

if (file.exists(file=peer.factors.f)) load(file=peer.factors.f) else {
  source('/home/t.cri.cczysz/thesis/scripts/peer.R', echo=T)
  save(peer.factors,file=peer.factors.f) 
}

### Calculating Residuals
# Call make_residuals function in apply over rows of probeset expression
dim(expression)
dim(peer.factors)
head(peer.factors)
MakeResiduals <- function(input.row) {
        #residuals(lm(input.row ~ 0 + peer.factors[, 2] + peer.factors[, 3] + peer.factors[, 4] + peer.factors[, 5] + peer.factors[, 6] + peer.factors[, 7] + peer.factors[, 8] + peer.factors[, 9] + peer.factors[, 10] + peer.factors[, 11] + peer.factors[, 12] + peer.factors[, 13] + peer.factors[, 14] + peer.factors[, 15] + peer.factors[, 16] + peer.factors[, 17] + peer.factors[, 18] + peer.factors[, 19] + peer.factors[, 20] + peer.factors[, 21]))
	residuals(lm(input.row ~ 0 + peer.factors[, -1]))
}

expr.residuals <- apply(expression, 1, MakeResiduals)
expr.residuals <- t(expr.residuals)
print(expr.residuals[1:5,1:10])
dim(expr.residuals)
# DE files

samples <- as.factor(data.sex.subset)
design <- model.matrix(~0 + samples)
colnames(design) <- c("Female","Male")
fit <- lmFit(expr.residuals, design)
contrast.matrix <- makeContrasts(mf = Male - Female, levels=design)
contrast.fit <- contrasts.fit(fit, contrast.matrix)
eb.fit <- eBayes(contrast.fit, robust=TRUE)
print(topTable(eb.fit, number=100))

g = ggplot(data=data.frame(eb.fit),aes(x=coefficients,y=lods)) 

pdf('/home/t.cri.cczysz/volcano.pdf')
g + geom_point() + xlab("fold change") + ylab("log odds")
dev.off() 
