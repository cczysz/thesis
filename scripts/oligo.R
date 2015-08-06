# Initially, testing probe removal for Caucasian individuals in the CD14 celltype
## Will be expanded after testing

# This file will source other scripts which each perform probe removal, normalization, PEER analysis, and differential expression

Sys.setenv(R_THREADS=4)
library(oligo)
library(limma)

files.dir = "/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD14/"
# out.dir = paste(files.dir,'out',sep='')
out.dir = '/scratch/t.cczysz/'
setwd(files.dir)

# Normalization files
phenotype.file = "CD14.Samples.POSTQC.ImmVarFinal.txt"
# raw_exp_rfile = "cd4_cau_raw.RData"
# norm_exp_rfile = "cd4_cau_expr.RData"

probe.info = "/home/t.cczysz/HuGeneProbeInfo.csv"

samples <- read.csv(phenotype_file, header=T)

data.cau <- samples[samples$Race == 'Caucasian', ]
data.subset <- data.cau[sample(1:length(data.cau), 10), ]

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
sex <- as.numeric(data_subset_sex == 'Male')
expression <- exprs(data.rma.subset)
source('/home/t.cri.cczysz/thesis/scripts/peer.r', echo=T)

### Calculating Residuals
# Call make_residuals function in apply over rows of probeset expression
MakeResiduals <- function(in.row) {
        residuals(lm(in.row ~ 0 + peer_factors[, 1] + peer_factors[, 2] + peer_factors[, 3] + peer_factors[, 4] + peer_factors[, 5] + peer_factors[, 6] + peer_factors[, 7] + peer_factors[, 8] + peer_factors[, 9] + peer_factors[, 10] + peer_factors[, 11] + peer_factors[, 12] + peer_factors[, 13] + peer_factors[, 14] + peer_factors[, 15] + peer_factors[, 16] + peer_factors[, 17] + peer_factors[, 18] + peer_factors[, 19] + peer_factors[, 20]))
}
residuals <- apply(as.matrix(expression), 1, MakeResiduals)
# DE files

samples <- as.factor(sex)
design <- model.matrix(~0 + samples)
colnames(design) <- c("Female","Male")
fit <- lmFit(exp, design)
contrast.matrix <- makeContrasts(mf = Male - Female, levels=design)
contrast.fit <- contrasts.fit(fit, contrast.matrix)
eb.fit <- eBayes(contrast.fit, robust=TRUE)
