# Initially, testing probe removal for Caucasian individuals in the CD14 celltype
## Will be expanded after testing

# This file will source other scripts which each perform probe removal, normalization, PEER analysis, and differential expression

Sys.setenv(R_THREADS=4)
library(oligo)

files_dir = "/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD14/"
# out_dir = paste(files_dir,'out',sep='')
out_dir = '/scratch/t.cczysz/'
setwd(files_dir)

# Normalization files
phenotype_file = "CD14.Samples.POSTQC.ImmVarFinal.txt"
# raw_exp_rfile = "cd4_cau_raw.RData"
# norm_exp_rfile = "cd4_cau_expr.RData"

probeinfo = "/home/t.cczysz/HuGeneProbeInfo.csv"

samples <- read.csv(phenotype_file,header=T)

data_cau <- samples[samples$Race == 'Caucasian',]
data_subset <- data_cau[sample(1:length(data_cau),10),]

data_subset_ids <- data_subset[,1]
data_subset_files <- data_subset[,2]
data_subset_sex <- data_subset$Sex

rawData <- read.celfiles(as.character(data_subset_files))

#row.names(phenoData(rawData)) <- data_subset_ids
# fit1 <- fitProbeLevelModel(rawData,normalize=F,background=F)
# RMA defaults to background subtraction, quantile normalization, and summarization via median polish
data_rma <- rma(rawData)

# Remove Probes
	# Write code to remove given probes from raw data, perform normalization
# PEER files
# Required for PEER:
	# Raw Expression as `expression`
	# Vector of sex covariates as `sex`
# Outputs:
	# Factors: #(Genes)x#(Factors) matrix of PEER factors for use in linear model
sex <- as.numeric(data_subset_sex == 'Male')
expression <- exprs(data_rma)
source('/home/t.cri.cczysz/thesis/scripts/peer.r',echo=T)

### Calculating Residuals
# Call make_residuals function in apply over rows of probeset expression
make_residuals <- function(in_row) {
        residuals(lm(in_row ~ 0 + peer_factors[,1] + peer_factors[,2] + peer_factors[,3] + peer_factors[,4] + peer_factors[,5] + peer_factors[,6] + peer_factors[,7] + peer_factors[,8] + peer_factors[,9] + peer_factors[,10] + peer_factors[,11] + peer_factors[,12] + peer_factors[,13] + peer_factors[,14] + peer_factors[,15] + peer_factors[,16] + peer_factors[,17] + peer_factors[,18] + peer_factors[,19] + peer_factors[,20]))
}
residuals <- apply(as.matrix(expression),1,make_residuals)
# DE files
