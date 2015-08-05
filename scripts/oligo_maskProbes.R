files_dir = "/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD4/"
out_dir = paste(files_dir,'out',sep='')

# Normalization files
phenotype_file = "CD4.Samples.POSTQC.ImmVarFinal.txt"
raw_exp_rfile = "cd4_cau_raw.RData"
norm_exp_rfile = "cd4_cau_expr.RData"

# PEER files
peer_factors_file = "cd4_cau_factors.txt"
residuals_file = "cd4_lm.peer_residuals.txt"

# DE files
probeinfo = "/home/t.cczysz/HuGeneProbeInfo.csv"

Sys.setenv(R_THREADS=4)
library(oligo)

setwd(files_dir)

samples <- read.csv(phenotype_file,header=T)

data_cau <- samples[samples$Race == 'Caucasian',]
data_subset <- data_cau[sample(1:length(data_cau),10),]

data_ids <- data_subset[,1]
data_subset_files <- data_subset[,2]

rawData <- read.celfiles(as.character(data_subset_files))

row.names(phenoData(rawData)) <- data_ids
fit1 <- fitProbeLevelModel(rawData,normalize=F,background=F)

RLE(fit1,type='values')
NUSE(fit1,type='values')
dev.off()
