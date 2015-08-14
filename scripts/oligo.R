# Initially, testing probe removal for Caucasian individuals in the CD14 celltype
## Will be expanded after testing

# This file will source other scripts which each perform probe removal, normalization, PEER analysis, and differential expression

# Sys.setenv(R_THREADS=4)
library(oligo)
library(limma)
library(ggplot2)

LoadData <- function(population,cell.type) {
	files.dir = "/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD14/"
	setwd(files.dir)

	out.dir = '/scratch/t.cczysz/'

	phenotype.file = "CD14.Samples.POSTQC.ImmVarFinal.txt"
	probe.info = "/home/t.cri.cczysz/HuGeneProbeInfo.csv"

	# Save files
	peer.factors.f = "/scratch/t.cczysz/peer_factors.Robj"
	residual.exp.f = "/scratch/t.cczysz/residuals.Robj"
	# raw_exp_rfile = "cd4_cau_raw.RData"
	# norm_exp_rfile = "cd4_cau_expr.RData"

	samples <- read.csv(phenotype.file, header=T)

	data.cau <- samples[samples$Race == 'Caucasian', ]
	#data.subset <- data.cau[sample(1:nrow(data.cau), 10), ]
	data.subset <- data.cau

	data.ids.subset <- data.subset[,1]
	data.files.subset <- data.subset[,2]
	data.sex.subset <- data.subset$Sex

	raw.data.subset <- read.celfiles(as.character(data.files.subset))

	# Remove Probes
		# Write code to remove given probes from raw data, perform normalization
	if (F) {
	bg <- oligo::backgroundCorrect(raw.data.subset)
	normalized2 <- normalize(bg)
	}
}

source('/home/t.cri.cczysz/thesis/scripts/generate_exp_object.R')

raw.data <- load.cel.files("Caucasian","CD14")

source('/home/t.cri.cczysz/thesis/scripts/removeProbes.R')

# RMA defaults to background subtraction, quantile normalization, and summarization via median polish
# data.rma.subset <- rma(raw.data.subset,target=NULL)
# PEER files
# Required for PEER:
	# Raw Expression as `expression`
	# Vector of sex covariates as `sex`
# Outputs:
	# Factors: #(Genes)x#(Factors) matrix of PEER factors for use in linear model
sex <- as.numeric(data.sex.subset == 'Male')
expression <- exprs(data.rma.subset)

if (file.exists(file=peer.factors.f)) load(file=peer.factors.f) else {
  source('/home/t.cri.cczysz/thesis/scripts/peer.R')
  peer.factors <- RunPeer(expression,k=20,sex)
  save(peer.factors,file=peer.factors.f) 
}

### Calculating Residuals
# Call make_residuals function in apply over rows of probeset expression
if (file.exists(file=residual.exp.f)) load(file=residual.exp.f) else {
	MakeResiduals <- function(input.row) {
		residuals(lm(input.row ~ 0 + peer.factors[, -1]))
	}
	expr.residuals <- apply(expression, 1, MakeResiduals)
	expr.residuals <- t(expr.residuals)
}

samples <- as.factor(data.sex.subset)
design <- model.matrix(~0 + samples)
colnames(design) <- c("Female","Male")
fit <- lmFit(expr.residuals, design)
contrast.matrix <- makeContrasts(mf = Male - Female, levels=design)
contrast.fit <- contrasts.fit(fit, contrast.matrix)
eb.fit <- eBayes(contrast.fit, robust=TRUE)

top.probesets <- row.names(topTable(eb.fit, number=100))

GetProbesetInfo <- function(probeset){
	cmd <- paste("grep", as.character(probeset), probe.info)
	system(command=cmd, intern=T)
}

# probeset.info <- apply(as.matrix(top.probesets),1,GetProbesetInfo)
probeset.info <- read.csv(file=probe.info,header=T,row.names=1,as.is=T,strip.white=T)
probeset.info <- as.matrix(probeset.info)

de.genes <- probeset.info[rownames(topTable(eb.fit, number=100)), ]
head(de.genes)
#write(probeset.info[rownames(topTable(eb.fit, number=100)), ], file="/scratch/t.cczysz/de_probe.info")
#write(probeset.info,file="de_probe.info")

g = ggplot(data=data.frame(eb.fit),aes(x=coefficients,y=lods)) 

pdf('/home/t.cri.cczysz/volcano.pdf')
g + geom_point() + xlab("fold change") + ylab("log odds")
dev.off() 
