library(affy)
library(BiocInstaller)

biocLite("simpleaffy")
library(simpleaffy)
biocLite("affyQCReport")
library(affyQCReport)
# library(limma)

setwd("/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD4")
# phen <- read.csv('CD4.Samples.POSTQC.ImmVarFinal.txt',sep="\t",header=TRUE)

phen <- read.csv("CD4.Samples.Caucasian.txt",header=TRUE) # Import Caucasian data, file list made with grep and awk'ing 2nd column

# file_list <- c() #initalize
# file_list <- rbind(file_list,as.character(files[,1]) # save filenames to character array for use in ReadAffy()

# file_list <- c(as.character(files[,1]),recursive=TRUE)

if (file.exists("r_affy_data.RData")) {

	load("r_affy_data.RData")

} else {
	cau_data <- ReadAffy(filenames=as.character(phen$File.Name),phenoData = data.frame(phen))

	save(list=c("cau_data"),file="r_affy_data.RData")

}


if (file.exists("r_affy_data.RData")) {

	load("r_affy_expr.RData")

} 

# QCReport(cau_data,file="rawQC.pdf")
# QCReport(eset,file="normQC.pdf")

samps <- sample(1:ncol(exprs(eset)),10,replace=FALSE)
h <- hist(exprs(eset)[,samps],freq=FALSE,main="Intensity of 10 random, normalized arrays")
# pdf('normhist.pdf')
# plot(density(h))
# dev.off()
pdf('normbox.pdf')
boxplot(exprs(eset)[,samps])
dev.off()
# QCReport(eset,file="normQC.pdf")
