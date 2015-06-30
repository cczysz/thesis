# Using R 3.1.0
#library('biocInstaller')
biocLite('oligoData')
library('oligo')

data(affyExpressionFS)

sns <- sampleNames(affyExpressionFS)
sns <- gsub('1521','1251',sns)
sns <- gsub('r\\.CEL$','\\.CEL',sns)

wafer <- substr(sns,1,4)
experiment <- substr(sns,5,5)
tmp <- substr(sns,6,7)
complex <- rep('+',length(tmp))
complex[tmp == '00'] <- '-'

info <- data.frame(wafer=wafer,experiment=experiment,complex=complex)
rownames(info) <- sns
metadata <- data.frame(labelDescription=c('wafer','experiment','complex'))
sampleNames(affyExpressionFS) <- sns
pd <- new('AnnotatedDataFrame',data=info,varMetadata=metadata)
phenoData(affyExpressionFS) <- pd

# rm(tmp,wafer,experiment,complex,pd,metadata)
# Create list of CEL files
# CELlist <- list.celfiles('myCELs',full.names=T)
# First argument can be list?
# full.names means absolute path

# If files in WD:
# celList <- list.celfiles(getwd())

# Read CEL files into oligo format
# rawData <- read.celfiles(celList)
