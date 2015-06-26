# Using R 3.1.0
library('oligo')

# Create list of CEL files
CELlist <- list.celfiles('myCELs',full.names=T)
# First argument can be list?
# full.names means absolute path

# If files in WD:
celList <- list.celfiles(getwd())

# Read CEL files into oligo format
rawData <- read.celfiles(celList)
