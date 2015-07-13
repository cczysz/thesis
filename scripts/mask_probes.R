probe_dir = '/group/stranger-lab/forCharles/probes_mapping/intersections/'
r_obj_dir = '/group/stranger-lab/forCharles/probes_mapping/Robjects/'

load(paste(r_obj_dir,'maskBADscoresDF.Robj',sep=''))
row.names(maskBADscoresDF) <- maskBADscoresDF$probes

probes_geneMapping <- read.table(file=paste(probe_dir,'probes.gencode.v22.exon.intersect_stats_ensembl.txt',sep=''),
	header=F,
	as.is=T,
	sep=',',
	skip=1)

probeIds <- probes_geneMapping[,1]

numMultiMap <- probes_geneMapping[,2]
multiMappingIdx <- probes_geneMapping[,2] - 1
multiMappingIdx[multiMappingIdx > 0] <- 1

probesOverlappingSnps <<- read.table(file='/home/t.cczysz/probes_SNP_overlap.txt',header=F,as.is=T)

# ProbeSNP <- function(probeID,probeIds){
	#idx <- grep(probeID,probesOverlappingSnps)
	# if (length(idx)==0) {
		# variant = 0
		# return(variant)
	# } else { variant = 1
		#return(variant) }
# }
ProbeSNP <- function(probeID){
	if (is.element(probeID,probesOverlappingSnps)) {
		variant = 0
		return(variant)
	} else { variant = 1
		return(variant) 
	}
}

#SNP_overlap <- apply(as.matrix(probeIds,ncol=1),1,ProbeSNP)
SNP_overlap <- is.element(probeIds,probesOverlappingSnps$V1)
SNP_overlap <- as.numeric(SNP_overlap)

# badScores = runif(length(probeIds),min=0,max=1)
# load(file='/group/stranger-lab/forCharles/probes_mapping/Robjects/exmask.Robj')

corr_table = data.frame(probeIds,multiMappingIdx,numMultiMap,SNP_overlap)
row.names(corr_table) <- corr_table$probeIds

# corr_table <- merge(corr_table,maskBADscoresDF$quality.score,by=0)
