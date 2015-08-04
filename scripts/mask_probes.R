probe_dir = '/group/stranger-lab/forCharles/probes_mapping/intersections/'
r_obj_dir = '/group/stranger-lab/forCharles/probes_mapping/Robjects/'
home = '/home/t.cczysz'

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

SNP_overlap <- is.element(probeIds,probesOverlappingSnps$V1)
SNP_overlap <- as.numeric(SNP_overlap)

snp_in_probes <- read.table(file=paste(probe_dir,'probes.dbsnp142.intersect.bed',sep=''),header=F,as.is=T)

maf_snps <- read.table(file=paste(probe_dir,'dbsnp142.overlappingprobes_MAF0.1AFR._EAS_EUR.txt',sep=''),header=F,as.is=T)

probe_to_snp <- cbind(snp_in_probes$V4,snp_in_probes$V10)

maf_probes <- probe_to_snp[is.element(probe_to_snp[,2],c(maf_snps$V1)),1]

maf_probes <- unique(sort(maf_probes))

maf_probes_logical <- as.numeric(is.element(probeIds,maf_probes))

# badScores = runif(length(probeIds),min=0,max=1)
# load(file='/group/stranger-lab/forCharles/probes_mapping/Robjects/exmask.Robj')

corr_table = data.frame(probeIds,multiMappingIdx,numMultiMap,SNP_overlap,maf=maf_probes_logical)
row.names(corr_table) <- corr_table$probeIds

probes_to_remove <- sort(row.names(subset(corr_table,multiMappingIdx == 1 | maf == 1)))

probe_list <- read.table(file='/home/t.cczysz/probe.info',header=T,as.is=T,fill=T)

kill_list <- probe_list[is.element(probe_list[,1],probes_to_remove),]
