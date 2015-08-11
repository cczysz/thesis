library(oligo)

###### MAPPING BASED FILTERING

load("probes_mapping/Robjects/merge_probes_DF.Robj")
par.genes=as.character(read.table("probes_mapping/annotations/par.genes.txt")[,1])
Y.degenerate=as.character(read.table("probes_mapping/annotations/chrY.degenerate.txt")[,1])
Y.degenerate_Xhomolog=as.character(read.table("probes_mapping/annotations/chrY.degenerate_Xhomolog.txt")[,1])
Y.all=as.character(read.table("probes_mapping/annotations/chrY.genes.txt")[,2])

### filter probes overlapping snps MAF>0.1 in either EUR, EAS, AFR

probes_snps=merge_probes_DF[merge_probes_DF$overlapping_snps >0,"probeId"]
length(probes_snps)
#[1] 43657

### Filter multimapping probes

probes_multimapping=merge_probes_DF[grep(":",merge_probes_DF$gene_ensembl),"probeId"]
length(probes_multimapping)
#[1] 111453

### Filtering of probes belonging to probesets with cross-hybridization, except PAR

probes_cross=merge_probes_DF[merge_probes_DF$crosshyb_type == 3 & !merge_probes_DF$gene_ensembl%in%par.genes,"probeId"]
length(probes_cross)
# [1] 85790

filt_probes=unique(c(probes_snps,probes_multimapping,probes_cross))
merge_probes_DF_filt=merge_probes_DF[!merge_probes_DF$probeId%in%filt_probes,]
length(filt_probes)
#[1] 117406

###### EXPRESSION BASED FILTERING
### Background-correction, normalization of probe-level expression values.
### Step already done. Load Robj

load("Robjects/oligo.bgSubstractedNormalized.geneMappingProbes.Robj")

### Eliminate filtered probes

normalized2=normalized2[as.character(rownames(normalized2)[!rownames(normalized2)%in%filt_probes]),]


### Calculate densities gene expression. Identify global minimum, which will be the gene expression threshold.
### Previously, store ImmVarID2 for males and females on vectors of the same name

load("Robjects/phen.Robj")
males=colnames(normalized2)[colnames(normalized2)%in%as.character(phen[phen$Sex%in%"Male","ImmVarID2"])]
females=colnames(normalized2)[colnames(normalized2)%in%as.character(phen[phen$Sex%in%"Female","ImmVarID2"])]

pdf("ProbeExps.pdf")
hist(log2(normalized2[,males]),prob=T)
dM <- density(log2(normalized2[,males]))
min_males=optimize(approxfun(dM$x,dM$y),interval=c(3,4))$minimum
abline(v=min_males)
lines(dM)
legend("topright",legend=paste(c("min f(x) = ",min_males,collapse="")))

hist(log2(normalized2[,females]),prob=T)
dF <- density(log2(normalized2[,females]))
min_females=optimize(approxfun(dF$x,dF$y),interval=c(3,4))$minimum
abline(v=min_females)
lines(dF)
legend("topright",legend=paste(c("min f(x) = ",min_females,collapse="")))
dev.off()

median_genexp_females_probes=apply(normalized2[,females],1,median)
median_genexp_males_probes=apply(normalized2[,males],1,median)

###  Filter Y-linked non-PAR probes for which median expression in females is above threshold and distribution of expressions in male and females are equal (wilxon test pval > 10-5)

filt_Y_probe_exp=vector()
for (gene in Y.all[Y.all[!Y.all%in%par.genes]%in%merge_probes_DF_filt$gene_ensembl]) {
	probes=as.character(merge_probes_DF_filt[grep(gene,merge_probes_DF_filt$gene_ensembl),"probeId"]);
	for (probe in probes) {
                pval=wilcox.test(log2(normalized2[as.character(probe),females]),log2(normalized2[as.character(probe),males]))$p.value
		if (pval > 0.00001 && median_genexp_females_probes[probe] > min_females) {
			filt_Y_probe_exp=c(filt_Y_probe_exp,probe)
		}
        }
}
length(filt_Y_probe_exp)
# [1] 318

merge_probes_DF_filt=merge_probes_DF_filt[!merge_probes_DF_filt$probeId%in%filt_Y_probe_exp,]
normalized2=normalized2[as.character(rownames(normalized2)[!rownames(normalized2)%in%filt_Y_probe_exp]),]

exp_genes <- basicRMA(normalized2, as.character(merge_probes_DF_filt$gene_ensembl), normalize=F, background=F)

save(exp_genes,file="Robjects/exp_genes.Robj")
