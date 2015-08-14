library(oligo)
setwd('/group/stranger-lab/moliva/ImmVar/')

###### MAPPING BASED FILTERING

load("probes_mapping/Robjects/merge_probes_DF.Robj")

par.genes=as.character(read.table("probes_mapping/annotations/par.genes.txt")[,1])
Y.degenerate=as.character(read.table("probes_mapping/annotations/chrY.degenerate.txt")[,1])
Y.degenerate_Xhomolog=as.character(read.table("probes_mapping/annotations/chrY.degenerate_Xhomolog.txt")[,1])
Y.all=as.character(read.table("probes_mapping/annotations/chrY.genes.txt")[,2])

### filter probes overlapping snps MAF>0.1 in either EUR, EAS, AFR

merge_probes_DF=unique(merge_probes_DF[,c("probeId","gene_ensembl","overlapping_snps","crosshyb_type")])

## Different cross-hybridization inexes for some probes. To pick a single value per probe, assign the most stringent (highest multimapping index) one:

dup=as.character(merge_probes_DF$probeId[duplicated(merge_probes_DF$probeId)])
index_rows=seq(nrow(merge_probes_DF))

# Turning into apply

RmXHybrid <- function(d,merge_probes_DF=merge_probes_DF,index_rows=index_rows) {

	row=index_rows[merge_probes_DF$probeId%in%dup[d]]; 
	cross=as.character(merge_probes_DF[row,"crosshyb_type"]); 
	if (sum(as.numeric(cross%in%c("1","3"))) == 2) { merge_probes_DF[row[1],"crosshyb_type"] = "3"} 
	else if (sum(as.numeric(cross%in%c("1","0"))) == 2) { merge_probes_DF[row[1],"crosshyb_type"] = "0";} 
	else if (sum(as.numeric(cross%in%c("3","0"))) == 2) { merge_probes_DF[row[1],"crosshyb_type"] = "0";} 
	else {break}; 
	return(row[2])
}

if (file.exists('/scratch/t.cczysz/toremove.txt')) to.remove <- read.table(file='/scratch/t.cczysz/toremove.txt',header=F) else {
	to.remove <- apply(as.matrix(seq(length(dup))),1,RmXHybrid,merge_probes_DF=merge_probes_DF,index_rows=index_rows)
	write(file='/scratch/t.cczysz/toremove.txt',as.matrix(to.remove))
}

if (F) {
for (d in seq(length(dup))) { 
	row=index_rows[merge_probes_DF$probeId%in%dup[d]]; 
	cross=as.character(merge_probes_DF[row,"crosshyb_type"]); 
	if (sum(as.numeric(cross%in%c("1","3"))) == 2) { merge_probes_DF[row[1],"crosshyb_type"] = "3"} 
	else if (sum(as.numeric(cross%in%c("1","0"))) == 2) { merge_probes_DF[row[1],"crosshyb_type"] = "0";} 
	else if (sum(as.numeric(cross%in%c("3","0"))) == 2) { merge_probes_DF[row[1],"crosshyb_type"] = "0";} 
	else {break}; 
	merge_probes_DF=merge_probes_DF[-row[2],]; 
	length(index_rows) <- length(index_rows)-1;
}
}

to.remove <- as.numeric(unlist(to.remove))
merge_probes_DF <- merge_probes_DF[-to.remove,]
## Counts after pre-filtering

length(unique(merge_probes_DF$gene_ensembl))
length(unique(merge_probes_DF$probeId))

probes_snps=unique(merge_probes_DF[merge_probes_DF$overlapping_snps >0,"probeId"])
length(probes_snps)
# [1] 39846

### Filter multimapping probes

probes_multimapping=unique(merge_probes_DF[grep(":",merge_probes_DF$gene_ensembl),"probeId"])
length(probes_multimapping)
#[1] 73674

### Filtering of probes belonging to probesets with cross-hybridization, except PAR

probes_cross=unique(merge_probes_DF[merge_probes_DF$crosshyb_type != 1  & !merge_probes_DF$gene_ensembl%in%par.genes,"probeId"])
length(probes_cross)
# [1] 47699

filt_probes=unique(c(probes_snps,probes_multimapping,probes_cross))
merge_probes_DF_filt=merge_probes_DF[!merge_probes_DF$probeId%in%filt_probes,]
length(filt_probes)
#[1] 122062

## Counts after 1st filtering

length(unique(merge_probes_DF_filt$gene_ensembl))
length(unique(merge_probes_DF_filt$probeId))

###### EXPRESSION BASED FILTERING
### Background-correction, normalization of probe-level expression values.
### Step already done. Load Robj

#load("Robjects/oligo.bgSubstractedNormalized.geneMappingProbes.Robj")
source('/group/stranger-lab/moliva/ImmVar/scripts/generate_exp_object.R')

#raw_data <- load.cel.files("Caucasian","CD14")
#normalized2 <- backcorrect.normalize.probe.level(load.cel.files("Caucasian","CD14"))
normalized2 <- load.cel.files("Caucasian","CD14")

### Eliminate filtered probes

#normalized2=unique(normalized2)
normalized2=normalized2[as.character(rownames(normalized2)[!rownames(normalized2)%in%filt_probes]),]

normalized2 <- backcorrect.normalize.probe.level(normalized2)
### Calculate densities gene expression. Identify global minimum, which will be the gene expression threshold.
### Previously, store ImmVarID2 for males and females on vectors of the same name

load("Robjects/phen.Robj")
#males=colnames(normalized2)[colnames(normalized2)%in%as.character(phen[phen$Sex%in%"Male","ImmVarID2"])]
males <- phen[phen$Race=="Caucasian", ]$Sex == "Male"
females <- phen[phen$Race=="Caucasian", ]$Sex == "Female"
#females=colnames(normalized2)[colnames(normalized2)%in%as.character(phen[phen$Sex%in%"Female","ImmVarID2"])]

dM <- density(log2(exprs(normalized2)[,males]))
min_males=optimize(approxfun(dM$x,dM$y),interval=c(3,4))$minimum

dF <- density(log2(exprs(normalized2)[,females]))
min_females=optimize(approxfun(dF$x,dF$y),interval=c(3,4))$minimum

if (F) {
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
}

median_genexp_females_probes=apply(log2(exprs(normalized2)[,females]),1,median)
median_genexp_males_probes=apply(log2(exprs(normalized2)[,males]),1,median)

###  Filter Y-linked non-PAR probes for which median expression in females is above threshold and distribution of expressions in male and females are equal (wilxon test pval > 10-5)

filt_Y_probe_exp=vector()
for (gene in unique(as.character(merge_probes_DF_filt[merge_probes_DF_filt$gene_ensembl%in%Y.all[!Y.all%in%par.genes],"gene_ensembl"]))) {
	probes=as.character(merge_probes_DF_filt[grep(gene,merge_probes_DF_filt$gene_ensembl),"probeId"]);
	for (probe in probes) {
                pval=wilcox.test(log2(exprs(normalized2)[probe,females]),log2(exprs(normalized2)[probe,males]))$p.value
		if (pval > 0.00001 && median_genexp_females_probes[probe] > min_females) {
			filt_Y_probe_exp=c(filt_Y_probe_exp,probe)
		}
        }
}
filt_Y_probe_exp=unique(filt_Y_probe_exp)
length(filt_Y_probe_exp)
# [1] 411

### Filter probes below exp threshold

#filt_probe_min_exp=vector()
#for (gene in unique(as.character(merge_probes_DF_filt$gene_ensembl))) {
#        probes=as.character(merge_probes_DF_filt[grep(gene,merge_probes_DF_filt$gene_ensembl),"probeId"]);
#        for (probe in probes) {
#                if (median_genexp_males_probes[probe] < min_males && median_genexp_females_probes[probe] < min_females) {
#                        filt_probe_min_exp=c(filt_probe_min_exp,probe)
#                }
#        }
#}
#length(filt_probe_min_exp)

#filt_probes=c(filt_Y_probe_exp,filt_probe_min_exp)
filt_probes=c(filt_Y_probe_exp)

merge_probes_DF_filt=merge_probes_DF_filt[!merge_probes_DF_filt$probeId%in%filt_probes,]
normalized2=normalized2[as.character(rownames(normalized2)[!rownames(normalized2)%in%filt_probes]),]

## Counts after 2nd filtering

length(unique(merge_probes_DF_filt$gene_ensembl))
length(unique(merge_probes_DF_filt$probeId))

exp_genes <- basicRMA(exprs(normalized2), as.character(merge_probes_DF_filt$gene_ensembl), normalize=F, background=F)
probes_per_gene=as.data.frame(unlist(table(as.character(merge_probes_DF_filt$gene_ensembl))))

if (F) {
pdf("Probes_per_gene.pdf")
hist(probes_per_gene$Freq,breaks=seq(max(probes_per_gene$Freq)))
abline(v=8)
dev.off()
}

### Filter for genes with < 8 probes

genes=as.character(probes_per_gene$Var1[probes_per_gene$Freq > 7])
length(genes)
length(merge_probes_DF_filt$probeId[merge_probes_DF_filt$gene_ensembl%in%genes])
exp_genes=exp_genes[genes,]

save(exp_genes,file="/scratch/t.cczysz/exp_genes.Robj")
