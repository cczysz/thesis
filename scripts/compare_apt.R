home = '/home/t.cczysz/'

kill_summary = read.table(file=paste(home,'apt_out/rma.summary.txt',sep=''),header=T,row.names=1)

nokill_summary = read.table(file=paste(home,'apt_nk_out/rma.summary.txt',sep=''),header=T,row.names=1)

# Wilcox comparisons

wilcox_pvals <- numeric()

for (i in 1:dim(kill_summary)[1]) {
	pval <- wilcox.test(as.numeric(kill_summary[i,]),as.numeric(nokill_summary[i,]))$p.value
	wilcox_pvals <- rbind(wilcox_pvals,pval)	
}

pdf(file='wilcox.pdf')
hist(wilcox_pvals)
dev.off()

# Boxplots for QC
pdf(file='kill_bp.pdf')
boxplot(as.matrix(kill_summary))
dev.off()

pdf(file='nokill_bp.pdf')
boxplot(as.matrix(nokill_summary))
dev.off()
