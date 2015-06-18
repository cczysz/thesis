library(peer)

setwd('/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD4')

a <- read.table('cd4_expr.txt',header=TRUE,sep="\t",row.names=1,as.is=TRUE)
# a <- read.table('mod_resid.txt',header=TRUE,row.names=1)
dim(a)
#a <- a[,-1]
# head(a)
covs <- read.csv("cd4.cov",header=TRUE,row.names=1)
covs <- as.matrix(covs)

model = PEER()

PEER_setCovariates(model,covs[,1])
# Don't need covariates

PEER_setPhenoMean(model,t(as.matrix(a)))

PEER_setNk(model,20)

# PEER_getNk(model)

PEER_update(model)

factors = PEER_getX(model)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)

pdf('cd4_mod_peer_model.pdf')
PEER_plotModel(model)
dev.off()

write.table(residuals, file="cd4_mod_residuals.txt", sep="\t")
write.table(factors, file="cd4_mod_factors.txt", sep="\t")
write.table(weights, file="cd4_mod_weights.txt", sep="\t")

pdf('cd4_mod_peer_precision.pdf')
plot(precision)
dev.off()
