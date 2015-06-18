library(peer)

setwd('/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD14')

a <- read.table('cd14_expr.txt',header=TRUE,row.names=1)
head(a)
covs <- read.csv("cd14.cau.cov",header=TRUE)
covs <- as.matrix(covs)

model = PEER()

#call covariates, comment out this also if no covariates
PEER_setCovariates(model,covs[,2])

PEER_setPhenoMean(model,t(as.matrix(a)))

PEER_setNk(model,20)

PEER_getNk(model)

PEER_update(model)

factors = PEER_getX(model)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)

pdf('cd14_cau_peer_model.pdf')
PEER_plotModel(model)
dev.off()

write.table(residuals, file="cd14_cau_residuals.txt", sep="\t")
write.table(factors, file="cd14_cau_factors.txt", sep="\t")
write.table(weights, file="cd14_cau_weights.txt", sep="\t")

pdf('cd14_cau_peer_precision.pdf')
plot(precision)
dev.off()
