setwd('/home/t.cczysz')

load('mask.Robj')

mean_ps_scores <- function(probeset_id) {
	mean(exmask[[1]]$quality.score[exmask[[1]]$probeset == probeset_id])
}

mean_scores <- apply(as.matrix(unique(as.numeric(exmask[[1]]$probeset)),ncol=1),1,mean_ps_scores)

pdf('mean_ps_scores')
plot(density(mean_scores))
dev.off()
