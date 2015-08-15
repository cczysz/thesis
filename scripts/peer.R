library(peer)

RunPeer <- function(expression, k=20, covs) {
	model = PEER()
	PEER_setCovariates(model,as.matrix(sex))
	
	# Expression must be in NxG. N number of samples, G number of genes
	PEER_setPhenoMean(model,t(expression))
	PEER_setNk(model,k)
	PEER_update(model)
	peer.factors = PEER_getX(model)
	return(peer.factors)
}
