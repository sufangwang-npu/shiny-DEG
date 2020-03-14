library(cluster)
# this part is calculating distance based on different choices of distance metric
cmcdistance <- function (dist.choice,zscore.mat) {
	if (dist.choice == "Pearson correlation distance") {
		cluster.cor <- cor(zscore.mat, use="pairwise.complete.obs", method="pearson")
		cluster.dist <- as.dist(1 - cluster.cor)
	} else if (dist.choice == "Euclidean") {
		cluster.dist <- daisy(t(zscore.mat), metric="euclidean", stand=F)
	} else {
		cluster.dist <- daisy(t(zscore.mat), metric="manhattan", stand=F)
	}
	return (cluster.dist)	
}
		
		