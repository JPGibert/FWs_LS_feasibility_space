setwd("~/Desktop/JP/Papers_in_review_submitted/JP_FW_constrained_space/Analyses/Food_webs")
library(igraph)


FWs <- dir()
modularity <- rep(0,length(FWs))
for(i in 1:length(FWs)){
	i=5
	mat <- as.matrix(read.table(FWs[i]))
	g <- graph_from_adjacency_matrix(mat, mode="directed")
	wtc <- cluster_walktrap(g)
	modularity[i] <- modularity(wtc)
}
modularity


