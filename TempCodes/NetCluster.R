#### Latest Modified: 2015-05-18 BY WDF

##### Description: Construct signatures network and make communities

##### Parameters: node_file, edge_file, version

##### Output: node_cluster

args <- commandArgs(TRUE)
node_file_name <- args[1]
edge_file_name <- args[2]
prefix <- args[3]

options(stringsAsFactors = FALSE)
library(igraph)

edge.frame <- read.table( edge_file_name,header = T,sep = '\t')
mygraph <- graph_from_data_frame( edge.frame,directed = F )

node.frame <- read.table(node_file_name,header = T,sep = '\t', row.names = 1)
node.attr <- node.frame[ names(V(mygraph)), ]

vertex_attr( mygraph,'color' ) <- node.attr$Dataset
vertex_attr( mygraph,'size' ) <- node.attr$Size

# res <- cluster_edge_betweenness( mygraph,weights = E(mygraph)$weight_hyper_geometric_log_p, directed = F )

res <- cluster_walktrap( mygraph, weights = E(mygraph)$weight_hyper_geometric_log_p )

# pdf( paste0('../results/NETWORK/network_',version,'.pdf'), width=20,height = 20)
# plot(res,mygraph)
# dev.off()

#cluster <- cut_at(res,no=3)
cluster <- membership( res )

print(cluster)

cluster.assign <- as.numeric(cluster)
nodes <- names(cluster)


data.out <- data.frame(
  nodes = names(V(mygraph)),#names(cluster),
  cluster = cluster.assign
)

write.table( data.out,file = paste0(prefix,"node_cluster.txt"),quote = F,sep = '\t',row.names = F )

pdf( paste0(prefix,'dendrogram.pdf'),width = 20,height = 10 )
plot_dendrogram(res)
dev.off()

