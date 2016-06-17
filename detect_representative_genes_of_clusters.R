### Latest Modified: 2016_06_17 by WDF

### 0617 modified: add new parameters: all candidate genes

### Description: Given number of clusters, cluster assignment and gene exp. matrix, 

###	  detect representative genes of every single cluster by:

###   For a specific cluster k ,compare it with all others seperately 

###   1. Detection of significantly differential expressed genes:
###    (1). t.test (fdr<1e-5) and
###    (2). wilcoxon.test (fdr<1e-5)
###   2. Higher expression than other clusters
###
###   The intersection of all comparision results as final signals
options( stringsAsFactors = FALSE )
###	Parameters:
args <- commandArgs(TRUE)
###	  1. EXPRESSION_MATRIX_S1: i.e. processed matrix 
file <- args[1]
print(file)
data <- read.table( file, header = T, sep = '\t' )     #### read in files
data_matrix <- as.matrix( data[,c(-1,-2)] )            #### save only exp data in matrix form
###   2. LABEL_OF_DATA (the same as step1,such as hccdb1,dnameth.iqr10 etc)
output <- args[2]
###   3. NUMBER_OF_CLUSTERS
k <- as.numeric( args[3] )   #### factorization rank of NMF
###   4. Cluster assignment (contains two named columns: Samples___Cluster)
cluster.file <- args[4]
data_cluster <- read.table( cluster.file,header = T,sep = "\t",row.names = 1 )

###	  5. VERSION (the same as step1, such as a date v0616)
version <- args[5]

###   6. Gene file: we hope it contains a colname to specify whether gene entrez id (col: Entrez) or symbol (col: Symbol) is given
gene.file <- args[6]
check_gene <- FALSE
if( !is.na(gene.file) ){
  gene_name <- read.table( gene.file,header = T,sep = "\t" )
  if( "Entrez" %in% colnames(gene_name) ){
    check_gene = TRUE
    check_name = "Entrez"
  }
  if( "Symbol" %in% colnames(gene_name) ){
    check_gene = TRUE
    check_name = "Symbol"
  } 
} else{
  messgae( "NO appropriate candidate genes file is provided." )
}

###
###	examples: 
###
###	Output:
####  Create output folder
folder <- paste0( "~/HCCDB_Subtypes/results/",version,"/" )
system( paste0('mkdir ',folder) )

subfolder <- paste0( "~/HCCDB_Subtypes/results/",version,"/",output,"/" )
system( paste0('mkdir ',subfolder) )

subsubfolder <- paste0( "~/HCCDB_Subtypes/results/",version,"/",output,'/signal_gene_boxplot/' )
system( paste0('mkdir ',subsubfolder) )

###   1. Signals(GENE_ID) of each cluster in folder './results/LABEL_OF_DATA'
###   2. LABEL_OF_DATA.cluster i.e. cluster assignment of each sample
###   3. mkdir subfolder './signal_gene_boxplot/', and in this subfolder:
###			plot boxplot of selected genes in each dataset and labeled 
###			(1). red,cornflowerblue,antiquewhite4,darkmagenta,darkgoldenrod2,forestgreen
###


library(ggplot2)
colors <- c( "red","cornflowerblue","antiquewhite4","darkmagenta","darkgoldenrod2","forestgreen" )

########################### delete duplicate rows ####################################
entrez <- data$Entrez_ID
symbol <- data$Symbol

entrez.count <- table(entrez)
dup <- as.numeric( names( which( entrez.count > 1 ) ) )
delete.index <- c()
for( i in 1:length(dup) ){
  temp <- dup[i]
  temp.index <- which( entrez == temp )
  delete.index <- c(delete.index,-temp.index[-1])
}

if( length(delete.index)>0 ){
  data_matrix <- data_matrix[delete.index,]
  entrez <- entrez[delete.index]
  symbol <- data$Symbol[delete.index]
}

rownames(data_matrix) <- entrez

index <- which( is.na(data_matrix) )
data_matrix[index] <- 0

#####
## rows: genes
## columns: patients
#####

samples <- rownames( data_cluster )
samples.cluster <- data_cluster$Cluster
data_matrix_for_cluster_analysis <- data_matrix[ samples ]

if( check_gene ){
  if( check_name == "Entrez"){
    gene_index <- which( entrez %in% gene_name$Entrez )
    data_matrix_for_cluster_analysis <- data_matrix_for_cluster_analysis[gene_index,]
  }
  if( check_name == "Symbol"){
    gene_index <- which( symbol %in% gene_name$Symbol )
    data_matrix_for_cluster_analysis <- data_matrix_for_cluster_analysis[gene_index,]
  }
}

###################################Obtain signatures####################################
#### Now suppose there are k clusters###############
#### t.test + wilcoxon.test

gene.num <- dim(data_matrix_for_cluster_analysis)[1]

for( cluster in 1:k ){
  cluster1 <- which( cluster.index == cluster )
  genes_ID <- entrez
  genes_Symbol <- symbol
  
  for( comparison in 1:k ){
    if( comparison != cluster ){
      cluster2 <- which( cluster.index == comparison  )
      ### t.test
      t.pval <- array( 0,dim=gene.num )
      names( t.pval ) <- rownames( data_matrix_for_cluster_analysis )
      avg.exp <- matrix( 0,nrow = gene.num,ncol = 2 )
      rownames(avg.exp) <- row.names( data_matrix_for_cluster_analysis )
      ### colnames(avg.exp) <- c( "MeanOfCluster1","MeanOfCluster2" )
      
      for( i in 1:gene.num ){
        t.stat <- t.test( data_matrix_for_cluster_analysis[i,cluster1],data_matrix_for_cluster_analysis[i,cluster2] )
        t.pval[i] <- t.stat$p.value
        avg.exp[i,1] <- t.stat$estimate[1]
        avg.exp[i,2] <- t.stat$estimate[2]
      }
      
      t.pval.adjust <- p.adjust( t.pval,method = "BH" )
      names(t.pval.adjust) <- rownames( data_matrix_for_cluster_analysis )
      
      ### wilcox.test
      wilcox.pval <- array( 0,dim=gene.num )
      names( wilcox.pval ) <- rownames( data_matrix_for_cluster_analysis )
      
      for( i in 1:gene.num ){
        wil.stat <- wilcox.test( data_matrix_for_cluster_analysis[i,cluster1],data_matrix_for_cluster_analysis[i,cluster2] )
        wilcox.pval[i] <- wil.stat$p.value
      }
      
      wilcox.pval.adjust <- p.adjust( wilcox.pval,method = "BH" )
      names(wilcox.pval.adjust) <- rownames( data_matrix_for_cluster_analysis )
      
      ### output 
      data.diff <- data.frame( 
        Genes = rownames(data_matrix_for_cluster_analysis),
        Symbol = symbol,
        MeanOfCluster1 = avg.exp[,1],
        MeanOfCluster2 = avg.exp[,2],
        t.pval = t.pval,
        t.pval.adjust = t.pval.adjust,
        wilcox.pval = wilcox.pval,
        wilcox.pval.adjust = wilcox.pval.adjust
      )
      index1 <- which( data.diff$t.pval.adjust < 1e-5 & data.diff$wilcox.pval.adjust < 1e-5 & data.diff$MeanOfCluster1 > data.diff$MeanOfCluster2 )
      genes_ID <- intersect( data.diff$Genes[index1],genes_ID )
      genes_symbol <- intersect( data.diff$Symbol[index1],genes_symbol )
    }
  }
  
  output_sig <- data.frame(
    genes_ID = genes_ID,
    genes_Symbol = genes_Symbol
  )
  
  write.table( output_sig,file = paste0(subfolder,output,".",cluster,".sig"),quote = F,row.names = F,col.names = F )
  
  number_of_sigs <- length( genes_ID )
  
  print( number_of_sigs )
  
  if( number_of_sigs > 0 ){
    ### plot boxplot
    for( i in 1:number_of_sigs ){
      temp.symbol <- output_sig$genes_Symbol[i]
      temp.id <- output_sig$genes_ID[i]
      pdf( paste0( subsubfolder,temp.id,'_',strsplit(temp.symbol,split = "/")[[1]][1],'.pdf' ),height=3,width=5 )
      data_plot <- data.frame(
        exp = data_matrix_for_cluster_analysis[temp.id,],
        cluster = samples.cluster
      )
      p <- ggplot( data_plot,aes(x=factor(cluster),y=exp,fill=colors[cluster]) )
      p2 <- p+geom_boxplot()+geom_jitter( width=.2 )+theme_bw()+theme(legend.position="none")+xlab("Cluster")+ylab("Expression")+ggtitle(temp.symbol)
      print(p2)
      dev.off()
    }
  }
}






