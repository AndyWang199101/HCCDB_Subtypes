### Latest Modified: 2016_06_17 by WDF

### 0617 modified: add new parameters: all candidate genes

### 0621 modified: add anova analysis

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
data <- read.table( file, header = T, sep = ',' )     #### read in files
data_matrix <- as.matrix( data[,c(-1,-2)] )            #### save only exp data in matrix form

message( "------ Dimension of original data matrix ------" )
print( dim(data_matrix) )

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
} 
if( !check_gene){
  message( "NO appropriate candidate genes file is provided." )
}
  
##### DEFINE: the fdr threshold #####
fdr.threshold <- 1e-5
ANOVA <- FALSE

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

if( length(delete.index)>0 ){
  data_matrix <- data_matrix[delete.index,]
  entrez <- entrez[delete.index]
  symbol <- data$Symbol[delete.index]
}

symbol.count <- table(symbol)
dup <- names( which( symbol.count > 1 ) ) 
delete.index.symbol <- c()
for( i in 1:length(dup) ){
  temp <- dup[i]
  temp.index <- which( symbol == temp )
  delete.index.symbol <- c(delete.index.symbol,-temp.index[-1])
}

if( length(delete.index.symbol)>0 ){
  data_matrix <- data_matrix[delete.index.symbol,]
  entrez <- entrez[delete.index.symbol]
  symbol <- data$Symbol[delete.index.symbol]
}

############################  delete sd = 0 rows ######################################
index <- which( is.na(data_matrix) )
data_matrix[index] <- 0

delete.index.sd <- c()
for( i in 1:dim(data_matrix)[1] ){
  if( sd( data_matrix[i,] ) == 0 ){
    delete.index.sd <- c( delete.index.sd,-i )
  }
}
print( length(delete.index.sd) )
if( length(delete.index.sd)>0 ){
  data_matrix <- data_matrix[delete.index.sd,]
  entrez <- entrez[delete.index.sd]
  symbol <- data$Symbol[delete.index.sd]
}

rownames(data_matrix) <- entrez

#####
## rows: genes
## columns: patients
#####

samples <- rownames( data_cluster )
samples.cluster <- data_cluster$Cluster
data_matrix_for_cluster_analysis <- data_matrix[ ,samples ]

message( "------ Dimension of data matrix after selecting patients ------" )
print( dim(data_matrix_for_cluster_analysis) )


if( check_gene ){
  if( check_name == "Entrez"){
    gene_index <- which( entrez %in% gene_name$Entrez )
    #print(gene_index)
    data_matrix_for_cluster_analysis <- data_matrix_for_cluster_analysis[gene_index,]
    entrez <- entrez[ gene_index ]
    symbol <- symbol[ gene_index ]
  }
  if( check_name == "Symbol"){
    gene_index <- which( symbol %in% gene_name$Symbol )
    #print(gene_index)
    data_matrix_for_cluster_analysis <- data_matrix_for_cluster_analysis[gene_index,]
    entrez <- entrez[ gene_index ]
    symbol <- symbol[ gene_index ]
  }
}

message( "------ Dimension of data matrix after selecting genes ------" )
print( dim(data_matrix_for_cluster_analysis) )

detectSig <- function( exp.matrix,cluster.num,cluster.assignment,method = "anova",fdr.threshold = 1e-5 ){
  ### description: 
  ### exp.matrix: row(genes) col(patients)
  ### return the representative index for each cluster (a list)
  representive.genes <- list()
  
  gene.num <- dim(exp.matrix)[1]
  patient.num <- dim(exp.matrix)[2]
  ### Option: method aov("anova"), t.test("t"), wilcoxon.test("wilcox"), t&wilcox("t+wilcox")
  avg.exp <- matrix( 0,nrow = gene.num,ncol = cluster.num )
  rownames(avg.exp) <- row.names( exp.matrix )
  for( cluster in 1:cluster.num ){
    temp.index <- which( cluster.assignment == cluster )
    for( i in 1:gene.num ){
      avg.exp[i,cluster] <- mean( exp.matrix[i,temp.index] )
    }
  }
  
  if( method == "t+wilcox" ){
    message( "Method adopted: t+wilcox test" )
    for( cluster in 1:cluster.num ){
      cluster1 <- which( cluster.assignment == cluster )
      gene.index <- seq( gene.num )
      for( comparison in 1:cluster.num ){
        if( comparison != cluster ){
          cluster2 <- which( cluster.assignment == comparison  )
          ### t.test
          t.pval <- array( 0,dim=gene.num )
          names( t.pval ) <- rownames( exp.matrix )
          for( i in 1:gene.num ){
              obj <- try( t.test( exp.matrix[i,cluster1],exp.matrix[i,cluster2] ), silent = TRUE )
              if( is( obj,'try-error' ) ){
                t.pval[i] <- NA
                message( "The gene has happened an t.test-error" )
                print( rownames(exp.matrix)[i] )
              }
              else{
                t.stat <- obj 
                t.pval[i] <- t.stat$p.value
              }
              t.pval.adjust <- p.adjust( t.pval,method = "BH" )
              names(t.pval.adjust) <- rownames( exp.matrix )
          }
              ### wilcox.test
          wilcox.pval <- array( 0,dim=gene.num )
          names( wilcox.pval ) <- rownames( exp.matrix )
              
          for( i in 1:gene.num ){
            wil.stat <- wilcox.test( exp.matrix[i,cluster1],exp.matrix[i,cluster2] )
            wilcox.pval[i] <- wil.stat$p.value
          }
              
          wilcox.pval.adjust <- p.adjust( wilcox.pval,method = "BH" )
          names(wilcox.pval.adjust) <- rownames( exp.matrix )
          
              
          ### output 
          index1 <- which( t.pval.adjust < fdr.threshold & 
                             wilcox.pval.adjust < fdr.threshold & 
                             avg.exp[,cluster] > avg.exp[,comparison] )
          
          gene.index <- intersect( index1,gene.index )
        }
      }
        
      number_of_sigs <- length( genes_ID )
      print( number_of_sigs )
      if( number_of_sigs == 0 ){
        representive.genes <- c(representive.genes,NA)
      }
        
      if( number_of_sigs > 0 ){
        representive.genes <- c( representive.genes,list(gene.index) )
      }
    }
  }
  
  if( method == "anova" ){
    message("Method adopted: anova")
    aov.pvalue <- array( 0,dim = gene.num )
    for( i in 1:gene.num ){
      data_aov <- data.frame(
        exp = exp.matrix[i,],
        group = factor( cluster.assignment )
      )
      data.fit <- lm( exp ~ group, data = data_aov )
      aov.result <- anova( data.fit )
      aov.pvalue[i] <- aov.result$`Pr(>F)`[1]
    }
    aov.pvalue.adjust <- p.adjust( aov.pvalue,method = "BH" )
    aov.index <- which( aov.pvalue.adjust < fdr.threshold )
    exp.index <- apply( avg.exp,1,which.max )
    for( cluster in 1:cluster.num ){
      cluster.exp.index <- which( exp.index == cluster )
      index1 <- intersect( aov.index,cluster.exp.index )
      if( length(index1) > 0 ){
        representive.genes <- c( representive.genes,list(index1) )
      }
      else{
        representive.genes <- c(representive.genes,NA)
      }
    }
  }
  return( representive.genes )
}

###################################Obtain signatures####################################
#### Now suppose there are k clusters###############
#### t.test + wilcoxon.test

gene.num <- dim(data_matrix_for_cluster_analysis)[1]
representative.genes <- detectSig( exp.matrix = data_matrix_for_cluster_analysis, cluster.num = k,cluster.assignment = samples.cluster )

message( "------ Number of signatures in each cluster ------" )
for( i in 1:length(representative.genes) ){
  gene.index <- representative.genes[[ i ]]
  if( is.na(gene.index[1]) ){
    message( "No representative genes were detected in cluster:" )
    print(i)
  }
  else{
    genes_ID <- entrez[gene.index]
    genes_Symbol <- symbol[gene.index]
    output_sig <- data.frame(
      genes_ID = genes_ID,
      genes_Symbol = genes_Symbol
    )
    print(length(gene.index))
    write.table( output_sig,file = paste0(subfolder,output,".",i,".sig"),quote = F,row.names = F,col.names = F )
      ### plot boxplot
    for( i in 1:length(gene.index) ){
      temp.symbol <- output_sig$genes_Symbol[i]
      temp.id <- output_sig$genes_ID[i]
      pdf( paste0( subsubfolder,temp.id,'_',strsplit(temp.symbol,split = "/")[[1]][1],'.pdf' ),height=3,width=5 )
      data_plot <- data.frame(
        exp = data_matrix_for_cluster_analysis[which( rownames(data_matrix_for_cluster_analysis) == temp.id ),],
        cluster = samples.cluster
      )
      p <- ggplot( data_plot,aes(x=factor(cluster),y=exp,fill=colors[cluster]) )
      p2 <- p+geom_boxplot()+geom_jitter( width=.2 )+theme_bw()+theme(legend.position="none")+xlab("Cluster")+ylab("Expression")+ggtitle(temp.symbol)
      print(p2)
      dev.off()
    }
  }
}
 
  






