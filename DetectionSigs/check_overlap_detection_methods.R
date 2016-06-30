##### Latest Modified: 0630 BY WDF

##### Description: draw venn diagram to analyse the overlap of different detection method

args <- commandArgs(TRUE)
options(stringsAsFactors = FALSE)


##### Parameters:

##### 1. prefix1 (define the first set) e.g. "/data/home/dwang/HCCDB_Subtypes/results/v0616_1vall_hsa05200_sigs/"
##### 2. prefix2 (define the second set) e.g. "/data/home/dwang/HCCDB_Subtypes/results/v0616_aov_hsa05200_sigs/"
##### 3. prefix3 (define the third set) e.g. "/data/home/dwang/HCCDB_Subtypes/results/v0616_1v1_hsa05200_sigs/"

prefix1 <- args[1]
prefix2 <- args[2]
prefix3 <- args[3]

label1 <- strsplit( prefix1,"/" )[[1]]
label1 <- label1[ length(label1) ]
temp1 <- strsplit( label1,split = "_" )[[1]]
method1 <- temp1[2]
geneset <- temp1[3]

label2 <- strsplit( prefix2,"/" )[[1]]
label2 <- label2[ length(label2) ]
temp2 <- strsplit( label2,split = "_" )[[1]]
method2 <- temp2[2]

label3 <- strsplit( prefix3,"/" )[[1]]
label3 <- label3[ length(label3) ]
temp3 <- strsplit( label3,split = "_" )[[1]]
method3 <- temp3[2]

### define output folder
folder <- paste0( "~/HCCDB_Subtypes/results/detection_sig_overlap/",geneset,"/" )
system( paste0( "mkdir ",folder ) )

library("VennDiagram")

sig <- read.table( "/data/home/dwang/HCCDB_Subtypes/parameters/sig.file.0630",header = F )
sig <- sig$V1

data.sig <- c()
data.method1 <- c()
data.method2 <- c()
data.overlap <- c() ### format: num1;num2;overlap

colors <- c("cornflowerblue", "green", "yellow", "darkorchid1") 

for( temp in sig ){
  #data.sig <- c(data.sig,temp)
  
  filename <- strsplit( temp,split = "/" )[[1]][2]
  
  file1 <- paste0( prefix1,temp )
  check.file1 <- file.exists( file1 )
  
  file2 <- paste0( prefix2,temp )
  check.file2 <- file.exists( file2 )
  
  file3 <- paste0( prefix3,temp )
  check.file3 <- file.exists( file3 )
  
  count <- 1
  gene.list <- list()
  if( check.file1 ){
    sig.try <- try( read.table( file1,header = F ) )
    if( !inherits(sig.try,"try-error") ){
      sig.file <- sig.try
      gene.list <- c( gene.list,oneVSothers = list( sig.file$V1 ) )
      count <- count + 1
    }
  }
  
  if( check.file2 ){
    sig.try <- try( read.table( file2,header = F ) )
    if( !inherits(sig.try,"try-error") ){
      sig.file <- sig.try
      gene.list <- c( gene.list,ANOVA = list( sig.file$V1 ) )
      count <- count + 1
    }
  }
  
  if( check.file3 ){
    sig.try <- try( read.table( file3,header = F ) )
    if( !inherits(sig.try,"try-error") ){
      sig.file <- sig.try
      gene.list <- c( gene.list,oneVSone = list( sig.file$V1 ) )
      count <- count + 1
    }
  }
  if( count > 1 ){
    venn.diagram( gene.list,paste0( folder,filename,".tiff" ),col = "transparent", fill = colors[1:length(gene.list)], alpha = 0.50)
  }
}
