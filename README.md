# HCCDB_Subtypes

Based on our HCCDB database, which collects public available HCC gene expression datasets, we perform a signature network clustering-based methods, to integrate all these datasets to make HCC subtyping.

The datasets we used currently:

/data/home/dwang/HCCDB_Subtypes/DATASETS/HCCDB1/mRNA/HCCDB_hcc.txt
/data/home/dwang/HCCDB_Subtypes/DATASETS/HCCDB11/HCCDB_hcc.txt
/data/home/dwang/HCCDB_Subtypes/DATASETS/HCCDB12/HCCDB_hcc.txt
/data/home/dwang/HCCDB_Subtypes/DATASETS/HCCDB13/HCCDB_hcc.txt
/data/home/dwang/HCCDB_Subtypes/DATASETS/HCCDB14/HCCDB_hcc.txt
/data/home/dwang/HCCDB_Subtypes/DATASETS/HCCDB15/mRNA/HCCDB_hcc.txt
/data/home/dwang/HCCDB_Subtypes/DATASETS/HCCDB3/HCCDB_hcc.txt
/data/home/dwang/HCCDB_Subtypes/DATASETS/HCCDB4/HCCDB_hcc.txt
/data/home/dwang/HCCDB_Subtypes/DATASETS/HCCDB6/HCCDB_hcc.txt
/data/home/dwang/HCCDB_Subtypes/DATASETS/HCCDB8/mRNA/HCCDB_hcc.txt
/data/home/dwang/HCCDB_Subtypes/DATASETS/HCCDB9/HCCDB_hcc.txt

The main steps contain:

(1). Single datasets clustering and signature detection

(2). Signature network construction and clustering

(3). Re-classify single datasets based on consensus signature classes
