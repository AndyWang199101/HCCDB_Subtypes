# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 19:57:41 2016

@author: wdf

Latest Modified: 2016-06-04 BY WDF

Modified: 
(1). update the storage of results
???? (2). Use 5,000 most variated genes (overall 19600) as backgroud gene for hypergeometric test

Description: This script perform the following process
    Given signature results of s2, construct the signatures-similarity network whose:
        (1). Nodes: every signature set of each single cluster in each dataset (number of genes > 0)
        (2). Edges: weighed by either Jaccard similarity of two nodes or statistical significance of two nodes' overlap
    

Parameters:
    e.g. s3.para
    1. Signature files: one filename in which list all gene sigs. 
    2. Threshold: minimum size of nodes
    3. VERSION: as s1 & s2
    4. background gene set used to perform hypergeometric test

Output:
    1. node_file.txt with cols: ["Node","Dataset","Cluster_single","Size"]
    2. edge_file.txt with cols: ["Node1","Node2","weight_jaccard_index","weight_hyper_geometric_log_p"]

"""

from math import log
import scipy.stats as stats

import sys

sig_files = sys.argv[1]
node_threshold = int(sys.argv[2])
edge_threshold = int(sys.argv[3])
prefix = sys.argv[4]
background = sys.argv[5]



node_file_name = prefix+'node_file_'+str(node_threshold)+'_'+str(edge_threshold)+'.txt'
edge_file_name = prefix+'edge_file_'+str(node_threshold)+"_"+str(edge_threshold)+'.txt'
#def _hyper_geometric_log_pvalue( N,D1,D2 ):
#    intersect = len( set(D1).intersection( set(D2) ) )
#    cumulative_p = 0
#    for k in range( intersect+1,min( len(D1),len(D2) ) ):
#        temp1 = len( [i for i in it.combinations(D1,k) ] ) 
#        temp2 = len( [i for i in it.combinations( N-D1,len(D2)-k ) ] )
#        temp3 = len( [i for i in it.combinations( N,len(D2) ) ] )
#        current_p = exp( log( len(temp1) ) + log( len(temp2) )- log( len(temp3) ) )
#        cumulative_p += current_p
#    return log(cumulative_p)
     

files_file = open( sig_files )
files = [line.rstrip() for line in files_file]
files_file.close()

node_file = open( node_file_name ,"w" )
head = ["Node","Dataset","Cluster_single","Size"]
node_file.write( '\t'.join(head)+'\n' )

nodes = {}
for file in files:
    temp = file.split('/')[-1]
    dataset,cluster,sig = temp.split('.')
    fi = open( file )
    genes = [line.rstrip().split()[1] for line in fi]
    #print genes
    size = len(genes)
    if size >= node_threshold:   ######################################################
        nodes[temp] = genes
        node_file.write( '\t'.join([temp,dataset,cluster,str(size)])+'\n' )
    fi.close()
    
node_file.close()

edge_file = open( edge_file_name,"w" )
head = ["Node1","Node2","weight_jaccard_index","weight_hyper_geometric_log_p"]
edge_file.write( '\t'.join( head )+'\n' )

nodes_list = nodes.keys()

hsa_file = open( background )
all_genes = [ line.rstrip().split()[0] for line in hsa_file ]
all_genes = set( all_genes )

for i in range(len(nodes_list)-1):
    for j in range( i+1,len(nodes_list) ):
        Node1 = nodes_list[i]
        Node2 = nodes_list[j]
        genes1 = set( nodes[Node1] )
        genes2 = set( nodes[Node2] )
        jaccard_sim = 0
        if len(genes1)>0 and len(genes2)>0:
            inter = genes1.intersection(genes2)
            uni = genes1.union(genes2)
            jaccard_sim = float( len(inter) ) / float( len(uni) )
            #print len(inter),len(genes1),len(genes2)
            if (len(inter) == len(genes1) or len(inter) == len(genes2)):
                continue
            hp = stats.hypergeom.sf( len(inter),len(all_genes),len(genes1),len(genes2)  )
            hp = -log( hp,10 )
            if jaccard_sim > 0 and hp > edge_threshold:
                line = [Node1,Node2,str(jaccard_sim),str(hp)]
                edge_file.write( '\t'.join( line ) +"\n" )

edge_file.close()
