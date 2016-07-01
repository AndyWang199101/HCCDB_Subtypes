# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 18:14:42 2016

@author: wdf
"""

import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

#### define gene sets for every class
s1 = "ACP5,ACTA2,ADAM15,ADAM8,ADAM9,AEBP1,AIF1,AKT3,ALDOA,ALOX5AP,ANXA1,ANXA4,ANXA5,AP2S1,AQP1,ARF5,ARHGDIB,ARPC1B,ARPC2,ASAH1,ATP1B3,ATP6AP1,ATP6V0B,ATP6V1B2,ATP6V1F,BCL2A1,BLVRA,C1QB,C3AR1,CAPZA1,CBFB,CCL3,CCL5,CCND2,CCR7,CD151,CD37,CD3D,CD47,CD48,CD53,CD74,CD8A,CDC20,CDC25B,CDH11,CDK2AP1,CELF2,CHN1,COL11A1,COL15A1,COL1A2,COL3A1,COL4A1,COL4A2,COL5A2,COL6A1,COL6A2,CORO1A,CRABP2,CRIP2,CSRP3,CTGF,CTSC,CTSS,CXCL1,CXCR4,CYBA,CYBB,CYFIP1,CYP1B1,CYR61,DAB2,DCTN2,DDR1,DDR2,DDX11,DGKA,DGKZ,DNM2,DPYSL2,DUSP5,DUT,EFEMP1,EFNB1,F13A1,FBN1,FBRS,FCGBP,FCGR2A,FGL2,FHL3,FLNA,FUT4,FYB,GEM,GLIPR1,GNAI2,GNS,GPNMB,GRN,GSTP1,GUCY1A3,GYPC,HCLS1,HEXA,HIF1A,HK1,HLA-DMA,HLA-DOA,HLA-DPB1,HLA-DQA2,HLA-DQB1,ID3,IER3,IFI16,IFI30,IGFBP5,IGKC,IGLL1,IKBKE,IL15RA,IL2RB,IL2RG,IL7R,IQGAP1,IRF1,ITGB2,ITPR3,KIAA0101,KLC1,KLF5,LAMB1,LAPTM5,LCP1,LDHB,LGALS1,LGALS3BP,LGALS9,LGMN,LHFPL2,LITAF,LMO4,LOX,LSP1,LTBP2,LTBP3,LTF,LUM,LYN,M6PR,MAP1B,ME1,MFAP1,MGP,MPHOSPH6,MSN,MTHFD2,MYCBP2,NBL1,NPC2,NSMAF,OAZ1,PAK1,PAM,PAPSS1,PEA15,PFN1,PGK1,PIM2,PKD2,PKMYT1,PKN1,PLAUR,PLD3,PLP2,PNP,POLD3,POSTN,PPIC,PPP1CB,PPP4C,PPP4R1,PRKD2,PRMT5,PROCR,PSMD2,PTPRC,PYGB,QSOX1,RAB31,RALB,RALGDS,RCC1,RHOA,RIN2,RIT1,RNASE1,RNASE6,RPA2,RSU1,S100A10,S100A11,S100A13,SLA,SLBP,SLC1A5,SLC2A1,SLC2A5,SLC39A6,SLC7A5,SMAD2,SMARCD1,SPAG8,SRGN,SRI,STK38,STX3,TAGLN,TAX1BP3,TCF4,TFF3,THY1,TIMP2,TMSB4X,TNFRSF1B,TP53BP1,TPM2,TRA@,TRAF3,TRAF5,TRIP10,TSPAN3,TUBA4A,VCAN,ZNF384"
s2 = "ABCB10,ABCD3,ACP1,ADD3,AFP,AHCY,ARHGAP35,ARID3A,ATF2,ATM,ATP2B1,ATP2B2,ATP5F1,ATXN10,BCAM,BCLAF1,BRD3,BTG3,C5orf13,CASC3,CD46,CDK6,CHKA,CLK2,COL2A1,CPD,CSE1L,CSNK2A1,CSNK2A2,CTNNB1,CUL4A,CXADR,DDX1,DDX18,DEK,EIF4A2,EIF4B,ENPP1,EP300,ERBB3,FBL,FGFR3,FGFR4,FLNB,GBF1,GCN1L1,GLUD1,GNAI1,GPC3,GTF2I,GTF3C2,H1F0,HELZ,HMGCR,HNRNPA2B1,HNRNPC,HNRNPU,IDI1,IGF2,IGF2R,ITIH2,KLF3,LBR,MAPK6,MEST,NCOA4,NET1,NR2C1,NR5A2,NT5E,NUP153,PEG3,PHF3,PHKA2,PIGC,PLXNB1,PNN,POFUT1,PPARG,PPP2R1A,PRDX3,PTOV1,RAB4A,RBM39,RPL24,RPL27,RPL31,RPS19,RPS24,RPS25,RPS27,RPS5,RRP1B,SEPHS1,SLC6A2,SLC6A5,SMARCA1,SMARCC1,SNRPE,SNTB1,SREBF2,SSB,SUMO1,SUZ12,TARBP1,TBCE,TFIP11,TIA1,TIAL1,TM9SF4,TP53BP2,TPR,TRIM26,TTC3,UBE2K"
s3 = "ABCB4,ABCC6,ABHD2,ABP1,ACAA2,ACADM,ACADS,ACADSB,ACADVL,ACO1,ACOX1,ACOX2,ACSL1,ACY1,ADH4,ADH6,ADK,AGL,AGXT,AKR1C1,ALAS1,ALDH1A1,ALDH1B1,ALDH2,ALDH3A2,ALDH4A1,ALDH6A1,ALDH7A1,ALDOB,ALPL,AMFR,AMT,ANXA6,APCS,APOA1,APOC2,APOC4,APOH,AQP7,ARG1,ARHGEF12,ARSA,ASCL1,ASGR1,ASGR2,ASL,ASS1,ATOX1,ATP5D,ATP5J,AZGP1,BAAT,BDH1,BHMT,BLOC1S1,BLVRB,BPHL,BTD,C1R,C1S,C4A,C4BPA,C8B,CA2,CAT,CBR1,CD14,CD302,CD81,CECR1,CES1,CFB,CFH,CGREF1,CNGA1,COL18A1,COX5B,CP,CPA3,CPA4,CPB2,CPS1,CRABP1,CRYAA,CRYM,CSTB,CTH,CTSO,CXCL2,CYB5A,CYFIP2,CYP21A2,CYP27A1,CYP2C9,CYP2J2,CYP3A7,DAO,DCAF8,DECR1,DNASE1L3,DPAGT1,DRG2,ECHS1,ECI1,EDNRB,EGFR,EHHADH,EMP2,EPAS1,EPHX1,ETS2,F11,F2,F5,FAH,FANCA,FGB,FGG,FH,FKBP2,FLT4,FMO4,FOXO1,FXR2,GCH1,GCHFR,GCKR,GGH,GHR,GJB1,GLYAT,GOT2,GPT,GPX2,GPX3,GSTA2,GSTO1,GSTZ1,HAAO,HADH,HGD,HMGCS2,HMOX2,HPD,HRG,HRSP12,HSD17B10,HSD17B4,ICAM3,IDH2,IDH3A,IFIT1,IGF1,IL13RA1,IL32,IL6R,IMPA1,INSR,IQGAP2,ISG15,ITIH1,ITIH3,ITIH4,ITPR2,IVD,KCNJ8,KLKB1,KMO,KNG1,LCAT,LONP1,LPIN1,LPIN2,MAOA,MAOB,MAPRE3,MGST2,MME,MSMO1,MT2A,MTHFD1,MTHFS,MUT,MYLK,MYO1E,NDUFV2,NFIB,NFIC,NFKBIA,NNMT,NRG1,PAH,PAPSS2,PCCA,PCCB,PCK1,PCK2,PDK4,PGM1,PGRMC1,PIK3R1,PKLR,PLA2G2A,PLCG2,PLG,PLGLB2,PNPLA4,POLD4,PON3,PPP2R1B,PROS1,PTGR1,PTS,QDPR,RARRES2,RBP5,RGN,RHOB,RNASE4,SBDS,SDC1,SDHB,SDS,SELENBP1,SEPP1,SERPINA3,SERPINA6,SERPINC1,SERPING1,SHB,SHMT1,SLC10A1,SLC16A2,SLC23A1,SLC23A2,SLC2A2,SLC35D1,SLC6A1,SLC6A12,SLC7A2,SLC9A3R2,SLCO2A1,SLPI,SMARCA2,SOAT1,SOD1,SOD2,SORL1,SPAM1,SPARCL1,SRD5A1,SREBF1,SULT2A1,TCEA2,TDO2,TGFBR3,TINAGL1,TJP2,TMBIM6,TMOD1,TOB1,TPMT,TST,UQCRB,VSIG2,ZNF160"

s1 = s1.split(",")
s2 = s2.split(",")
s3 = s3.split(",")

print len(s1),len(s2),len(s3)

#gene_file = open( "../gene_train.txt" )
gene_file = open( "../genes.txt" )
genes_TCGA = [line.rstrip() for line in gene_file] 
gene_file.close()

s1 = [gene for gene in s1 if gene in genes_TCGA]
s2 = [gene for gene in s2 if gene in genes_TCGA]
s3 = [gene for gene in s3 if gene in genes_TCGA]



print len(s1),len(s2),len(s3)

### construct template
number_in_s1 = len(s1)
number_in_s2 = len(s2)
number_in_s3 = len(s3)

s1_template = np.zeros( number_in_s1+number_in_s2+number_in_s3 ) -1
s1_template[0:number_in_s1] = 1

s2_template = np.zeros( number_in_s1+number_in_s2+number_in_s3 )-1
s2_template[number_in_s1:(number_in_s1+number_in_s2)] = 1

s3_template = np.zeros( number_in_s1+number_in_s2+number_in_s3 )-1
s3_template[ (number_in_s1+number_in_s2): ] = 1

### compute the nearest template
file = open( "../genomic_tumor.txt" )
#head1 = file.readline()
#head2 = file.readline()

genomic = {}
head = file.readline().rstrip().split("\t")[1:]
for patient in head:
    genomic[patient] ={}
while 1:
    line = file.readline()
    if not line:
        break
    temp = line.rstrip().split("\t")
    gene  = temp[0]
    genes_expr = np.array( [float(x) for x in temp[1:]] )
    genes_expr = ( genes_expr - genes_expr.mean() ) / (genes_expr.std())
    
    if gene in s1 or gene in s2 or gene in s3:
        for i in range( len(temp)-1 ):
            genomic[ head[i] ][ gene ] = float( genes_expr[i] )

out = open( "S1S2S3.cluster","w" )
out2 = open( "TCGA_S123_expr.txt","w" )
out.write( "\t".join( ["Sample","Cluster"] )+"\n" )
out2.write( "Sample\tCluster\t"+"\t".join( s1 )+"\t"+"\t".join( s2 )+"\t"+"\t".join(s3)+"\n" )
out2.write( "Label\tNA\t"+"\t".join( [str(1) for dummy in s1] )+"\t"+"\t".join( [str(2) for dummy in s2] )+"\t"+"\t".join( [str(3) for dummy in s3] )+"\n" )
for patient in head: 
    out_line = [patient]
    out2_line = [patient]
    expr = []
    for gene in s1:
        expr.append( genomic[patient][gene] )
        
    for gene in s2:
        expr.append( genomic[patient][gene] )
        
    for gene in s3:
        expr.append( genomic[patient][gene] )
     
    expr = np.array( expr )
    #print expr
    cos1 = cosine_similarity( expr,s1_template )
    cos2 = cosine_similarity( expr,s2_template )
    cos3 = cosine_similarity( expr,s3_template )
    #print cos1,cos2,cos3
    if cos1 > cos2 and cos1 > cos3:
        out2_line.append( str(1) )
        out_line.append(str(1))
    elif cos2 > cos1 and cos2 > cos3:
        out_line.append(str(2))
        out2_line.append( str(2) )
    elif cos3 > cos1 and cos3 > cos2:
        out_line.append( str(3) )
        out2_line.append( str(3) )
    else:
        out_line.append( "NA" )
    for e in expr:
        out2_line.append( str(e) )
    out.write( "\t".join( out_line )+"\n" )
    out2.write( "\t".join( out2_line )+"\n" )
    
out2.close()

out.close()
file.close()




