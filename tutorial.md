#  Tutorial
  The R package __scDisProcema__ integrates Weighted Gene Co-Expression Network Analysis (WGCNA) pipline for tracing the dynamic change of gene expression patterns during disease progression using single-cell data.
  
  The workflow of the scDisProcema can be divided into three steps: 1) gene-state matrix as input to WGCNA for gene co-expression network module searching; 2) module dynamic analysis for degree of dynamic change and relevance to cells over the course of disease; 3) key module identification for different cell types.

## General usage
  The input of the scDisProcema is a Cell state-Genes matrix. The sample data can be obtained [here](/data/mean.csv).
  
    #read the file
    mean<-read.csv("mean.csv",row.names = 1)
    
    #check the file
    mean[1:4,1:4]
    # >        IGLV3.19    IGHV4.34 IGLV3.32     JCHAIN
    # H_B      0.897019627 1.711897262        0 7.54906712
    # H_CD4+ T 0.001269036 0.005076142        0 0.01895622
    # H_CD8+ T 0.000819104 0.002079264        0 0.10282906
    # H_mDC    0.096563011 0.001636661        0 0.05891980
    
    #replace the "." to "-" in the gene names
    colnames(mean)<-str_replace_all(colnames(mean),"\\.","-")

  Then filter the data:
      
      #if mad = TRUE, screen genes with median absolute deviation (MAD) larger than the threshold.
      #if goodsamples = TRUE, goodSamplesGenes function will be performed to screen out missing entries and genes with zero-variance
      
      data_Expr<-filter_data(mean,mad = TRUE,mad.thre = 0.01,goodsamples = TRUE)
    
  Third, construct gene co-expression network (GCN) and infer gene modules by Mod_Infer(). 
    
      #RsquaredCut: the scale-free topology fitting index R2, set to pick an appropriate soft-thresholding power for scale-free network construction
      #maxBlockSize and minModuleSize: the maximum block size and minimum module size for module detection
      #deepSplit: split sensitivity of the shear tree, ranging from 1 to 4. The greater, the more sensitive the shear tree
      #mergeCutHeight: dendrogram cut height for module merging. The greater, the less modules can be detected
      
      net<-Mod_Infer(data_Expr = data_Expr,RsquaredCut = 0.9,maxBlockSize = "total",
               deepSplit = 4,minModuleSize = 20,mergeCutHeight = 0.1)
 
 
  The net list contains information of color indexes, members and module eigengene values of all inferred modules. The Mod_Extra function can be convert these information to a format that is easy to use later.
