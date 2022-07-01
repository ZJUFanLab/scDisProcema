#  Tutorial
  The R package __scDisProcema__ integrates Weighted Gene Co-Expression Network Analysis (WGCNA) pipline for tracing the dynamic change of gene expression patterns during disease progression using single-cell data.
  
  The workflow of the scDisProcema can be divided into three steps: 1) gene-state matrix as input to WGCNA for gene co-expression network module searching; 2) module dynamic analysis for degree of dynamic change and relevance to cells over the course of disease; 3) key module identification for different cell types.

## General usage
  The input of the scDisProcema is a Cell state-Genes matrix. The sample data can be downloaded [here](/data/mean.csv).
  
    #read the file
    mean<-read.csv("mean.csv",row.names = 1)
    
    #check the file
    mean[1:4,1:4]
    # >        IGLV3.19    IGHV4.34 IGLV3.32     JCHAIN
    # H_B      0.897019627 1.711897262        0 7.54906712
    # H_CD4+ T 0.001269036 0.005076142        0 0.01895622
    # H_CD8+ T 0.000819104 0.002079264        0 0.10282906
    # H_mDC    0.096563011 0.001636661        0 0.05891980
    
    #cell states are obtained from 9 cell types of 5 disease states
    rownames(mean)
    # >[1] "H_B"           "H_CD4+ T"      "H_CD8+ T"      "H_mDC"        
    # [5] "H_Mega"        "H_Mono/Macro"  "H_NK"          "H_pDC"        
    # [9] "H_γδ T"        "MP_B"          "MP_CD4+ T"     "MP_CD8+ T"    
    # [13] "MP_mDC"        "MP_Mega"       "MP_Mono/Macro" "MP_NK"        
    # [17] "MP_pDC"        "MP_γδ T"       "SP_B"          "SP_CD4+ T"    
    # [21] "SP_CD8+ T"     "SP_mDC"        "SP_Mega"       "SP_Mono/Macro"
    # [25] "SP_NK"         "SP_pDC"        "SP_γδ T"       "SC_B"         
    # [29] "SC_CD4+ T"     "SC_CD8+ T"     "SC_mDC"        "SC_Mega"      
    # [33] "SC_Mono/Macro" "SC_NK"         "SC_pDC"        "SC_γδ T"      
    # [37] "MC_B"          "MC_CD4+ T"     "MC_CD8+ T"     "MC_mDC"       
    # [41] "MC_Mega"       "MC_Mono/Macro" "MC_NK"         "MC_pDC"       
    # [45] "MC_γδ T" 
    
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
 
 
  The net list contains information of color indexes, members and module eigengene (ME) values of all inferred modules. The Mod_Extra function can be convert these information to a format that is easy to use later.
      
      #If Cyt = TRUE, the function will export files available for Cytoscape to visualize the network
      #sftpower: the power of the weight between genes. The default value is the soft threshold obtained in Mod_Infer
      
      mod<-Mod_Extra(data_Expr = data_Expr,net = net,Cyt = TRUE,sftpower = net$sftpower)
      
      #the ME values of each module
      MEs_table <- mod$MEs_table
      
      #the members of each module
      modProbes<-mod$modProbes

  Fourth, calculate the module significant score (S(MS)). This score comprehensively measures the correlation between the dynamic expression change of the modules and the dynamic change of the cell type during disease development. Here, we use the change of cell propotion to show the dynamic change of the cell types. The cell propotion file can be downloaded [here](/data/cell_ratio_total.csv).
      
      #read process the cell propotion file
      cell_propotion<-read.csv("cell_ratio_total.csv",row.names = 1)
      #Let the order of cell types in cell_propotion match that of MEs_table
      rownames(cell_propotion)<-paste(cell_propotion$disease_state,cell_propotion$cell_type,sep = "_")
      cell_propotion<-cell_propotion[rownames(MEs_table),]
      
      #calculate S(MS)
      sms<-MSScore(MEs_table = MEs_table,comp_group = cell_propotion$cell_type,values = cell_propotion$cell_ratio)
      
      #visualization
      New_score_scale_round2<-round(sms,2)
      NMF::aheatmap(New_score_scale_round2,filename = "heatmap-New_score_scale.pdf",
              Colv=NA,Rowv=NA,cexRow = 1,cexCol = 1,width =10,height = 5,
              txt = as.matrix(apply(New_score_scale_round2,2,as.character)),
              color = c("#332288","#88CCEE","white","#DDCC77","#882255"))

      
      

