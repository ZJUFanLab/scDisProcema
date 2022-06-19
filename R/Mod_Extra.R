#' @name Mod_Extra
#' @title Extract the information of inferred modules
#'
#' @description This function is able to extract the eigengene values and the members of the inferred gene modules. Besides, network list containing nodes, edges and degrees in a format suitable for importing to Cytoscape would also be exported if users want to
#' 
#' @param data_Expr A data.frame generated from Mean_Expre
#' @param net A list obtained in Mod_Infer
#' @param modules The modules users want to extract. Default is NULL, which indicates all of the modules
#' @param Cyt Logical. If TRUE, list files suitable for Cytoscape will be generated in the default path
#' @param sftpower The soft-threshold inferred from Mod_Infer
#' @details The edge weight threshold was set to 0.5^sftpower
#' @returns a list contains the following components:
#' @returns modProbes A list of members of each module
#' @returns MEs_table A dataframe of ME values of each module
#' @returns CytoscapeInput-edges-module_all.csv A dataframe of gene network of each module, suitable for Cytoscape. Nodes of the network can be seen in "fromNode" And "toNode"; "weight" represents the edge weight; "fromAltName" and "fromAltName" are the optional alternate names for the nodes; kWithin represents the intermodule connectiviey degree of each genes
#' @import WGCNA 
#' @import stringr
#' @import magrittr
#' @importFrom utils write.csv
#' @importFrom dplyr setdiff
#' @export
Mod_Extra<-function(data_Expr,net, modules = NULL,Cyt = TRUE,sftpower = NULL){
  if(is.null(data_Expr)){
    print("aggument 'data_Expr' is missing, with no default")
    break
  }
  if(is.null(net)){
    print("aggument 'net' is missing, with no default")
    break
  }
  if(is.null(sftpower)){
     sftpower<-net$sftpower
  }
  probes<-colnames(data_Expr)
  inModule<-""
  modProbes<-list()
  moduleColors = labels2colors(net$colors)
  MEs<-net$MEs
  MEs_col = MEs
  colnames(MEs_col) = paste0("ME", labels2colors(
    as.numeric(str_replace_all(colnames(MEs),"ME",""))))
  MEs_col = orderMEs(MEs_col)
  
  if(is.null(modules)){
    for (m in 1:length(unique(moduleColors))) {
      module = unique(moduleColors)[m];
      inModule = (moduleColors==module);
      modProbes[[module]] = probes[inModule]; 
    }
    MEs_col<-MEs_col
  }else{
    for (m in 1:length(unique(modules))) {
      module = unique(moduleColors)[m];
      inModule = (moduleColors==module);
      modProbes[[module]] = probes[inModule]; 
    }
    MEs_col<-MEs_col[,paste0("ME",unique(modules))]
  }
  if(Cyt){
    if(is.null(sftpower)){
      sftpower<-net$sftpower
    }
    for (module in names(modProbes)) {
      module_gene<-modProbes[[module]]
      ADJ1=abs(cor(data_Expr[,moduleColors==module],use="p"))^sftpower
      Alldegrees1=intramodularConnectivity(ADJ1, rep(module,length(module_gene)))
      
      modProbes.inter<-intersect(colnames(data_Expr),module_gene)
      inModule = (moduleColors==module)
      modTOM = ADJ1
      dimnames(modTOM) = list(modProbes.inter, modProbes.inter)
      modTOM = ADJ1
      
      cyt = exportNetworkToCytoscape(
        modTOM,
        edgeFile = NULL,
        nodeFile = NULL,
        weighted = TRUE,
        threshold = 0.5^sftpower,
        nodeNames = modProbes.inter, 
        nodeAttr = moduleColors[inModule]
      )
      
      kWithin<-data.frame(kWithin=Alldegrees1[modProbes.inter,]$kWithin,fromNode=modProbes.inter)
      Cyto<-cyt$edgeData
      diff<-setdiff(Cyto$toNode,Cyto$fromNode)
      diff_link<-Cyto[Cyto$toNode %in% diff,]
      diff_link$fromNode<-Cyto[Cyto$toNode %in% diff,][,2]
      diff_link$toNode<-Cyto[Cyto$toNode %in% diff,][,1]
      Cyto<-rbind(Cyto,diff_link)
      Cyto<-merge(Cyto,kWithin,by="fromNode")
      
      fromNode_diff<-setdiff(Cyto$fromNode,colnames(data_Expr))
      Cyto$fromNode[Cyto$fromNode %in% fromNode_diff]<-str_replace_all(Cyto$fromNode[Cyto$fromNode %in% fromNode_diff],"\\.","-")
      toNode_diff<-setdiff(Cyto$toNode,colnames(data_Expr))
      Cyto$toNode[Cyto$toNode %in% toNode_diff]<-str_replace_all(Cyto$toNode[Cyto$toNode %in% toNode_diff],"\\.","-")
      write.csv(Cyto,paste("CytoscapeInput-edges-",module,"_all.csv",sep = ""),row.names = F)
      
    }
  }
  mod<-list(modProbes,MEs_col)
  names(mod)<-c("modProbes","MEs_table")
  return(mod)
}

  