#' @name Ide_Hub
#' @title Hub genes identification
#'
#' @description This function is able to identify hub genes for key modules
#' @aliases Ide_Hub
#' 
#' @param net A list obtained in \code{\link{Mod_Infer}}
#' @param data_Expr A data.frame generated from \code{\link{Mean_Expre}}
#' @param KeyMod The dataframe of identified key modules obtained in \code{\link{Ide_KeyMod}}
#' @param markers A dataframe of marker genes of each given cell groups
#' @param mod A list generated from \code{\link{Mod_Extra}}
#' @param sftpower The soft-threshold inferred from Mod_Infer
#' 
#' @details If markers is NULL, we will calculate the logFC values between the query cell type and others;
#'          If mod id NULL, this function will generated a list containing members and ME values of each module  
#' @returns a dataframe of hub genes for each key module of each cell group/type
#' @import stringr
#' @import WGCNA
#' @import magrittr
#' @importFrom stats aggregate
#' @importFrom Seurat FindAllMarkers
#' @importFrom dplyr filter
#' @importFrom dplyr top_n
#' @export
Ide_Hub<-function(net,data_Expr,KeyMod,markers = NULL,mod = NULL,sftpower = NULL){
  if(is.null(net)){
    print("aggument 'net' is missing, with no default")
    break
  }
  if(is.null(data_Expr)){
    print("aggument 'data_Expr' is missing, with no default")
    break
  }
  if(is.null(KeyMod)){
    print("aggument 'KeyMod' is missing, with no default")
    break
  }
  if(is.null(mod)){
    mod<-Mod_Extra(data_Expr = data_Expr,net = net,Cyt = FALSE,sftpower = NULL)
  }
  if(is.null(sftpower)){
    sftpower<-net$sftpower
  }
  comp_group<-rownames(KeyMod)
  modProbes<-mod$modProbes
  Hubs<-c()
  for (c in rownames(KeyMod)) {
    module<-KeyMod$Keymodule[rownames(KeyMod) == c]
    module_gene<-modProbes[[module]]
    moduleColors = labels2colors(net$colors)
    ADJ1=abs(cor(data_Expr[,moduleColors==module],use="p"))^sftpower
    Alldegrees1=intramodularConnectivity(ADJ1, rep(module,length(module_gene))) 
    
    modProbes.inter<-intersect(colnames(data_Expr),module_gene)
    
    kWithin<-data.frame(kWithin=Alldegrees1[modProbes.inter,]$kWithin,fromNode=modProbes.inter,row.names=modProbes.inter)
    if(!is.null(markers)){
      hubs<-intersect(markers[which(markers$cluster == c),]$gene,module_gene)
    }else{
      hubs<-c()
    }
    hubs<-intersect(markers[which(markers$cluster == c),]$gene,module_gene)
    if(!length(hubs) == 0){
      k<-top_n(kWithin[hubs,],n=5,wt=kWithin)
      hubs<-rownames(k)
    }else{
      data<-data_Expr[,moduleColors==module]
      if(is.null(comp_group)){
        print("aggument 'comp_group' is missing, with no default")
        break
      }
      data<-cbind(data,comp_group)
      colnames(data)[ncol(data)]<-"group"
      comp_mean<-aggregate(data[,1:(ncol(data)-1)],list(data$group),mean)
      rownames(comp_mean)<-comp_mean[,1]
      comp_mean<-comp_mean[,-1]
      
      other_types<-setdiff(rownames(comp_mean),c)
      other_types_mean<-apply(comp_mean[other_types,],2,mean)
      logFC_mean<-log2(comp_mean[c,]/other_types_mean)
      logFC_mean<-t(logFC_mean)
      logFC_mean<-as.data.frame(logFC_mean)
      colnames(logFC_mean)<-"logFC_mean"
      
      logFC<-cbind(logFC_mean,kWithin[rownames(logFC_mean),"kWithin"])
      colnames(logFC)[2]<-"kWithin"
      logFC<-logFC[order(abs(logFC$logFC_mean),logFC$kWithin,decreasing = TRUE),]
      
      k<-logFC %>% filter(abs(logFC) > 0.5) %>% top_n(5,wt=kWithin)
      
      hubs<-rownames(k)
    }
    
    Hub<-hubs[1]
    for(i in 2:length(hubs)){
      Hub<-c(paste(Hub,hubs[i],sep = ", "))
    }
    Hubs<-c(Hubs,Hub)
  }
  Hubs<-data.frame(Hubs = Hubs,row.names = rownames(KeyMod))
  return(Hubs)
}
