#' @name Ide_Hub
#' @title Mean expression of features in given groups
#' 
#' @description Process the Seurat Object. According to the group information provided by users, cells were divided into cell states and the average expression value of features(all the features of Seurat object or provided by the user) in each cell state was obtained for subsequent gene module inference
#' 
#' @param object A Seurat object.
#' @param group.1 The group information provided by users, which must be an identity class of the Seurat object
#' @param group.2 The second group information provided by users, which must be an identity class of the Seurat object
#' @param group.3 The third group information provided by users, which must be an identity class of the Seurat object
#' @param features Genes to used.Default is to use all genes
#' @param slot Slot to pull data from..Default slot will be "counts"
#' @details The order of groups is interchangeable, but group.3 cannot be set if group.2 is Null. The as_matrix function was used to transform a sparse matrix into an ordinary one in a faster and more memory efficient manner
#' @return a data.frame of mean expression of features in each cell states
#' @import stringr
#' @import Seurat
#' @importFrom stats aggregate
#' @importFrom dplyr setdiff
#' @export
Mean_Expre<-function(object,group.1,group.2 = NULL,group.3 = NULL,features = NULL,slot = "count"){
  if(is.null(object)){
    print("argument 'object' is missing, with no default")
    break
  }
  if(is.null(group.1)){
    print("There is no group designated")
    break
  }
  if(is.null(features)){
    features = rownames(object)
  }
  
  count<-t(Seurat::GetAssayData(object = object, 
                        slot = slot)[features,])
  
  count<-as.matrix(count)
  
  if(is.null(group.2)){
    group.2<-group.1
    if(is.null(group.3)){
      group.3<-group.1
    }else{
      print("argument 'group.2' is skipped")
      break
    }
  }else{
    if(is.null(group.3)){
      group.3<-group.2
    }
  }
  
  group<-data.frame(group1 = object@meta.data[,group.1],
                    group2 = object@meta.data[,group.2],
                    group3 = object@meta.data[,group.3],
                    row.names = colnames(object))
  
  if(group.2 == group.1){
    if(group.3 == group.2){
      g<-group$group1
    }
  }else{
    if(group.3 == group.2){
      g<-paste(group$group1,group$group2,sep = "_")
    }else{
      g<-paste(group$group1,group$group2,group$group3,sep = "_")
    }
  }
  
  group$group<-g
  
  object.mean<-data.frame()
  for (g1 in as.character(unique(object@meta.data[,group.1]))) {
    for (g2 in as.character(unique(object@meta.data[,group.2]))) {
      for (g3 in as.character(unique(object@meta.data[,group.3]))) {
        
        rownames<-rownames(group)[which(group$group1 == g1 & group$group2 == g2 & group$group3 == g3)]  
        if(length(rownames) == 0){
          next
        }
        mat<-count[rownames,]
        mat<-data.frame(group=group[rownames,]$group,mat)
        mean<-aggregate(mat[,2:ncol(mat)],list(mat$group),mean)   
        rownames(mean)<-mean$Group.1
        mean<-mean[,-1]
        #rm(count)
        #gc()
        object.mean<-rbind(object.mean,mean)
      }
    }
  }
  gc()
  
  #Put a dot in the line name and change it to a minus
  features_diff<-setdiff(colnames(object.mean),features)
  colnames(object.mean)[colnames(object.mean) %in% features_diff]<-str_replace_all(features_diff,"\\.","-")
  return(object.mean)
}
