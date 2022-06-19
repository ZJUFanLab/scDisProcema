#' @name CCScore
#' @title Calculate the correlation coefficient score
#'
#' @description  Assessment of the relationship between the trait of cell type and the module activity
#' 
#' @aliases CCScore
#' @param MEs_table A data.frame contains the eigengene values of each module, obtained in \code{\link{Mod_Extra}}
#' @param comp_group The vector contains cell group/type class identity. The length of the vector must be equal to the number of rows of MEs_table
#' @param values The trait values with which to compare, corresponding with comp_group. The value must be numerical and the length of which must be equal to the length of comp_group
#' @param scale Logical. If TRUE, the score will be scaled according to group.

#' @details The order of the comp_group, values and rownames of the MEs_table should correspond to each other
#' @returns a dataframe of CCscores for each cell group/type
#' @import stringr
#' @importFrom stats aggregate
#' @importFrom Hmisc rcorr
#' @export

CCScore<-function(MEs_table,comp_group,values,scale = FALSE){
  if(is.null(MEs_table)){
    print("aggument 'MEs_table' is missing, with no default")
    break
  }
  if(is.null(comp_group)){
    print("aggument 'comp_group' is missing, with no default")
    break
  }
  if(is.null(values)){
    print("aggument 'values' is missing, with no default")
    break
  }
  if(length(comp_group) != length(values)){
    print("Group information does not match the values")
    break
  }
  if(length(comp_group) != nrow(MEs_table)){
    MEs_col<-t(MEs_table)
    if(length(comp_group) != nrow(MEs_table)){
      print("Group information does not match the original ME table")
      break
    }
  }
  MEs_table$group<-comp_group
  
  MEs_cor<-data.frame(row.names = colnames(MEs_table)[1:(ncol(MEs_table)-1)])
  for (c in unique(comp_group)) {
    MEs_table_c<-subset(MEs_table,MEs_table$group==c)
    cell_ratio_c<-values[which(comp_group == c)]
    
    cor<-""
    pvalue<-""
    for (i in 1:(ncol(MEs_table_c)-1)) {
      MEc<-rcorr(MEs_table_c[,i],cell_ratio_c,type = "pearson")
      cor[i]<-MEc$r[1,2]
    }
    MEs_cor<-cbind(MEs_cor,cor)
    colnames(MEs_cor)[ncol(MEs_cor)]<-c
  }
  
  score<-t(MEs_cor)
  score<-as.data.frame(apply(score,2,as.numeric),row.names = rownames(score))
  
  if(scale){
    scale_score<-scale(t(score))
    scale_score<-data.frame(t(scale_score),row.names = rownames(score))
    return(scale_score)
  }else{
    return(score)
  }
}