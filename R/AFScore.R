#' @name AFScore
#' @title Calculate the correlation activity fluctuation score
#'
#' @description  This function can assess the degree to which the module's activity fluctuates
#' 
#' @aliases AFScore
#'
#' @param MEs_table A data frame contains the eigengene values of each module, obtained in \code{\link{Mod_Extra}} 
#' @param comp_group The vector contains cell group/type class identity. The length of the vector must be equal to the number of rows of MEs_table
#' @param scale Logical. If TRUE, the score will be scaled according to group
#'
#' @details The order of the comp_group, and rownames of the MEs_table should correspond to each other
#' @returns a dataframe of AFscores for each cell group/type
#' @import stringr
#' @importFrom stats aggregate
#' @export

AFScore<-function(MEs_table,comp_group,scale = FALSE){
  if(is.null(MEs_table)){
    print("aggument 'MEs_table' is missing, with no default")
    break
  }
  if(is.null(comp_group)){
    print("aggument 'comp_group' is missing, with no default")
    break
  }
  
  MEs_table$group<-comp_group
  sd<-aggregate(MEs_table[,1:(ncol(MEs_table)-1)],list(MEs_table$group),sd)
  rownames(sd)<-sd[,1]
  sd<-sd[,-1]
  score<-as.data.frame(apply(sd,2,as.numeric),row.names = rownames(sd))
  if(scale){
    scale_score<-scale(t(score))
    scale_score<-data.frame(t(scale_score),row.names = rownames(score))
    return(scale_score)
  }else{
    return(score)
  }
}