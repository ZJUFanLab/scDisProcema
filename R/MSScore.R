#' @name MSScore
#' @title Calculate the module significance score
#'
#' @description  For each cell group/type, this function can assess the comprehensive correlation of inferred gene modules by scoring them
#' 
#' @aliases MSScore
#' 
#' 
#' 
#' @param MEs_table A data.frame contains the eigengene values of each module, obtained in \code{\link{Mod_Extra}}.
#' @param comp_group The vector contains cell group/type class identity. The length of the vector must be equal to the number of rows of MEs_table
#' @param values The trait values with which to compare, corresponding with comp_group. The value must be numerical and the length of which must be equal to the length of comp_group
#' @param Scc_scale Logical. If TRUE, Scaling is performed when calculating the correlation coefficient score
#' @param Saf_scale Logical. If TRUE, Scaling is performed when calculating the activity fluctuation score
#' @details The correlation coefficient score and activity fluctuation score will be calculated together with module significance score with MSScore function. If you want to obtain them respectively, run \code{\link{CCScore}} and \code{\link{AFScore}} seperately
#'          The order of the comp_group, values and rownames of the MEs_table should correspond to each other
#' @returns a dataframe of MSscores for each cell group/type
#' @import stringr
#' @importFrom stats aggregate
#' @importFrom Hmisc rcorr
#' @export

MSScore<-function(MEs_table,comp_group,values,Scc_scale = FALSE,Saf_scale = FALSE){
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
  scc<-CCScore(MEs_table,comp_group = comp_group,values = values,scale = Scc_scale)
  saf<-AFScore(MEs_table,comp_group = comp_group,scale = Saf_scale)
  score<-scc*saf
  
  scale_score<-scale(t(score))
  scale_score<-data.frame(t(scale_score),row.names = rownames(score))
  return(scale_score)
}
