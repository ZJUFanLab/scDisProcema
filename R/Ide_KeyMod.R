#' @name Ide_KeyMod
#' @title Key modules identification
#'
#' @description This function is able to identify key modules for given cell groups according to the module significance score
#' @aliases Ide_KeyMod
#' 
#' @param MEs_table A data.frame contains the eigengene values of each module, obtained in \code{\link{Mod_Extra}}
#' @param MS_score A data.frame contains the module significance score from \code{\link{MSScore}}. If NULL, the score will be calculated
#' @param comp_group The vector contains cell group/type class identity. The length of the vector must be equal to the number of rows of MEs_table
#' @param values The trait values with which to compare, corresponding with comp_group. The value must be numerical and the length of which must be equal to the length of comp_group
#' @param Scc_scale Logical. If TRUE, Scaling is performed when calculating the correlation coefficient score
#' @param Saf_scale Logical. If TRUE, Scaling is performed when calculating the activity fluctuation score
#' @details This function will also export the absolute values of module significance score of corresponding key module
#' @returns a dataframe of key modules for each given group
#' @import stringr
#' @import WGCNA
#' @importFrom stats aggregate
#' @importFrom Seurat FindAllMarkers
#' @export
Ide_KeyMod<-function(MEs_table,MS_score = NULL,comp_group = NULL,values = NULL,Scc_scale = FALSE,Saf_scale = FALSE){
  if(is.null(MEs_table)){
    print("aggument 'MEs_table' is missing, with no default")
    break
  }
  if(!is.null(MS_score)){
    sms<-MS_score
  }else{
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
    sms<-MSScore(MEs_table,comp_group = comp_group,values = values,
                 Scc_scale = Scc_scale,Saf_scale = Saf_scale)
  }
  sms_abb<-abs(sms)
  module<-str_remove_all(colnames(sms_abb),"ME")
  sms_keymodule<-data.frame(Max = apply(sms_abb, 1, max),row.names = rownames(sms_abb))
  Keymodule<-c()
  for (m in 1:nrow(sms_abb)) {
    Keymodule[m]<-module[which(sms_abb[m,] == sms_keymodule$Max[m])]
  }
  sms_keymodule$Keymodule<-Keymodule
  
  return(sms_keymodule)
}
