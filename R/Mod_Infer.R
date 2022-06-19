#' @name Mod_Infer
#' @title WGCNA module inference
#' 
#' @description According to the input expression matrix, the WGCNA gene module was obtained
#' @aliases Mod_Infer
#'
#' @param data_Expr A data.frame generated from \code{\link{Mean_Expre}}
#' @param sftpower The soft-threshold to identify scale-free gene network
#' @param RsquaredCut The desired minimum scale free topology fitting index. Default is 0.9
#' @param maxBlockSize Maximum cluster gap given as the fraction.Default is all of genes of given dataset
#' @param minModuleSize Min cluster gap given as the fraction.Default is 20
#' @param deepSplit Parameter that control the resolution sensitivity of the module
#' @param mergeCutHeight Parameter for module merging.The larger this parameter is, the fewer modules there are
#' @details When sftpower is NULL, the appropriate soft threshold is inferred through scale-free network analysis. And the inferred or input sftpower would be in the name of 'sftpower' of the exported list
#'          
#' @return A list of inferred network information. Details can be seen in \code{\link[WGCNA]{blockwiseModules}}
#' @import WGCNA 
#' @export
Mod_Infer<-function(data_Expr,sftpower = NULL,RsquaredCut = 0.9,maxBlockSize = "total",
                    deepSplit = 4,minModuleSize = 20,
                    mergeCutHeight = 0.1)
  {
  if(is.null(data_Expr)){
    print("aggument 'data_Expr' is missing, with no default")
    break
  }
  if(is.null(sftpower)){
    sft <- pickSoftThreshold(data_Expr, RsquaredCut = RsquaredCut,verbose = 5)
    sftpower<-sft$powerEstimate
  }
  if(maxBlockSize == "total"){
    maxBlockSize<-ncol(data_Expr)
  }
  net <- blockwiseModules(
    data_Expr,
    power = sftpower,
    maxBlockSize = maxBlockSize,
    TOMType = "unsigned",
    deepSplit = deepSplit, minModuleSize = minModuleSize,
    mergeCutHeight = mergeCutHeight,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = TRUE,
    saveTOMFileBase = "FPKM-TOM",
    loadTOMs = TRUE,
    verbose = 3
  )
  net$sftpower<-sftpower
  return(net)
 }


