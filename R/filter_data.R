#' @name filter_data
#' @title Filter data
#' 
#' @description filter the data with Mean Absolute Deviation and goodsample function
#' 
#' @param data_Expr A data.frame generated from \code{\link{Mean_Expre}}
#' @param mad Logical. Whether to use the MAD values to filter data. Default is TRUE
#' @param mad.thre The threshold to filter the mad.Default is max(max(m.mad)*0.25,0.01)
#' @param goodsamples Logical. Whether to use the goodSamplesGenes function to filter data. Default is TRUE

#' @details MAD values are used to screen for features with variability, and the goodSamplesGenes function is used to check for missing values
#' @return a data.frame with filtered values
#' @import stringr
#' @importFrom WGCNA goodSamplesGenes
#' @importFrom dynamicTreeCut printFlush
#' @importFrom stats mad
#' @importFrom stats quantile
#' @export
filter_data<-function(data_Expr,mad = TRUE,mad.thre = NULL,goodsamples = TRUE){
  if(is.null(data_Expr)){
    print("aggument 'data_Expr' is missing, with no default")
    break
  }
  if(mad == T){
    m.mad <- apply(data_Expr,2,mad)
    if(!is.null(mad.thre)){
      data_Expr <- data_Expr[,which(m.mad > mad.thre)]
    }else{
      data_Expr <- data_Expr[,which(m.mad > max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01))]
    }
  }
  if(goodsamples == T){
    gsg = WGCNA::goodSamplesGenes(data_Expr, verbose = 3)
    gsg$allOK
    gsg$goodSamples
    
    if (!gsg$allOK){
      # Optionally, print the gene and sample names that were removed:
      if (sum(!gsg$goodGenes)>0) 
        dynamicTreeCut::printFlush(paste("Removing genes:", 
                         paste(names(data_Expr)[!gsg$goodGenes], collapse = ",")));
      if (sum(!gsg$goodSamples)>0) 
        dynamicTreeCut::printFlush(paste("Removing samples:", 
                         paste(rownames(data_Expr)[!gsg$goodSamples], collapse = ",")));
      # Remove the offending genes and samples from the data:
      data_Expr = data_Expr[gsg$goodSamples, gsg$goodGenes]
    }
  }
  return(data_Expr)
}
