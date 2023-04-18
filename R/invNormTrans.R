#' invNormTrans
#'
#' This function applies an inverse normal transformation (INT) on the given
#' gene metrics using the RankNorm() function from the RNOmni package, but is
#' more rapid due to the use of multi-threading when possible.
#'
#' @param genes matrix or data.frame containing gene metrics with genes in rows
#' and samples in columns. Gene and sample identifiers should be in rownames
#' and colnames, respectively.
#'
#' @return Returns a data.frame in the same format as the input gene data, but
#' with inverse normal transformed values.
#' @export
#'
#' @examples
#' \dontrun{
#' int_genes <- invNormTrans(genes)
#' }
invNormTrans <- function(genes){
  cl=parallel::makeCluster(future::availableCores())
  doParallel::registerDoParallel(cl)
  int_exp=list()
  int_exp=foreach::foreach(i=1:dim(genes)[1],.packages = 'RNOmni') %dopar%
    RNOmni::RankNorm(as.numeric(genes[i,]))
  int_exp=as.data.frame(t(structure(int_exp, row.names = c(NA, -length(int_exp[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
  parallel::stopCluster(cl)
  colnames(int_exp)=colnames(genes)
  rownames(int_exp)=rownames(genes)
  return(int_exp)
}
