#' zscale.genes
#'
#' This function applies Z-scaling for all genes across samples very rapidly
#' using multi-threading when possible.
#'
#' @param genes matrix or data.frame containing gene metrics with genes in rows
#' and samples in columns. Gene and sample identifiers should be in rownames
#' and colnames, respectively.
#'
#' @return Returns an expression set object (from DESeq2) with the gene metrics
#' Z-scaled across samples.
#' @export
#'
#' @examples
#' \dontrun{
#' z_genes <- zscale.genes(genes)
#' }
zscale.genes <- function(genes){
  cl=parallel::makeCluster(future::availableCores())
  doParallel::registerDoParallel(cl)
  genes_z=list()
  genes_z=foreach::foreach(i=1:dim(genes)[1]) %dopar%
    as.numeric(scale(as.numeric(genes[i,])))
  genes_z=as.data.frame(t(structure(genes_z, row.names = c(NA, -length(genes_z[[1]])), class = "data.frame"))) # Very rapidly convert a list to data.frame
  parallel::stopCluster(cl)
  colnames(genes_z)=colnames(genes)
  rownames(genes_z)=rownames(genes)
  z_set=Biobase::ExpressionSet(assayData=as.matrix(genes_z))
  return(z_set)
}
