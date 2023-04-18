#' corDensPlot
#'
#' This function generates a png image file of a density plot of regulators'
#' activity-to-expression Pearson correlations.
#'
#' @param expression matrix or data.frame of expression values
#' @param activity matrix or data.frame of inferred activities
#' @param file_prefix string for the file name prefix for the output png file
#'
#' @return Writes the density plot for all regulators' matched activities and
#' expression values Pearson correlations to a png image file.
#' @export
#'
#' @examples
#' \dontrun{
#' corDensPlot(reg_expression,reg_activity,"Regulators")
#' # Generates the file "Regulators_activity_expression_correlation.png" in
#' # your working directory.
#' }
corDensPlot <- function(expression,activity,file_prefix){
  cl=parallel::makeCluster(future::availableCores())
  doParallel::registerDoParallel(cl)
  expActCor=list()
  startTime=Sys.time()
  expActCor=foreach::foreach(i=1:dim(expression)[1]) %dopar%
    stats::cor(t(expression)[,i],t(activity)[,i])
  expActCor=unlist(expActCor)
  stopCluster(cl)
  Sys.time()-startTime
  names(expActCor)=rownames(expression)

  grDevices::png(paste(file_prefix,"_activity_expression_correlation.png",sep=""),width = 800,height = 600)
  plot(stats::density(expActCor),main="Activity-Expression Correlation",xlab="Pearson r")
  grDevices::dev.off()
}
