#' rfcvPlot
#'
#' This function determines summary statistics for a list of rfcv objects
#' across all different feature counts, and from them, creates a png image of
#' the random forest cross-validation plot showing prediction error as a
#' function of feature counts.
#'
#' @param rfcv_list list of length 2 or more of rfcv outputs.
#' @param phenotype string indicating phenotype of interest, which is just used
#' in plot title.
#' @param prefix string indicating the output file prefix.
#'
#' @return Writes a png file of the random forest cross-validation plot that
#' can indicate how many regulators are required to minimize prediction error.
#' @export
#'
#' @examples
#' \dontrun{
#' rfcvPlot <- rfcvPlot(cv_rndFor,"Sim_pheno","Toy_data")
#' }
rfcvPlot <- function(rfcv_list,phenotype,prefix){

  # Average across however many different seeds for the feature counts in the rfcv() analysis, and then plot it.
  cv_rndFor_error=rfcv_list[[1]]$error.cv
  for(i in 2:length(rfcv_list)){
    cv_rndFor_error=cbind(cv_rndFor_error,rfcv_list[[i]]$error.cv)
  }
  cv_avg=rowMeans(cv_rndFor_error)
  cv_sd=MatrixGenerics::rowSds(cv_rndFor_error)

  grDevices::png(paste(prefix,"_RF_cross-validation_error_across_feature_counts.png",sep=""),width = 800,height = 600)
  graphics::par(mfrow=c(1,1))
  plot(x=as.numeric(names(cv_avg)),
       y=cv_avg,
       ylim=range(c(cv_avg-cv_sd,cv_avg+cv_sd)),
       pch=19,
       xlab="Feature count",
       ylab="Mean cross-validation error",
       main=paste("Average",phenotype,"rfcv() results over 12 different sampling seeds",sep = " "))
  graphics::arrows(x0=as.numeric(names(cv_avg)),
         y0=cv_avg-cv_sd,
         x1=as.numeric(names(cv_avg)),
         y1=cv_avg+cv_sd,
         length = 0.05,
         angle = 90,
         code = 3) # This draws the error bars
  grDevices::dev.off()

}
