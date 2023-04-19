#' trainTestPlot
#'
#' This function generates a plot that compares the phenotype distributions of
#' the training and test sets to confirm they are comparable. A barplot is used
#' for binary phenotypes while a density plot is used for continuous phenotypes.
#'
#' @param train_pheno data.frame containing the training set phenotype data.
#' @param test_pheno data.frame containing the test set phenotype data.
#' @param phenotype string indicating the column name for the phenotype of
#' interest.
#' @param binary logical indicating whether the phenotype is binary. If FALSE,
#' the phenotype is assumed to be continuous. Categorical data with two levels
#' will be treated as binary, but more than two levels will not work. Defaults
#' to FALSE.
#' @param prefix string indicating the output file prefix that will be added to
#' the file name along with the column name for the phenotype of interest.
#'
#' @return Writes the phenotype distributions plot to png file.
#' @export
#'
#' @examples
#' \dontrun{
#' trainTestPlot(train_phenos,test_phenos,"Sim_pheno",prefix="Toy_data")
#' }
trainTestPlot <- function(train_pheno,test_pheno,phenotype,binary=FALSE,prefix){

  if(binary){ # If phenotype is binary...
    grDevices::png(paste(prefix,phenotype,"distributions_in_training_and_test_sets.png",sep="_"),width = 800,height = 600)
    graphics::par(mfrow=c(1,1))
    graphics::barplot(rbind(proportions(table(train_pheno)),proportions(table(test_pheno))),
            col = c("blue","red"),beside = T,legend = c("Train","Test"),
            xlab=phenotype,ylab="Proportion",main="Phenotype distributions")
    grDevices::dev.off()
  } else{ # If phenotype is continuous...
    grDevices::png(paste(prefix,phenotype,"distributions_in_training_and_test_sets.png",sep="_"),width = 800,height = 600)
    graphics::par(mfrow=c(1,1))
    plot(stats::density(train_pheno),xlab=phenotype,main="Phenotype distributions",col="blue")
    graphics::lines(density(test_pheno),col="red")
    graphics::legend(x="topright",legend=c("Training","Test"),lty=par("lty"),col=c("blue","red"),text.col=c("blue","red"))
    grDevices::dev.off()
  }

}
