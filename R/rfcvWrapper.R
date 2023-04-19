#' rfcv.wrapper
#'
#' This is a wrapper function for calling the rfcv function in the randomForest
#' package with a given seed. This used for multi-threading the rfcv analysis.
#'
#' @param train_x data.frame containing the activities for the top phenotype
#' associated regulators for the training set.
#' @param train_y data.frame containing the phenotype data for the training set.
#' @param seed numeric indicating the RNG seed.
#'
#' @return Returns a list with the following components:
#' list(n.var=n.var,error.cv=error.cv,predicted=cv.pred)
#' @export
#'
#' @examples
#' \dontrun{
#' cv_rndFor <- rfcv.wrapper(train_zact,train_pheno,seed=1)
#' }
rfcv.wrapper <- function(train_x,train_y,seed){

  set.seed(seed)
  temp=randomForest::rfcv(trainx=train_x, trainy=train_y, scale = F, step = -5, cv.fold = 5)
  return(temp)

}
