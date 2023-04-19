#' phenoAssoc
#'
#' This function determines the association between the given phenotype and all
#' input training set activities and then returns a list of the top
#' Bonferroni-significant regulators (within a given min/max).
#'
#' @param activities data.frame containing the training set activities
#' @param phenotype data.frame containing the training set phenotype data
#' @param binary logical indicating whether the phenotype is binary. If FALSE,
#' the phenotype is assumed to be continuous. Categorical data with two levels
#' will be treated as binary, but more than two levels will not work. Defaults
#' to FALSE.
#' @param min_regs numeric indicating the minimum number of best associated
#' regulators to use in the cross-validation random forest analysis. Defaults
#' to 100.
#' @param max_regs numeric indicating the maximum number of best associated
#' regulators to use in the cross-validation random forest analysis. Defaults
#' to 500.
#' @param Bonf_thresh numeric indicating the Bonferroni threshold for
#' significance that will determine how many of the best associated regulators
#' to use in the cross-validation random forest analysis as long as the number
#' lies within the given min/max. Defaults to 0.05.
#'
#' @return Returns a vector of names for the best associated regulators that
#' will be used in the cross-validation random forest analysis
#' @export
#'
#' @examples
#' \dontrun{
#' phenoAssoc(train_act,train_pheno,prefix="Toy_data")
#' }
phenoAssoc <- function(activities,phenotype,binary=FALSE,min_regs=100,max_regs=500,Bonf_thresh=0.05){

  if(binary){ # If phenotype is binary...
    a_pheno=list()
    for(i in 1:dim(activities)[1]){
      a_pheno[[i]]=summary(stats::glm((as.numeric(as.factor(phenotype))-1)~as.numeric(activities[i,]),family = stats::binomial())) # Have to coerce the phenotype data in this way to ensure it is interpreted as a binary input
    }
  } else{ # If phenotype is continuous...
    a_pheno=list()
    for(i in 1:dim(activities)[1]){
      a_pheno[[i]]=summary(stats::lm(phenotype~as.numeric(activities[i,])))
    }
  }
  names(a_pheno)=rownames(activities)
  sum_a_pheno=as.data.frame(matrix(nrow = length(a_pheno),ncol = 2))
  for(i in 1:length(a_pheno)){
    sum_a_pheno[i,1]=stats::coef(a_pheno[[i]])[2,1]
    sum_a_pheno[i,2]=stats::coef(a_pheno[[i]])[2,4]
  }
  colnames(sum_a_pheno)=c("Beta","P")
  rownames(sum_a_pheno)=names(a_pheno)
  sum_a_pheno=sum_a_pheno[order(sum_a_pheno$P),] # re-order regulators by significance of association w/ phenotype
  sum_a_pheno$BonferroniP=sum_a_pheno$P*dim(sum_a_pheno)[1]
  bonf_sig=sum(sum_a_pheno$BonferroniP<Bonf_thresh)
  print(paste(bonf_sig,"regulators associated with phenotype with Bonferroni-adjusted P<",as.character(Bonf_thresh),sep = " "))

  # Select regulators most significantly associated with phenotype for cross-validation random forest within the bounds of [min_regs,max_regs]
  if(bonf_sig>min_regs){
    if(bonf_sig<max_regs){
    cvrf_regs=rownames(sum_a_pheno)[1:bonf_sig]
    } else {
    cvrf_regs=rownames(sum_a_pheno)[1:max_regs]
    }
  } else{
    cvrf_regs=rownames(sum_a_pheno)[1:min_regs]
  }

  return(cvrf_regs)

}
