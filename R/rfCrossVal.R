#' rfCrossVal
#'
#' This is the main function for the random forest cross-validation analysis
#' wherein top phenotype-associated regulators are used in random forest models
#' over a range of feature counts to predict the phenotype of interest. By
#' inspecting the output plot, one may gain a reasonable estimation for the
#' number of features (in this case regulators) required to minimize prediction
#' error.
#'
#' @param act_file string indicating the path and file name for the regulator
#' activities previously inferred by VIPER using the prepActExp function in
#' this package.
#' @param pheno_file string indicating the path and file name for the phenotype
#' data with the same sample in the same order as the regulator activities file.
#' @param phenotype string indicating the column name for the phenotype of
#' interest in the phenotype file.
#' @param prefix string indicating a prefix to be added to all output file
#' names.
#' @param seed numeric indicating the RNG seed. Defaults to 123
#' @param train_pct numeric indicating the proportion of samples to use in the
#' training set, while the remaining samples will go to the test set. Defaults
#' to 0.7.
#' @param min_regs numeric indicating the minimum number of regulators to use
#' in the random forest cross-validation testing for feature counts. This
#' number of regulators will be used in rfcv if the number of regulators that
#' associate with the phenotype at a Bonferroni-adjusted P less than the given
#' threshold is less than min_regs. Defaults to 100.
#' @param max_regs numeric indicating the maximum number of regulators to use
#' in the random forest cross-validation testing for feature counts. This
#' number of regulators will be used in rfcv if the number of regulators that
#' associate with the phenotype at a Bonferroni-adjusted P less than the given
#' threshold is greater than max_regs. Defaults to 500.
#' @param Bonf_thresh numeric indicating the significance threshold for
#' Bonferroni-adjusted P used to identify regulators best associated with the
#' phenotype of interest. Defaults to 0.05.
#' @param replicates numeric indicating the number of independent random forest
#' cross-validation analyses to run. Higher numbers will tend to generate
#' cleaner plots (up to a point), but will also add to the analysis time.
#' Defaults to 12.
#'
#' @return This will output a png image file of the random forest cross-
#' validation plot and a .RData workspace image file that can be used in the
#' final random forest modeling used to identify representative phenotypic
#' master regulators (MRs).
#' @export
#'
#' @examples
#' \dontrun{
#' rfCrossVal("Toy_data_regulators_activities.txt",
#' "GEUVADIS_simulated_phenotype.txt",
#' "Sim_pheno",
#' prefix="Toy_data")
#' }
rfCrossVal <- function(act_file,pheno_file,phenotype,prefix,seed=123,train_pct=0.7,min_regs=100,max_regs=500,Bonf_thresh=0.05,replicates=12){

  # Send standard output to a log file.
  sink(paste(prefix,"_rfcv_std_out.log",sep = ""),type = "output",split = T)

  # Read in data
  print("Reading in data")
  vip=utils::read.table(act_file,sep = "\t",header = T,row.names = 1)
  phenos=utils::read.table(pheno_file,sep = "\t",header = T,row.names = 1)

  # Test for binary vs continuous phenotype. If neither, stop analysis.
  num_pheno=is.numeric(phenos[,phenotype])
  bin_pheno=(length(levels(as.factor(phenos[,phenotype])))==2)
  if(!num_pheno && !bin_pheno){
    print("Phenotype data must be either continuous numeric or binary (numeric or categorical)")
    q("n")
  }

  # Split data for training/test sets, but make sure the phenotype distribution is comparable between sets.
  print("Splitting samples into training and test sets")
  suppressWarnings(RNGkind(sample.kind = "Rounding")) # This is now necessary after R v3.6.0 for consistent results with pre-3.6.0 scripts because the default sampler changed.
  set.seed(seed)
  rnd_samples=sample(colnames(vip),size = ceiling(dim(vip)[2]*train_pct)) # train_pct of samples rounded up to nearest integer
  train_vip=vip[,rnd_samples]
  train_phenos=phenos[rnd_samples,phenotype]
  test_vip=vip[,!(colnames(vip) %in% rnd_samples)]
  test_phenos=phenos[!(rownames(phenos) %in% rnd_samples),phenotype]

  # These plots can indicate concerning deviations in the distributions of phenotype between the training and test sets
  print("Plot phenotype distributions for training and test sets. If concerning difference are observed, consider either re-running this analysis with a different sampling seed number (as a 5th argument), or performing outlier sample removal prior to re-running the analysis.")
  trainTestPlot(train_phenos,test_phenos,phenotype,binary=bin_pheno,prefix)

  # Identify the regulators whose expression or activities are best associated with the phenotype
  print("Identifying the regulators whose activities are most significantly associated with the phenotype")
  cvrf_regs=phenoAssoc(train_vip,train_phenos,binary=bin_pheno,min_regs,max_regs,Bonf_thresh)
  train_vip_top=train_vip[cvrf_regs,]
  test_vip_top=test_vip[cvrf_regs,]

  # Z-scale values prior to RF
  train_zvip <- zscale.genes(train_vip_top)
  test_zvip <- zscale.genes(test_vip_top)

  # Run cross-validation random forest to identify the number of features required to plateau the error
  # using the above identified associated regulators
  print("Running random forest cross-validation analyses across different feature counts with 12 different sampling seeds")
  cv_rndFor=list()
  cl=parallel::makeCluster(future::availableCores())
  doParallel::registerDoParallel(cl)
  if(bin_pheno){ # If phenotype is binary...
    cv_rndFor=foreach::foreach(i=1:replicates,.packages = 'MRaQTL') %dopar%
      MRaQTL::rfcv.wrapper(seed=i, train_x=t(train_zvip), train_y=as.factor(train_phenos))
  } else{ # If phenotype is continuous...
    cv_rndFor=foreach::foreach(i=1:replicates,.packages = 'MRaQTL') %dopar%
      MRaQTL::rfcv.wrapper(seed=i, train_x=t(train_zvip), train_y=train_phenos)
  }
  parallel::stopCluster(cl)

  # Average across multiple different seeds for the predictors counts in the rfcv() analysis, and then plot it.
  rfcvPlot <- rfcvPlot(cv_rndFor,phenotype,prefix)

  # Save workspace for easy loading in the final MR analysis
  print("Saving workspace that will be used in the final master regulator (MR) analysis")
  save(list=ls(all.names=T),file=paste(prefix,"_rfcv_workspace.RData",sep = ""),envir=environment())

  # RFCV analysis finished
  print("Random forest cross-validation analysis is complete. Look at the plot of mean cross-validation error across feature counts to determine how many features are needed to minimize prediction error while keeping the feature count as low as possible. This feature count will be used in the following master regulator analysis to pick out that number of putative master regulators (MRs) for your phenotype of interest.")

  closeAllConnections()

}
