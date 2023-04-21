#' grabMRs
#'
#' This function will grab the expression and activity data for anyn list of
#' genes, though it is implemented here for the purpose of grabbing inferred
#' phenotype MRs.
#'
#' @param exp_data string, data.frame or matrix. If a string, it must indicate
#' the path and file name for the expression data file. Otherwise, the given,
#' pre-existing data.frame or matrix containing the expression data will be
#' used.
#' @param act_data string, data.frame or matrix. If a string, it must indicate
#' the path and file name for the activity data file. Otherwise, the given,
#' pre-existing data.frame or matrix containing the activity data will be
#' used.
#' @param mr_list string or vector. If a string, it must indicate
#' the path and file name for the gene (MR) name list. Otherwise, the given,
#' pre-existing vector of gene (MR) names will be used.
#' @param prefix string indicating a prefix to be added to all output file
#' names.
#'
#' @return Writes to text file the expression and activity data for the given
#' list of genes/MRs.
#' @export
#'
#' @examples
#' \dontrun{
#' # Run with objects already in the workspace.
#' grabMRs(gene_exp,gene_act,mrs,"Toy_data")
#'
#' # Run with names of files to be read.
#' grabMRs("./Toy_data_expressed_genes_QNorm_INT_expression.txt",
#' "Toy_data_regulators_activities.txt",
#' "Toy_data_representative_Sim_pheno_MRs_list.txt",
#' "Toy_data")
#' }
grabMRs <- function(exp_data,act_data,mr_list,prefix){

  # Read in data
  if(length(exp_data)==1){
    gene_exp=utils::read.table(exp_data,sep = "\t",header = T,row.names = 1)
  } else {gene_exp=exp_data}

  if(length(act_data)==1){
    gene_act=utils::read.table(act_data,sep = "\t",header = T,row.names = 1)
  } else {gene_act=act_data}

  if(length(mr_list)==1){
    mrs=scan(mr_list,character())
  } else {mrs=mr_list}

  # Grab the MRs
  mrs_exp=gene_exp[mrs,]
  mrs_act=gene_act[mrs,]

  # Write data to file
  utils::write.table(cbind("Gene_name"=rownames(mrs_exp),mrs_exp),
                     paste(prefix,"_MRs_only_QNorm_INT_expression.txt",sep=""),
                     sep="\t",quote = F,row.names = F)
  utils::write.table(cbind("Gene_name"=rownames(mrs_act),mrs_act),
                     paste(prefix,"_MRs_only_activities.txt",sep=""),
                     sep="\t",quote = F,row.names = F)

}
