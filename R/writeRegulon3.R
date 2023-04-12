#' write.regulon3
#'
#' Write regulon object to text file. This is a modification of the write.regulon() function included in the VIPER that is much faster due to multi-threading.
#'
#' @param regulon A regulon object generated from the aracne2regulon() function in the VIPER package.
#' @param file string for the file name to which the regulon should be written.
#' @param sep string for the character(s) or regular expression escaped character to be used as the field delimiter in the text file. Default is "\t", which is the tab-delimiter.
#' @param header logical indicating whether the header should be printed. Default is TRUE.
#' @param n numeric indicating the number of interactions to print. Default is Inf.
#' @param regulator string specifying a particular regulator. Default is NULL.
#' @param cpus numeric indicating the number of CPUs to use in parallel processing. Defaults to all detected available cores.
#' @param toScreen logical indicating whether the output should also be printed in the console. Defaults to FALSE.
#'
#' @return Write the regulon object to the indicated file and returns a data.frame version of the regulon formatted as in the output file.
#' @export
#'
#' @examples
#' \dontrun{
#' regulon<-aracne2regulon(tissue_net,genes_set,format="3col") # NOTE: The aracne network file cannot have a header!!!
#' regTable<-write.regulon3(regulon,file="tissue_interactome.txt",sep=""))
#' }
write.regulon3 <- function(regulon,file="",sep="\t",header=TRUE,n=Inf,regulator=NULL,cpus=future::availableCores(),toScreen=F) {
  fullTab<-data.frame()
  if(is.null(regulator)){
    cl=makeCluster(cpus)
    on.exit(stopCluster(cl))
    registerDoParallel(cl)
    tab=list()
    tab=foreach(tf=names(regulon)) %dopar%
      cbind(rep(tf,length(names(regulon[[tf]]$tfmode))),names(regulon[[tf]]$tfmode),regulon[[tf]]$tfmode,regulon[[tf]]$likelihood)
    fullTab=as.data.frame(do.call("rbind",tab)) # Very rapidly convert a list of data.frames to data.frame
  } else {
    tf<-regulator
    x<-regulon[[tf]]
    targets<-names(x$tfmode)
    moas<-x$tfmode
    likelihoods<-x$likelihood
    tab<-cbind(rep(tf,length(targets)),targets,moas,likelihoods)
    fullTab<-rbind(fullTab,tab)
  }
  colnames(fullTab)=c("Regulator","Target","MoA","likelihood")
  rownames(fullTab)=NULL

  if(file!=""){
    write.table(fullTab,file = file,sep = sep, col.names = header, row.names = F,quote = F)
  }

  if(toScreen){
    print(fullTab)
  }

  return(fullTab)
}
