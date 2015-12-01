#' Generating dendrograms for the entire chromosome
#'
#' Uses the {finddistnace} function to cluster strains based on dissimilarity of SNPs. One drendrogram for every 100 snps
#'
#' @param snpdata SNP data downloaded from opensource databases such as <http://msub.csbio.unc.edu/>or <http://www.informatics.jax.org/javawi2/servlet/WIFetch?page=snpQF>
#' @param strain.names a character vector of strain names
#' @return a large list of dendrograms
#' @examples
#' chrX = SampSNP() #load SNP data
#' strain.names=c(paste("strain", 1:(ncol(chrX)-1), sep="")) #assign strain names, strain names must match between trait and SNP dendrograms
#' mysnptree=maketree(snpdata=chrX, strain.names=strain.names)
#' mysnptree
#' plot(mysnptree[1])
#'
#' @author Xingyao Chen
#' @export

maketree=function(snpdata, strain.names){
  snptree <- vector("list", (nrow(snpdata)-100)/50)
  for(i in  1:((nrow(snpdata)-100)/50)){
  snptree[[i]] <- finddistance(snpdata, start=i*50, strain.names)}
  return(snptree)
}
