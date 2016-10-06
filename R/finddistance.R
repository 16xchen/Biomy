#' Clustering strains based on dissimilarity of SNPs
#'
#' One drendrogram for every 100 snps
#'
#' @param snpdata SNP data downloaded from opensource databases such as <http://msub.csbio.unc.edu/>or <http://www.informatics.jax.org/javawi2/servlet/WIFetch?page=snpQF>
#' @param start starting row of SNP, default is 1
#' @param strain.names a character vector of matching strain names
#' @return an object of class dendrogram
#' @author Xingyao Chen
#' @export
#'
finddistance <- function(snpdata, start=1, strain.names){
  if(!require("dendextend"))
    install.packages('dendextend',  repos = 'http://cran.us.r-project.org')
  library(dendextend)
  dat <- snpdata[start:(start+100),]
  colnames(dat)<- c("position", strain.names)
  N <- ncol(dat)-1
  dist <- matrix(0, N, N)
  rownames(dist) <- colnames(dist) <- names(dat)[-1]
  for (i in 1:N)
    for(j in 1:N)
      dist[i,j] <- length(which((dat[,i+1]!=dat[,j+1])==TRUE))
  cl <- hclust(as.dist(dist))
  return(as.dendrogram(cl))
}
