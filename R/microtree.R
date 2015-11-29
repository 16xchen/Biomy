#' Clustering strains base on microbiome relative abundance data
#'
#' Calculates p-values for hierarchical clustering via multiscale bootstrap resampling.
#'
#' Hierarchical clustering is done for given data and p-values are computed for each of the clusters.
#'
#'
#' @param microdata A matrix of relative abundance data
#' @param nboot The number of bootstrap replications. The default is 1000.
#' @return An object of dendrgram class
#' @examples
#' set.seed(1234)  #simulate microbiome relative abundance data
#' x <- rnorm(n=120, mean=5, sd=20)
#' abundance=matrix(x,10,12)
#' rname=paste("bacteria",1:10,sep="")
#' cname=paste("species", 1:12, sep="")
#' rownames(abundance)=rname
#' colnames(abundance)=cname
#' #
#' microdend=microtree(microdata=abundance, nboot=10)
#' microdend
#' plot(microdend)
#'
#'@author Xingyao Chen
#' @export
#'
#'
#'
microtree=function(mircodata, nboot=1000){
  require("pvclust")  || install.packages("pvclust")
  library("pvclust")
  require("dendextend") || install.packages('dendextend')
  library(dendextend)
  a.clust = pvclust(abundance, method.dist="cor", method.hclust="complete", nboot=nboot)
  return(as.dendrogram(a.clust))
}
