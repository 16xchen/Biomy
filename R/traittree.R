#' Clustering strains base on quantitative trait data
#'
#' Calculates p-values for hierarchical clustering via multiscale bootstrap resampling.
#'
#' Hierarchical clustering is done for given data and p-values are computed for each of the clusters.
#'
#'
#' @param microdata A matrix of relative trait data
#' @param nboot The number of bootstrap replications. The default is 1000.
#' @return An object of dendrgram class
#' @examples
#' set.seed(1234)  #simulate microbiome relative trait data
#' x <- rnorm(n=(16*12), mean=10, sd=10)
#' trait=matrix(x,16,12)
#' rname=paste("trait",1:16,sep="")
#' cname=paste("strain", 1:12, sep="")
#' rownames(trait)=rname
#' colnames(trait)=cname
#'
#' microdend=microtree(microdata=trait, nboot=10)
#' microdend
#' plot(microdend)
#'
#'@author Xingyao Chen
#' @export
#'
#'
#'
traittree=function(mircodata, nboot=1000){
  require("pvclust")  || install.packages("pvclust")
  library("pvclust")
  require("dendextend") || install.packages('dendextend')
  library(dendextend)
  a.clust = pvclust(trait, method.dist="cor", method.hclust="complete", nboot=nboot)
  return(as.dendrogram(a.clust))
}
