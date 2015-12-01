#' Subestting for highest correlations
#'
#' Subsets dataframe from {snpcor} function to contain only the highly correlated SNP dendrograms
#'
#'
#' @param cordata A dataframe generated from the function snpcor
#' @param threshold Percent thrshold (0<threshold<1) to search for highest correlations, default is 0.95
#' @return a subsetted correlation dataframe
#' @examples
#' set.seed(1234) #simulate microbiome relative abundance dataframe
#' x <- rnorm(n=(16*12), mean=10, sd=10)
#' trait=matrix(x,16,12)
#' rname=paste("trait",1:16,sep="")
#' cname=paste("strain", 1:12, sep="")
#' rownames(trait)=rname
#' colnames(trait)=cname
#'
#' mymicro.dend=microtree(microdata=abundance, nboot=10) #cluster strains by microbiome
#'
#' chrX = SampSNP()  #load SNP data
#' strain.names=colnames(abundance) #assign the same strain names
#' mysnptree=maketree(snpdata=chrX, strain.names=strain.names) #cluster strains by every 100 SNPs
#'
#' mycordata=snpcor(snpdata=chrX, snptree=mysnptree, micro.dend=mymicro.dend) #calculate correlation coefficient for micribiome dendrogram and each SNP dendrogram
#' mycordata.best=snpcor.best(cordata=mycordata,threshold = 0.95) #find highly correlated SNP dendrograms
#' @author Xingyao Chen
#' @export

snpcor.best =function(cordata, threshold){
  cor=cordata[,4]
  best <- which(cor>(threshold*(max(na.omit(cor)))))
  cordata.best=cordata[c(best),]
  return(as.data.frame(cordata.best))
}
