#' Making tanglegrams of microbiome tree matched with snp trees
#'
#' @param snptree A list of dendrograms generated from snp data
#' @param micro.dend A dendrogram generated from microbiome relative abundance data
#' @param cordata A dataframe which includes positions and corrlation values of highly correlated snp trees
#' @param chr.num The chromosome number from which the snps are located
#' @return Plots highlt correlated tanglegrams with SNP location on chromosome in basepairs
#' @examples
#' set.seed(1234) #simulate microbiome relative abundance dataframe
#' x <- rnorm(n=120, mean=5, sd=20)
#' abundance=matrix(x,10,12)
#' rname=paste("bacteria",1:10,sep="")
#' cname=paste("species", 1:12, sep="")
#' rownames(abundance)=rname
#' colnames(abundance)=cname
#'
#' mymicro.dend=microtree(microdata=abundance, nboot=10) #cluster strains by microbiome
#'
#' chrX = read.csv("~/Biomy/chrX.csv")  #load SNP data from opensource database
#' strain.names=colnames(abundance) #assign the same strain names
#' mysnptree=maketree(snpdata=chrX, strain.names=strain.names) #cluster strains by every 100 SNPs
#'
#' mycordata=snpcor(snpdata=chrX, snptree=mysnptree, micro.dend=mymicro.dend) #calculate correlation coefficient for micribiome dendrogram and each SNP dendrogram
#' mycordata.best=snpcor.best(cordata=mycordata,threshold = 0.95) #find highly correlated SNP dendrograms
#' drawtangle(snptree=mysnptree, micro.dend=mymicro.dend, cordata.best=mycordata.best, chr.num="X") #plot out as tanglegrams
#'
#' @author Xingyao Chen
#' @export
#'
#'
drawtangle = function(snptree, micro.dend, cordata.best, chr.num)
{
  require("dendextend") || install.packages('dendextend')
  library(dendextend)
  if(ncol(cordata.best)==1){
    cordata.best=t(cordata.best)
  }
  for(j in 1:nrow(cordata.best)){
  tanglegram(untangle(intersect_trees(micro.dend, snptree[[cordata.best[j,1]]])),
             main="Correlations",cex_main=1,
             main_left = "Microbiome tree",
             main_right= paste("chr",chr.num, "SNP @", cordata.best[j,2],"to", cordata.best[j,3], "bp", sep=" "),
             sub=paste("correlation",round(cordata.best[j, 4],digit=4),
                       sep=":"), cex_sub=0.9, highlight_distinct_edges=F,
             lwd=1.5, columns_width=c(4.5,1.8,4.5),
             margin_outer = 3,
             margin_inner=6)
}
}
