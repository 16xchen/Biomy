#' Matching dendrograms
#'
#' Making tanglegrams of trait dendrograms matched up with snp dendrograms
#'
#' @param snptree A list of dendrograms generated from snp data
#' @param trait.dend A dendrogram generated from quantitative trait data
#' @param cordata A dataframe which includes positions and corrlation values of highly correlated snp trees
#' @param chr.num The chromosome number from which the snps are located
#' @return Plots highlt correlated tanglegrams with SNP location on chromosome in basepairs
#' @examples
#' set.seed(1234) #simulate quantitatie trait dataframe
#' x <- rnorm(n=(16*12), mean=10, sd=10)
#' trait=matrix(x,16,12)
#' rname=paste("trait",1:16,sep="")
#' cname=paste("strain", 1:12, sep="")
#' rownames(trait)=rname
#' colnames(trait)=cname
#'
#' mytrait.dend=traittree(traitdata=trait, nboot=10) #cluster strains by trait
#'
#' chrX = SampSNP()  #load SNP data
#' strain.names=colnames(trait) #assign the same strain names
#' mysnptree=maketree(snpdata=chrX, strain.names=strain.names) #cluster strains by every 100 SNPs
#'
#' mycordata=snpcor(snpdata=chrX, snptree=mysnptree, trait.dend=mytrait.dend) #calculate Baker's correlation coefficient for trait dendrogram and each SNP dendrogram
#' mycordata.best=snpcor.best(cordata=mycordata,threshold = 0.95) #find highly correlated SNP dendrograms
#' drawtangle(snptree=mysnptree, trait.dend=mytrait.dend, cordata.best=mycordata.best, chr.num="X") #plot out as tanglegrams
#'
#' @author Xingyao Chen
#' @export
#'
#'
drawtangle = function(snptree, trait.dend, cordata.best, chr.num)
{
  require("dendextend") || install.packages('dendextend')
  library(dendextend)
  if(ncol(cordata.best)==1){
    cordata.best=t(cordata.best)
  }
  for(j in 1:nrow(cordata.best)){
  tanglegram(untangle(intersect_trees(trait.dend, snptree[[cordata.best[j,1]]])),
             main="Correlations",cex_main=1,
             main_left = "traitbiome tree",
             main_right= paste("chr",chr.num, "SNP @", cordata.best[j,2],"to", cordata.best[j,3], "bp", sep=" "),
             sub=paste("correlation",round(cordata.best[j, 4],digit=4),
                       sep=":"), cex_sub=0.9, highlight_distinct_edges=F,
             lwd=1.5, columns_width=c(4.5,1.8,4.5),
             margin_outer = 3,
             margin_inner=6)
}
}
