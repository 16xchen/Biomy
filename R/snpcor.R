#' Correlating SNP dendrograms with the microbiome dendrogram
#'
#' Calculate Baker's Gamma correlation coefficient between two trees
#'
#' @param snpdata snpdata SNP data downloaded from opensource databases such as <http://msub.csbio.unc.edu/>or <http://www.informatics.jax.org/javawi2/servlet/WIFetch?page=snpQF>
#' @param snptree a list of dendrograms generated from snp data
#' @param micro.dend a dendrogram generated from microbiome relative abunance data
#' @return a data frame of positions, snp chromosomal location, and correlation coefficient
#' @examples
#' set.seed(1234) #simulate quantitative trait dataframe
#' x <- rnorm(n=(16*12), mean=10, sd=10)
#' trait=matrix(x,16,12)
#' rname=paste("trait",1:16,sep="")
#' cname=paste("strain", 1:12, sep="")
#' rownames(trait)=rname
#' colnames(trait)=cname
#'
#' mymicro.dend=microtree(microdata=trait, nboot=10) #cluster strains by microbiome
#'
#' chrX = SampSNP() #load SNP data
#' strain.names=colnames(trait) #assign the same strain names
#' mysnptree=maketree(snpdata=chrX, strain.names=strain.names) #cluster strains by every 100 SNPs
#'
#' mycordata=snpcor(snpdata=chrX, snptree=mysnptree, micro.dend=mymicro.dend) #calculate correlation coefficient for micribiome dendrogram and each SNP dendrogram
#'
#' @author Xingyao Chen
#' @export

snpcor=function(snpdata, snptree, micro.dend){
  require("dendextend") || install.packages('dendextend')
  require("corrplot") || install.packages("corrplot")
  library("dendextend")
  library("corrplot")

values <- numeric(length=(nrow(snpdata)-100)/50)
cor <- vector("list", (nrow(snpdata)-100)/50)
for(i in  1:((nrow(snpdata)-100)/50)){
  dend_both <- dendlist(intersect_trees(micro.dend, snptree[[i]]))
  cor[[i]] <- cor.dendlist(dend_both)
  values[i] <- cor[[i]][1,2]
}
names(values) <- c(seq(1,nrow(snpdata)-150, by=50))
cor <- as.numeric(values)
pos.start <- as.numeric(names(values))
pos.end <- as.numeric(names(values))+100
cordata <- cbind(1:length(values), snpdata[,1][c(pos.start)], snpdata[,1][c(pos.end)], cor)
colnames(cordata)=c("pos","start", "end", "cor")
return(cordata)
}




