% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snpcor.best.R
\name{snpcor.best}
\alias{snpcor.best}
\title{Subestting for highest correlations}
\usage{
snpcor.best(cordata, threshold)
}
\arguments{
\item{cordata}{q dataframe generated from the function snpcor}

\item{threshold}{percent threshold (0-1) to search for highest correlations, default is 0.95}
}
\value{
a subsetted correlation dataframe
}
\description{
Subsets dataframe from {snpcor} function to contain only the highly correlated SNP dendrograms
}
\examples{
set.seed(1234) #simulate quantitative trait dataframe
x <- rnorm(n=(16*12), mean=10, sd=30)
trait=matrix(x,16,12)
rname=paste("trait",1:16,sep="")
cname=paste("strain", 1:12, sep="")
rownames(trait)=rname
colnames(trait)=cname

mytrait.dend=traittree(traitdata=trait, nboot=10) #cluster strains by traitbiome

chrX = SampSNP()  #load SNP data
strain.names=colnames(trait) #assign the same strain names
mysnptree=maketree(snpdata=chrX, strain.names=strain.names) #cluster strains by every 100 SNPs

mycordata=snpcor(snpdata=chrX, snptree=mysnptree, trait.dend=mytrait.dend) #calculate correlation coefficient for micribiome dendrogram and each SNP dendrogram
mycordata.best=snpcor.best(cordata=mycordata,threshold = 0.95) #find highly correlated SNP dendrograms
}
\author{
Xingyao Chen
}

