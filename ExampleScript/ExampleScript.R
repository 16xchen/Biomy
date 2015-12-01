#Example script using package Biomy
#Correlating quantitative traits of strains to Single Nucleotide Polymorphisms of strains
#author: Xingyao Chen
#date: 11/29/15
#
install.packages("devtools")
library(devtools)
install_github("16xchen/Biomy")
library(Biomy)
rm(list=ls())
#
#load sample SNP data from UNC CSBio
chrX = SampSNP()
mysnptree=maketree(chrX, strain.names=c(paste("strain", 1:(ncol(chrX)-1), sep="")))
#
#simulate trait dataset
set.seed(1234)
x <- rnorm(n=(16*12), mean=10, sd=30)
trait=matrix(x,16,12)
rname=paste("trait",1:16,sep="")
cname=paste("strain", 1:12, sep="")
rownames(trait)=rname
colnames(trait)=cname
#
#run the correlation analysis
mytrait.dend=traittree(trait, nboot=10)
mycordata=snpcor(snpdata=chrX, snptree=mysnptree, trait.dend=mytrait.dend)
mycordata.best=snpcor.best(cordata=mycordata,threshold = 0.98)
drawtangle(snptree=mysnptree, trait.dend=mytrait.dend, indata=mycordata[44:46,], chr.num="X")
