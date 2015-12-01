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
chrX = SampSNP()
mysnptree=maketree(chrX, strain.names=c(paste("strain", 1:(ncol(chrX)-1), sep="")))
#
set.seed(1234)
x <- rnorm(n=(16*12), mean=10, sd=10)
trait=matrix(x,16,12)
rname=paste("trait",1:16,sep="")
cname=paste("strain", 1:12, sep="")
rownames(trait)=rname
colnames(trait)=cname
#
mymicro.dend=microtree(trait, nboot=10)
mycordata=snpcor(snpdata=chrX, snptree=mysnptree, micro.dend=mymicro.dend)
mycordata.best=snpcor.best(cordata=mycordata,threshold = 0.9)
drawtangle(snptree=mysnptree, micro.dend=mymicro.dend, cordata.best=mycordata.best, chr.num="X")
