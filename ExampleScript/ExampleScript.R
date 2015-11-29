#Example script using package Biomy
#Correlating relative abunadnce microbiome of host strain to Single Nucleotide Polymorphisms of host strain
#author: Xingyao Chen
#
install_github("16xchen/Biomy")
library(Biomy)
rm(list=ls())
#
chrX = read.csv("~/Biomy/chrX.csv")
mysnptree=maketree(chrX, strain.names=c(paste("strain", 1:(ncol(chrX)-1), sep="")))
#
set.seed(1234)
x <- rnorm(n=120, mean=5, sd=20)
abundance=matrix(x,10,12)
rname=paste("bacteria",1:10,sep="")
cname=paste("strain", 1:12, sep="")
rownames(abundance)=rname
colnames(abundance)=cname
#
mymicro.dend=microtree(abundance, nboot=10)
mycordata=snpcor(snpdata=chrX, snptree=mysnptree, micro.dend=mymicro.dend)
mycordata.best=snpcor.best(cordata=mycordata,threshold = 0.9)
drawtangle(snptree=mysnptree, micro.dend=mymicro.dend, cordata.best=mycordata.best, chr.num="X")
