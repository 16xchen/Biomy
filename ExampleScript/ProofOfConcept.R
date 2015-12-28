#Example script using package Biomy
#Correlating calorie intake of strains to Single Nucleotide Polymorphisms of strains
#author: Xingyao Chen
#date: 11/29/17
#
install.packages("devtools")
library(devtools)
install_github("16xchen/Biomy")
library(Biomy)
library(reshape2)
rm(list=ls())
#
#load sample SNP data from UNC CSBio
chr4=SampSNP()
#load sample Phenotype data from Mouse Pheome Dataset
pheno=SampPheno()
pheno.dat=pheno[,c(1,2,4)]
mypheno=dcast(pheno.dat, strain~varname)
rownames(mypheno)=mypheno[,1]
mypheno=mypheno[,-1]
mypheno=t(mypheno)
#
#define names, names must match between trees
new.names <- c("129S1/J","A/J","AKR/J","BALB/cJ","C3H/HeJ","C77BL/6J","CBA/J","DBA/2J","FVB/NJ","NOD/ShiLtJ","NU/J","SJL/J")
#
#run the correlation analysis
mysnptree=maketree(chr4, strain.names=new.names)
mytrait.dend=traittree(mypheno, nboot=100)
mycordata=snpcor(snpdata=chr4, snptree=mysnptree, trait.dend=mytrait.dend)
mycordata.best=snpcor.best(cordata=mycordata,threshold = 0.8)
drawtangle(snptree=mysnptree, trait.dend=mytrait.dend, indata=mycordata.best[1,], chr.num=4)
for(i in 1:nrow(mycordata.best)){
  drawtangle(snptree=mysnptree, trait.dend=mytrait.dend, indata=mycordata.best[i,], chr.num=4)
}





