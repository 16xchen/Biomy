## Biomy
##Description: 
  Biomy, an add-on package for the widely used statistical language R, is a novel tool for associating quantitative traits to genes. 
	By comparing the dissimilarity among strains based on both phenotypic and genetic data, Biomy links a set of quantitative traits to specific regions in the genome. If the strains vary phenotypically the same way they vary genetically, then Biomy assumes there is a high probability that the genome region interacts with the phenotype of interest. 
	The function {traittree} partitions a set of mouse strains by their quantitative phenotypes using clustering analysis. {traittree} implements {pvclust}, which is an R package that assesses the uncertainty in hierarchical clustering with p-values. A p-value of a cluster indicates how strongly the cluster is supported by data. For each cluster, p-values are calculated via multiscale bootstrap resampling. {trainttree} returns an object of the dendrogram class.
	The function {maketree} partitions the same set of strains based on their SNP data. Starting from position 1 on a chromosome, {maketree} generates a dissimilarity matrix for every 100 SNPs, with each one after the first overlapping the upstream by 50 SNPs. Using the dissimilarity matrices, {maketree} implements hierarchical clustering analysis to partition strains, returning a large list of dendrograms. 
	The function {snpcor} runs correlation calculations with the phenotype dendrogram against each SNP dendrogram. {snpcor} generates a correlation matrix between all of dendrograms (phenotype and SNP), using the function {cor.dendlist} from the R package {dendextend}. {snpcor} assumes the labels in the 2 trees match as least partially. The method implemented for calculating the correlation coefficient is called Baker’s Gamma, which is a measure of similarity between 2 dendrograms. It is calculated by taking 2 items and finding the highest possible level of k (number of cluster groups created when cutting the tree) for which the 2 items still belong to the same tree. Then k is returned, and the same is done for these 2 items for the second tree. These 2 sets of numbers are paired and a generic spearman correlation is calculated. The value ranges between -1 to 1; Values near 0 means that the two trees are not statistically similar. The Baker’s Gamma Correlation Coefficient is not affected by the height of a branch but only by its relative position compared with other branches. {snpcor} returns a dataframe of the correlation coefficient between the phenotype dendrogram and each SNP dendrogram, as well as the approximate location of the SNP in basepairs. 
	The function {drawtangle} visualizes the phenotype-SNP correlation by plotting out tanglegrams, which are 2 side by side dendrograms with lines connecting identical strains.
	A proof of concept script using Biomy is written in the ExampleScript folder. It uses a dataset of metabolic phenotypes of 12 mouse strains from Mouse Phenome Database and chromosome 4 SNP data of the same mouse strains from UNC. The region of highest correlated SNPs was calculated to be 53-54 Mbp, which  is near the abca1 gene that is associated with multiple metabolic diseases.
#
##Link: 
Visit this link to download the latest package: <https://github.com/16xchen/Biomy>
#
##Install:
To install the package into R, use the following commands:

install.packages("devtools")

library(devtools)

install_github("16xchen/Biomy")

library(Biomy)
#
##Contact: 
If there is any trouble with the package, contact me at <16xchen@mdirss.org>
#
#
