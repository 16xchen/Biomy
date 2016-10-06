#' Loading a sample phenotype data
#'
#' dataset originally downloaded from <http://phenome.jax.org/>
#'
#' @return a large dataframe of SNP data from various mouse strains
#' @examples
#' chrX=SampSNP()
#' head(chrX)
#'
#' @author Xingyao Chen
#' @export
#'
#'


SampPheno=function(){
  if(!require("RCurl"))
    install.packages("RCurl",  repos = 'http://cran.rstudio.com/')
  if(!require("foreign"))
    install.packages("foreign")
  library(RCurl)
  library(foreign)
  url= "https://raw.githubusercontent.com/16xchen/Biomy/master/phenodata.csv"
  pheno.data = getURL(url, .opts = list(ssl.verifypeer = FALSE))
  return(read.csv(textConnection((pheno.data))))
}
