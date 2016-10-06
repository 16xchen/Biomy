#' Loading a sample phenotype data
#'
#' dataset is originally downloaded from <http://msub.csbio.unc.edu/>
#'
#' dataset is caloric intake data publsihed by Bodnar1 (2006)
#'
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
SampSNP=function(){
  if(!require("RCurl"))
    install.packages("RCurl",  repos = 'http://cran.rstudio.com/')
  if(!require("foreign"))
    install.packages("foreign")
  library(RCurl)
  library(foreign)
  url= "https://raw.githubusercontent.com/16xchen/Biomy/master/chr4.csv"
  chr4.data = getURL(url, .opts = list(ssl.verifypeer = FALSE))
  return(read.csv(textConnection(chr4.data)))
}
