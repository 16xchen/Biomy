#' Loading a sample SNP data
#'
#' dataset originally downloaded from <http://msub.csbio.unc.edu/>
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
    install.packages("RCurl")
  if(!require("foreign"))
    install.packages("foreign")
  library(RCurl)
  library(foreign)
  url= "https://raw.githubusercontent.com/16xchen/Biomy/master/chrX.csv"
  chrX.data = getURL(url, .opts = list(ssl.verifypeer = FALSE))
  return(read.csv(textConnection(chrX.data)))
}
