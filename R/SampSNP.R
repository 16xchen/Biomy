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
  require("RCurl")||install.packages("RCurl")
  require("foreign") || install.packages("foreign")
  library(RCurl)
  library(foreign)
  url= "https://raw.githubusercontent.com/16xchen/Biomy/master/chrX.csv"
  chrX.data = getURL(url)
  return(read.csv(textConnection(chrX.data)))
}
