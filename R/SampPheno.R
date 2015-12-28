


SampSNP=function(){
  if(!require("RCurl"))
    install.packages("RCurl")
  if(!require("foreign"))
    install.packages("foreign")
  library(RCurl)
  library(foreign)
  url= "https://raw.githubusercontent.com/16xchen/Biomy/master/phenodata.csv"
  pheno.data = getURL(url, .opts = list(ssl.verifypeer = FALSE))
  return(read.csv(textConnection((pheno.data))))
}
