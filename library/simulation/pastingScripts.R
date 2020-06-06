read.txt("~/Dropbox/clarkiaSeedBanks/library/simulation/priors")

library(tidyverse)
#library(readr)
mystring <- read_file("~/Dropbox/clarkiaSeedBanks/library/simulation/priors")
modelLines = readLines("~/Dropbox/clarkiaSeedBanks/library/simulation/fullModel.txt")

fileConn<-file("~/Dropbox/clarkiaSeedBanks/library/simulation/fullModelScript.R")
writeLines(modelLines, fileConn)
close(fileConn)

initLines = "model {"
priorLines = readLines("~/Dropbox/clarkiaSeedBanks/library/simulation/priors.R")
likelihoodLines = readLines("~/Dropbox/clarkiaSeedBanks/library/simulation/likelihoods.R")
derivedLines = readLines("~/Dropbox/clarkiaSeedBanks/library/simulation/deriveds.R")
tailLines = "}"

fileConn<-file("~/Dropbox/clarkiaSeedBanks/library/simulation/fullModelScript.R")
writeLines(c(initLines,priorLines,tailLines), fileConn)
close(fileConn)

"Hyperprior" %in% priorLines

section <- function(variable = "Hyperprior",vec = priorLines){
  vals <- grep(variable,priorLines)
  tmp <- priorLines[(vals[1]+1):(vals[2]-1)]
  return(tmp)
}

fileConn<-file("~/Dropbox/clarkiaSeedBanks/library/simulation/fullModelScript.R")
writeLines(c(initLines,section(variable="Hyperprior",vec = priorLines),tailLines), fileConn)
close(fileConn)
