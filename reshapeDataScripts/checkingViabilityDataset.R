rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Loading required packages
# -------------------------------------------------------------------
library(dplyr)
library(reshape2)
library(tidyr)
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Import and organize viability trial data
# -------------------------------------------------------------------
# -------------------------------------------------------------------
# setwd and read data files
setwd("~/Dropbox/modelsF2019/viability/")
df <- read.csv(file="viability.csv",header=TRUE)

countsCorrect<-df %>%
  dplyr::mutate(germNot2=germStart-germCount) %>%
  dplyr::mutate(checkGerm=germNot2==germNot) %>%
  dplyr::filter(isFALSE(checkGerm)) %>%
  nrow

# if germCount is correct, nrow=0
print(countsCorrect)

df %>% 
  dplyr::filter(is.na(germStart))

df %>% 
  dplyr::filter(is.na(germCount))

df %>% 
  dplyr::filter(is.na(viabStart))

df %>% 
  dplyr::filter(is.na(viabStain)) %>%
  dplyr::select(viabStart) %>% sum

df %>%
  dplyr::filter(germStart - germCount < viabStart) 

df %>%
  dplyr::filter(germNot < viabStart) 

df %>%
  dplyr::filter(viabStart < viabStain) 

