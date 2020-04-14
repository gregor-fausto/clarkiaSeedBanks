rm(list=ls(all=TRUE)) # clear R environment
options(stringsAsFactors = FALSE)
# -------------------------------------------------------------------
# Loading libraryd packages
# -------------------------------------------------------------------
library(xlsx)
library(dplyr) # manipulate data frames
library(tidyr) # manipulate data frames
library(reshape2) # manipulate data frames
###########################################################################
# Read in data
###########################################################################
setwd("/Users/Gregor/Dropbox/Clarkia-LTREB/20_demography_sites")
seedBags<-read.xlsx(file="seed bag data.xls",sheetIndex = 1, header = TRUE)
setwd("/Users/Gregor/Dropbox/projects/clarkiaScripts/reshapeDataScripts")
viabilityEstimates<-read.xlsx(file="germ&viability.xls",sheetIndex=1,header=TRUE)

# rename variables
viabilityEstimates<-viabilityEstimates %>%
  dplyr::rename(site=population) %>%
  dplyr::rename(bagNo=bag..) %>%
  dplyr::rename(round=Round) %>%
  dplyr::rename(age=Age..yrs.) %>%
  dplyr::rename(block=Block) %>%
  dplyr::rename(germStart=Germination.starting..) %>%
  dplyr::rename(notGerm=X..not.germinated) %>%
  dplyr::rename(percGerm=X..Germination) %>%
  dplyr::rename(viabStart=Viability.starting..) %>%
  dplyr::rename(numStain=total...stained) %>%
  dplyr::rename(notGermIsViabStart=not.Germ...viab.start.) %>%
  dplyr::rename(percViab=X..Viability) %>%
  dplyr::rename(totEstViab=Total.Est..viability) %>%
  dplyr::select(-c(NA.))

viabilityEstimates$numStain<-as.numeric(as.character(viabilityEstimates$numStain))  
viabilityEstimates$percViab<-as.numeric(as.character(viabilityEstimates$percViab))  
viabilityEstimates$totEstViab<-as.numeric(as.character(viabilityEstimates$totEstViab))  

# calculate V1, viability of age(1) seeds
age1<-viabilityEstimates %>%
  dplyr::filter(age==1) %>%
  dplyr::group_by(site,bagNo,round) %>%
  dplyr::summarise(viabStart=sum(viabStart,na.rm=TRUE),
                   germStart=sum(germStart,na.rm=TRUE),
                   notGerm=sum(notGerm,na.rm=TRUE),
                   numStain=sum(numStain,na.rm=TRUE)) %>%
  dplyr::mutate(V1=ifelse(viabStart>0,
                            ((germStart-notGerm)+notGerm*(numStain/viabStart))/germStart,
                            (germStart-notGerm)/germStart))

age1<-age1 %>%
  dplyr::group_by(site,round) %>%
  dplyr::summarise(part1=sum(germStart*V1,na.rm=TRUE),part2=sum(germStart,na.rm=TRUE),weightedV1=part1/part2)

# calculate V2, viability of age(2) seeds
age2<-viabilityEstimates %>%
  dplyr::filter(age==2) %>%
  dplyr::group_by(site,bagNo,round) %>%
  dplyr::summarise(viabStart=sum(viabStart,na.rm=TRUE),
                   germStart=sum(germStart,na.rm=TRUE),
                   notGerm=sum(notGerm,na.rm=TRUE),
                   numStain=sum(numStain,na.rm=TRUE)) %>%
  dplyr::mutate(V2=ifelse(viabStart>0,
                          ((germStart-notGerm)+notGerm*(numStain/viabStart))/germStart,
                          (germStart-notGerm)/germStart))

age2<-age2 %>%
  dplyr::group_by(site,round) %>%
  dplyr::summarise(part1=sum(germStart*V2,na.rm=TRUE),part2=sum(germStart,na.rm=TRUE),weightedV2=part1/part2)

# calculate V3, viability of age(3) seeds
age3<-viabilityEstimates %>%
  dplyr::filter(age==3) %>%
  dplyr::group_by(site,bagNo,round) %>%
  dplyr::summarise(viabStart=sum(viabStart,na.rm=TRUE),
                   germStart=sum(germStart,na.rm=TRUE),
                   notGerm=sum(notGerm,na.rm=TRUE),
                   numStain=sum(numStain,na.rm=TRUE)) %>%
  dplyr::mutate(V3=ifelse(viabStart>0,
                          ((germStart-notGerm)+notGerm*(numStain/viabStart))/germStart,
                          (germStart-notGerm)/germStart))

age3<-age3 %>%
  dplyr::group_by(site,round) %>%
  dplyr::summarise(part1=sum(germStart*V3,na.rm=TRUE),part2=sum(germStart,na.rm=TRUE),weightedV3=part1/part2)

# bind all
age1$age="one"
age2$age="two"
age3$age="three"

age1<-age1 %>% 
  dplyr::select(c(site,round,weightedV1,age)) %>%
  dplyr::rename(viab=weightedV1)
age2<-age2 %>% 
  dplyr::select(c(site,round,weightedV2,age)) %>%
  dplyr::rename(viab=weightedV2)
age3<-age3 %>% 
  dplyr::select(c(site,round,weightedV3,age)) %>%
  dplyr::rename(viab=weightedV3)

viability<-rbind(age1,age2,age3)

#V.33
age1<-viability %>% filter(age=="one") %>%
  dplyr::mutate(V.33=viab^(1/3)) %>% 
  dplyr::rename(weightedV1=viab) %>% select(-age)

#V1.33
round1<-viability %>% filter(round==1&age!="three")
round2<-viability %>% filter(round==2&age!="three")

round1<-tidyr::spread(round1,age,viab) %>% 
  dplyr::mutate(V1.33=ifelse(two<one,one*(two/one)^(1/3),two^(1/3)))
round2<-tidyr::spread(round2,age,viab) %>%
  dplyr::mutate(V1.33=ifelse(two<one,one*(two/one)^(1/3),two^(1/3)))

age2<-rbind(round1,round2) %>% 
  dplyr::select(-one) %>% 
  dplyr::rename(weightedV2=two) 

#V2.33
round1<-viability %>% filter(round==1&age!="one")

round1<-tidyr::spread(round1,age,viab) %>%
  dplyr::mutate(V2.33=ifelse(three<two,two*(three/two)^(1/3),three^(1/3)))

age3<-round1 %>% dplyr::select(-two) %>% dplyr::rename(weightedV3=three)

# rename seed bag variables
seedBags<-seedBags %>%
  dplyr::rename(bagNo=bag..) %>%
  dplyr::rename(round=Round) %>%
  dplyr::rename(startYear=Start_yr.Oct.) %>%
  dplyr::rename(age=Age) %>%
  dplyr::rename(viabilityEstYear=Yr_germ.Jan..viability.Oct.) %>%
  dplyr::select(-c(viability.))

# age 1 bags
rates1<-seedBags %>% 
  dplyr::filter(age==1) %>%
  dplyr::left_join(age1,by=c("site","round"))

s1<-rates1 %>% 
  dplyr::filter(!is.na(Jan_intact) & !is.na(Jan_germ)) %>%
  dplyr::mutate(S=Jan_germ+Jan_intact*V.33,N=ifelse(!is.na(Jan_germ)&!is.na(Jan_intact),100,NA)) %>% 
  dplyr::group_by(site,round) %>%
  dplyr::summarise(s1=sum(S,na.rm=TRUE)/sum(N,na.rm=TRUE))

round(spread(s1,round,s1)[,2:4],digits=4)

g1 <- rates1 %>%
  dplyr::filter(!is.na(Jan_intact) & !is.na(Jan_germ)) %>%
  dplyr::mutate(S=Jan_germ,N=Jan_germ+Jan_intact*V.33) %>%
  dplyr::group_by(site,round) %>%
  dplyr::summarise(g1=sum(S,na.rm=TRUE)/sum(N,na.rm=TRUE))

# for calculating s2, need to filter more
# see file for 'Estimation of Clarkia vital rates.doc'
rates1$bound <- paste(rates1$site, rates1$bagNo, sep = "-")
viabilityEstimates$bound <- paste(viabilityEstimates$site,viabilityEstimates$bagNo, sep="-")

s2 <- rates1[ which( rates1$bound %in% viabilityEstimates$bound), ] %>%
  dplyr::filter(!is.na(Jan_intact) & !is.na(Oct_intact)) %>%
  dplyr::mutate(S=Oct_intact,N=Jan_intact*V.33) %>%
  dplyr::group_by(site,round) %>%
  dplyr::summarise(s2=sum(S,na.rm=TRUE)/sum(N,na.rm=TRUE))
  
s0<-s2
s0$s0 <- s0$s2^(3/8)
s0<-s0[,-3]

# age 2 bags
rates2<-seedBags %>% 
  dplyr::filter(age==2) %>% 
  dplyr::left_join(age2,by=c("site","round")) %>%
  dplyr::left_join(s0,by=c("site","round")) %>%
  dplyr::left_join(s1,by=c("site","round")) %>%
  dplyr::left_join(s2,by=c("site","round")) %>% 
  dplyr::left_join(g1,by=c("site","round"))

s3<-rates2 %>% 
  dplyr::filter(!is.na(Jan_intact) & !is.na(Jan_germ)) %>%
  dplyr::mutate(S=Jan_germ+Jan_intact*V1.33,N=100*s1*(1-g1)*s2) %>% 
  dplyr::group_by(site,round) %>%
  dplyr::summarise(s3=sum(S,na.rm=TRUE)/sum(N,na.rm=TRUE)) %>%
  dplyr::mutate(s3=ifelse(s3<=1,s3,1))

round(spread(s3,round,s3)[,2:4],digits=4)

g2 <- rates2 %>%
  dplyr::filter(!is.na(Jan_intact) & !is.na(Jan_germ)) %>%
  dplyr::mutate(S=Jan_germ,N=Jan_germ+Jan_intact*V1.33) %>%
  dplyr::group_by(site,round) %>%
  dplyr::summarise(g2=sum(S,na.rm=TRUE)/sum(N,na.rm=TRUE))

# for calculating s4, need to filter more
# ACTUALLY, check this with Bill
# see file for 'Estimation of Clarkia vital rates.doc'
rates2$bound <- paste(rates2$site, rates2$bagNo, sep = "-")
viabilityEstimates$bound <- paste(viabilityEstimates$site,viabilityEstimates$bagNo, sep="-")

s4 <- rates2[ which( rates2$bound %in% viabilityEstimates$bound), ] %>%
  dplyr::filter(!is.na(Jan_intact) & !is.na(Oct_intact)) %>%
  dplyr::mutate(S=Oct_intact,N=Jan_intact*V1.33) %>%
  dplyr::group_by(site,round) %>%
  dplyr::summarise(s4=sum(S,na.rm=TRUE)/sum(N,na.rm=TRUE)) %>%
  dplyr::mutate(s4=ifelse(s4<=1,s4,1))

# age 3 bags
rates3<-seedBags %>% 
  dplyr::filter(age==3) %>% 
  dplyr::left_join(age3,by=c("site","round")) %>%
  dplyr::left_join(s0,by=c("site","round")) %>%
  dplyr::left_join(s1,by=c("site","round")) %>%
  dplyr::left_join(s2,by=c("site","round")) %>% 
  dplyr::left_join(g1,by=c("site","round")) %>%
  dplyr::left_join(s3,by=c("site","round")) %>%
  dplyr::left_join(s4,by=c("site","round")) %>% 
  dplyr::left_join(g2,by=c("site","round"))

s5<-rates3 %>% 
  dplyr::filter(!is.na(Jan_intact) & !is.na(Jan_germ)) %>%
  dplyr::mutate(S=Jan_germ+Jan_intact*V2.33,N=100*s1*(1-g1)*s2*s3*(1-g2)*s4) %>% 
  dplyr::group_by(site,round) %>%
  dplyr::summarise(s5=sum(S,na.rm=TRUE)/sum(N,na.rm=TRUE)) %>%
  dplyr::mutate(s5=ifelse(s5<=1,s5,1))

g3 <- rates3 %>%
  dplyr::filter(!is.na(Jan_intact) & !is.na(Jan_germ)) %>%
  dplyr::mutate(S=Jan_germ,N=Jan_germ+Jan_intact*V2.33) %>%
  dplyr::group_by(site,round) %>%
  dplyr::summarise(g3=sum(S,na.rm=TRUE)/sum(N,na.rm=TRUE))

# for calculating s6, need to filter more
# ACTUALLY, check this with Bill
# see file for 'Estimation of Clarkia vital rates.doc'
rates3$bound <- paste(rates3$site, rates3$bagNo, sep = "-")
viabilityEstimates$bound <- paste(viabilityEstimates$site,viabilityEstimates$bagNo, sep="-")

s6 <- rates3[ which( rates3$bound %in% viabilityEstimates$bound), ] %>%
  dplyr::filter(!is.na(Jan_intact) & !is.na(Oct_intact)) %>%
  dplyr::mutate(S=Oct_intact,N=Jan_intact*V2.33) %>%
  dplyr::group_by(site,round) %>%
  dplyr::summarise(s6=sum(S,na.rm=TRUE)/sum(N,na.rm=TRUE)) %>%
  dplyr::mutate(s6=ifelse(s6<=1,s6,1))

###########################################################################
# create summary table
###########################################################################
summarizeParams<-function(df){
  name <- substitute(df) %>%
    deparse()
  names(df)[3] <- "estimate"
  estimateVec <- df %>% 
    dplyr::group_by(site) %>%
    dplyr::summarise(estimate=mean(estimate,na.rm=TRUE))
  colnames(estimateVec)[2] <- name[1]
  estimateVec<-as.data.frame(estimateVec)
  return(estimateVec)
}
names<-c("s0","s1","s2","s3","s4","s5","s6","g1","g2","g3")
l<-list(s0,s1,s2,s3,s4,s5,s6,g1,g2,g3)
for(i in 1:10){
  l[[i]]<-summarizeParams(l[[i]])
  colnames(l[[i]])[2] <- names[i]
}

bound<-do.call(cbind, l) 
bound<-cbind(bound[1],bound[2],bound[4],bound[6],bound[8],bound[10],bound[12],bound[14],bound[16],bound[18],bound[20])

setwd("/Users/Gregor/Dropbox/projects/clarkiaScripts/reshapeData")
siteAbiotic<-read.csv("siteAbiotic.csv",header=TRUE) %>%
  select(site,easting)

estimates<-bound %>% left_join(siteAbiotic,by="site") %>%
  arrange(easting)

setwd("/Users/Gregor/Dropbox/projects/clarkiaScripts")
write.csv(estimates,"estimates.csv")

###########################################################################
# compare estimates and original data
###########################################################################
setwd("/Users/Gregor/Dropbox/projects/clarkiaScripts/originalAnalysis")
tableA2<-read.csv(file="tablea2.csv")
tableA2<-tableA2[,1:11]

estimates<-cbind(site=estimates[,1],round(estimates[,2:11],digits=2))

tableA2<-tableA2[order(match(tableA2[,1],estimates[,1])),]

percentDiff<-cbind(site=estimates[,1],abs((estimates[,2:11]-tableA2[,2:11])/estimates[,2:11]*100))
write.csv(percentDiff,"difference.csv")

# CSV files
setwd("~/Dropbox/projects/clarkiaScripts/reshapeData")
raw.count.list<-list(read.csv("s1.csv",header=TRUE),
                             read.csv("s2.csv",header=TRUE),
                             read.csv("s3.csv",header=TRUE),
                             read.csv("s4.csv",header=TRUE),
                             read.csv("g1.csv",header=TRUE),
                             read.csv("g2.csv",header=TRUE))

list.of.data.frames<-list()

estimateList<-c("s1","s2","s3","s4","g1","g2")

for(i in 1:length(raw.count.list)){
  demData<-raw.count.list[[i]]
  list.of.data.frames[[i]]<-demData %>%
    dplyr::select(site, year, S, N) %>%
    dplyr::mutate(estimate=S/N) %>%
    dplyr::select(-c(S,N))
}

seedVitalEstimates = Reduce(function(...) dplyr::left_join(...,by=c("site","year")), list.of.data.frames)

names(seedVitalEstimates)[3:8]<-c("s1","s2","s3","s4","g1","g2")

rawCounts <- seedVitalEstimates %>% 
  dplyr::group_by(site) %>%
  dplyr::summarise(s1=mean(s1,na.rm=TRUE),
                   s2=mean(s2,na.rm=TRUE),
                   s3=mean(s3,na.rm=TRUE),
                   s4=mean(s4,na.rm=TRUE),
                   g1=mean(g1,na.rm=TRUE),
                   g2=mean(g2,na.rm=TRUE)) %>%
  as.data.frame()

rawCounts<-rawCounts[order(match(rawCounts[,1],estimates[,1])),]

estimatesParsed<-cbind(site=estimates[,1],estimates[,3:6],estimates[,9:10])

percentDiff2<-cbind(site=estimatesParsed[,1],abs((estimatesParsed[,2:7]-rawCounts[,2:7])/estimatesParsed[,2:7]*100))
write.csv(percentDiff2,"difference2.csv")

###########################################################################
# create summary Rdata objects
###########################################################################
setwd("~/Dropbox/projects/clarkiaScripts/reshapeData")
save(s0, file = "s0.RData")
save(s1, file = "s1.RData")
save(s2, file = "s2.RData")
save(s3, file = "s3.RData")
save(s4, file = "s4.RData")
save(s5, file = "s5.RData")
save(s6, file = "s6.RData")

save(g1, file = "g1.RData")
save(g2, file = "g2.RData")
save(g3, file = "g3.RData")