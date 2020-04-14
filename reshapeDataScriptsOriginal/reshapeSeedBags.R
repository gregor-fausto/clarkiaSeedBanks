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
viabilityEst<-read.xlsx(file="germ&viability.xls",sheetIndex=1,header=TRUE)

# calculating V1
viabilityAgeOne<-viabilityEst %>%
  dplyr::filter(Age..yrs.==1) %>% 
  dplyr::group_by(population, bag..) %>%
  dplyr::summarise(germStart=sum(Germination.starting..),
                     noGerm=sum(X..not.germinated),
                     viabStart=sum(Viability.starting..),
                     stained=sum(as.numeric(as.character(total...stained)))) %>%
  dplyr::mutate(stained=ifelse(is.na(stained),0,stained)) %>%
  dplyr::mutate(viab=ifelse(viabStart>0,
                            ((germStart-noGerm)+noGerm*(stained/viabStart))/germStart,
                            (germStart-noGerm)/germStart)) %>%
  dplyr::summarise(part1=sum(germStart*viab),part2=sum(germStart),weightedV1=part1/part2) %>%
  dplyr::mutate(V.33=(weightedV1)^(1/3))

weights<-viabilityEst %>%
  dplyr::group_by(population,Round,Age..yrs.,bag..) %>%
  dplyr::summarise(germStart=sum(Germination.starting..),
                   noGerm=sum(X..not.germinated),
                   viabStart=sum(Viability.starting..),
                   stained=sum(as.numeric(as.character(total...stained)))) %>%
  dplyr::mutate(stained=ifelse(is.na(stained),0,stained)) %>%
  dplyr::mutate(viab=ifelse(viabStart>0,
                            ((germStart-noGerm)+noGerm*(stained/viabStart))/germStart,
                            1)) %>%
  dplyr::rename(site=population) %>%
  dplyr::rename(bagNo=bag..) %>%
  dplyr::rename(Age=Age..yrs.) %>%
  dplyr::mutate(bagNo=as.numeric(bagNo)) %>%
  dplyr::group_by(site,Round,Age) %>%
  dplyr::summarise(part1=sum(germStart*viab),part2=sum(germStart),V=part1/part2)

weights<-weights %>%
  dplyr::select(-c(part1,part2)) %>% 
  #tidyr::spread(site,V)
  dplyr::group_by(site,Round,Age) %>%
  dplyr::summarise(viability=mean(V,na.rm=TRUE))

weights$startYear=ifelse(weights$Round==1,2006,ifelse(weights$Round==2,2007,2008))
 
#rename variables
seedBags<-seedBags %>%
  dplyr::rename(bagNo=bag..) %>%
  dplyr::rename(viabilityEstYear=Yr_germ.Jan..viability.Oct.) %>%
  dplyr::rename(startYear=Start_yr.Oct.) %>%
  dplyr::left_join(weights) %>%
  dplyr::select(-viability.)

###########################################################################
# (V1) Viability at the end of year 1 
###########################################################################
V1<-seedBags %>%
  dplyr::filter(Age==1) %>% 
  dplyr::group_by(site,startYear) %>%
  dplyr::summarise(V1=mean(viability,na.rm=TRUE))

###########################################################################
# (V2) Viability at the end of year 2
###########################################################################
V2<-seedBags %>%
  dplyr::filter(Age==2) %>%
  dplyr::group_by(site,startYear) %>%
  dplyr::summarise(V2=mean(viability,na.rm=TRUE))

###########################################################################
# (V3) Viability at the end of year 3
###########################################################################
V3<-seedBags %>%
  dplyr::filter(Age==3) %>%
  dplyr::group_by(site,startYear) %>%
  dplyr::summarise(V3=mean(viability,na.rm=TRUE))

###########################################################################
# Interpolation
###########################################################################
viabs<-full_join(V1,V2) %>%
  full_join(V3)

viabs$V1<-ifelse(is.na(viabs$V2),viabs$V1,ifelse(viabs$V2<viabs$V1,viabs$V1*(viabs$V2/viabs$V1)^(1/3),(viabs$V2)^(1/3)))
viabs$V1<-ifelse(is.na(viabs$V3),viabs$V1,ifelse(viabs$V3<viabs$V1,viabs$V1*(viabs$V3/viabs$V1)^(1/6),(viabs$V3)^(1/6)))
viabs$V2<-ifelse(is.na(viabs$V3),viabs$V2,ifelse(viabs$V3<viabs$V2,viabs$V2*(viabs$V3/viabs$V2)^(1/2),(viabs$V3)^(1/2)))

V1<-select(viabs,-c(V2,V3)) %>% filter(!is.na(V1))
V2<-select(viabs,-c(V1,V3)) %>% filter(!is.na(V2))
V3<-select(viabs,-c(V1,V2)) %>% filter(!is.na(V3))

###########################################################################
# (s1) counts for seed survivorship from Oct (year 0) to Jan (year 1)
###########################################################################
s1<-seedBags %>%
  dplyr::filter(Age==1) %>% 
  dplyr::mutate(initialSeedNo=100) %>%
  dplyr::left_join(V1,by=c("site","startYear")) %>% 
  dplyr::mutate(S=(((V1)^(1/3))*Jan_intact)+Jan_germ) %>% 
  dplyr::group_by(site,startYear) %>%
  dplyr::summarise(N=sum(initialSeedNo),S=ceiling(sum(S,na.rm=TRUE)))
  
###########################################################################
# (g1) counts for germination from Oct (year 0) to Jan (year 1)
###########################################################################
g1 <- seedBags %>%
  dplyr::filter(Age==1) %>%
  dplyr::group_by(site,startYear) %>%
  dplyr::summarise(S=sum(Jan_germ,na.rm=TRUE)) %>%
  dplyr::left_join(s1,by=c("site","startYear")) %>%
  dplyr::select(-N) %>% 
  dplyr::rename(S=S.x,N=S.y) %>%
  dplyr::select(site,startYear,N,S)

###########################################################################
# (s2) counts for survival from Jan (year 1) to Oct (year 1)
###########################################################################  
s2 <- seedBags %>%
  dplyr::filter(Age==1) %>%
  dplyr::left_join(V1,by=c("site","startYear")) %>%
  dplyr::mutate(N=(((V1)^(1/3))*Jan_intact),S=(V1)*(Oct_intact)) %>%
  dplyr::group_by(site,startYear) %>%
  dplyr::summarise(N=ceiling(sum(N,na.rm=TRUE)),S=ceiling(sum(S,na.rm=TRUE)))

###########################################################################
# (s0) counts for survival from July (year 0) to Oct (year 0)
###########################################################################  
s0<-s2 %>%
  dplyr::mutate(s0=S/N) %>%
  dplyr::group_by(site) %>%
  dplyr::summarise(s0Mean = mean(s0), s0SD = sd(s0))

###########################################################################
# (s3) counts for seed survivorship from Oct (year 1) to Jan (year 2)
###########################################################################
s3<-seedBags %>%
  dplyr::filter(Age==2) %>% 
  dplyr::left_join(V1,by=c("site","startYear")) %>% 
  dplyr::left_join(V2,by=c("site","startYear")) %>% 
  dplyr::mutate(S=(((ifelse(V2/V1>1,1,V2/V1))^(1/3))*(V1)*Jan_intact)+Jan_germ) %>% 
  dplyr::group_by(site,startYear) %>%
  dplyr::summarise(S=ceiling(sum(S,na.rm=TRUE))) %>%
  dplyr::left_join(s2,by=c("site","startYear")) %>%
  dplyr::select(-N) %>% 
  dplyr::rename(S=S.x,N=S.y) %>%
  dplyr::select(site,startYear,N,S)
  
# adjust cases where S > N
s3<-s3 %>%
  dplyr::mutate(N=ifelse(N<S,S,N))

###########################################################################
# (g2) counts for germination from Oct (year 1) to Jan (year 2)
###########################################################################
g2 <- seedBags %>%
  dplyr::filter(Age==2) %>%
  dplyr::group_by(site,startYear) %>%
  dplyr::summarise(S=sum(Jan_germ,na.rm=TRUE)) %>%
  dplyr::left_join(s3,by=c("site","startYear")) %>%
  dplyr::select(-N) %>% 
  dplyr::rename(S=S.x,N=S.y) %>%
  dplyr::select(site,startYear,N,S)

###########################################################################
# (s4) counts for survival from Jan (year 2) to Oct (year 2)
###########################################################################  
s4 <- seedBags %>%
  dplyr::filter(Age==2) %>%
  dplyr::group_by(site,startYear) %>%
  dplyr::left_join(V1,by=c("site","startYear")) %>% 
  dplyr::left_join(V2,by=c("site","startYear")) %>%
  dplyr::mutate(N=(((ifelse(V2/V1>1,1,V2/V1))^(1/3))*(V1)*Jan_intact),S=(V2)*(Oct_intact)) %>%
  dplyr::summarise(N=ceiling(sum(N,na.rm=TRUE)),S=ceiling(sum(S,na.rm=TRUE)))

# adjust cases where S > N
s4<-s4 %>%
  dplyr::mutate(N=ifelse(N<S,S,N))

###########################################################################
# (s5) counts for seed survivorship from Oct (year 1) to Jan (year 2)
###########################################################################
s5<-seedBags %>%
  dplyr::filter(Age==3) %>% 
  dplyr::left_join(V2,by=c("site","startYear")) %>% 
  dplyr::left_join(V3,by=c("site","startYear")) %>% 
  dplyr::mutate(S=(((ifelse(V3/V2>1,1,V3/V2))^(1/3))*(V2)*Jan_intact)+Jan_germ) %>% 
  dplyr::group_by(site,startYear) %>%
  dplyr::summarise(S=ceiling(sum(S,na.rm=TRUE))) %>%
  dplyr::left_join(s4,by=c("site","startYear")) %>%
  dplyr::select(-N) %>% 
  dplyr::rename(S=S.x,N=S.y) %>%
  dplyr::select(site,startYear,N,S)

# adjust cases where S > N
s5<-s5 %>%
  dplyr::mutate(N=ifelse(N<S,S,N))

###########################################################################
# (g3) counts for germination from Oct (year 1) to Jan (year 2)
###########################################################################
g3 <- seedBags %>%
  dplyr::filter(Age==3) %>%
  dplyr::group_by(site,startYear) %>%
  dplyr::summarise(S=sum(Jan_germ,na.rm=TRUE)) %>%
  dplyr::left_join(s5,by=c("site","startYear")) %>%
  dplyr::select(-N) %>% 
  dplyr::rename(S=S.x,N=S.y) %>%
  dplyr::select(site,startYear,N,S)

###########################################################################
# (s6) counts for survival from Jan (year 2) to Oct (year 2)
###########################################################################  
s6 <- seedBags %>%
  dplyr::filter(Age==3) %>%
  dplyr::group_by(site,startYear) %>%
  dplyr::left_join(V2,by=c("site","startYear")) %>% 
  dplyr::left_join(V3,by=c("site","startYear")) %>%
  dplyr::mutate(N=(((ifelse(V3/V2>1,1,V3/V2))^(1/3))*(V2)*Jan_intact),S=(V3)*(Oct_intact)) %>%
  dplyr::summarise(N=ceiling(sum(N,na.rm=TRUE)),S=ceiling(sum(S,na.rm=TRUE)))

# adjust cases where S > N
s6<-s6 %>%
  dplyr::mutate(N=ifelse(N<S,S,N))

###########################################################################
# create summary table
###########################################################################
summarizeParams<-function(df){
 name <- substitute(df) %>%
    deparse()
 estimateVec <- df %>% 
   mutate(estimate = S/N) %>% 
   group_by(site) %>%
   summarise(estimate=mean(estimate))
colnames(estimateVec)[2] <- name[1]
estimateVec<-as.data.frame(estimateVec)
    return(estimateVec)
}
names<-c("s1","s2","s3","s4","s5","s6","g1","g2","g3")
l<-list(s1,s2,s3,s4,s5,s6,g1,g2,g3)
for(i in 1:9){
  l[[i]]<-summarizeParams(l[[i]])
  colnames(l[[i]])[2] <- names[i]
}

bound<-do.call(cbind, l) 
bound<-cbind(bound[1],bound[2],bound[4],bound[6],bound[8],bound[10],bound[12],bound[14],bound[16],bound[18])

setwd("/Users/Gregor/Dropbox/projects/clarkiaScripts/reshapeData")
siteAbiotic<-read.csv("siteAbiotic.csv",header=TRUE) %>%
  select(site,easting)

estimates<-bound %>% left_join(siteAbiotic,by="site") %>%
  arrange(easting)

setwd("/Users/Gregor/Dropbox/projects/clarkiaScripts")
write.csv(estimates,"estimates.csv")

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

