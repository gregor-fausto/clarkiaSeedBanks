rm(list=ls(all=TRUE)) # clear R environment

pops=c('BG','BR','CF','CP3','DEM','DLW','EC','FR','GCN','KYE','LCE','LCW','LO',
  'MC','OKRE','OKRW','OSR','S22','SM','URS');
#pops="S22"
numpops=length(pops);

## read basic seed bag data (same as AmNat)
# seedBagMaster<-readxl::read_excel("~/Downloads/seed bag MasterG_updated.xlsx")
# 
# seedBagMaster <- seedBagMaster %>% 
#   dplyr::filter(!is.na(`intactJ`) & !is.na(`Oct total`))
# 
# names(seedBagMaster)
# Pop = pull(seedBagMaster,SITE)
# Bag = pull(seedBagMaster,`bag #`)
# Round = pull(seedBagMaster,Round)
# Age = pull(seedBagMaster,`Age at retrieval`)
# JanGerm = pull(seedBagMaster,`Jan germ`)
# JanIntact = pull(seedBagMaster,`Jan intact`)
# JanTot = pull(seedBagMaster,`Jan total`)
# OctIntact = pull(seedBagMaster,`Oct total`)

# read subset seed bag data to compare to JAGS fit
seedBagMaster <- readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/seedBagMaster.RDS")

names(seedBagMaster)
Pop = pull(seedBagMaster,site)
Bag = pull(seedBagMaster,bagNo)
Round = pull(seedBagMaster,round)
Age = pull(seedBagMaster,age)
JanGerm = pull(seedBagMaster,seedlingJan)
JanIntact = pull(seedBagMaster,intactJan)
JanTot = pull(seedBagMaster,totalJan)
OctIntact = pull(seedBagMaster,intactOct)

## read data from viability trials
# viabFinal <-readxl::read_excel("~/Downloads/germ_viab_data_final_updated.xlsx")
# 
# viabFinal <- viabFinal %>% 
#   dplyr::filter(!is.na(`#germinating`) & !is.na(`#tested_for_viability`))
# 
# PopV = pull(viabFinal,population)
# BagV = pull(viabFinal,`bag#`)
# RoundV = pull(viabFinal,Round)
# AgeV = pull(viabFinal,`Age`)
# GermStart = pull(viabFinal,`#tested_for_germination`)
# NoGerm = pull(viabFinal,`#not germinated`)
# ViabStart = pull(viabFinal,`#tested_for_viability`)
# NumStain = as.numeric(pull(viabFinal,`#viable(stained_red)`))

# read data from viability trials to compare to JAGS
viabFinal <- readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/viabilityFinal.RDS")

PopV = pull(viabFinal,site)
BagV = pull(viabFinal,bagNo)
RoundV = pull(viabFinal,round)
AgeV = pull(viabFinal,age)
GermStart = pull(viabFinal,germStart)
NoGerm = pull(viabFinal %>% dplyr::mutate(noGerm=germStart-germCount),noGerm)
ViabStart = pull(viabFinal,viabStart)
NumStain = as.numeric(pull(viabFinal,viabStain))

# convert NaNs to zeros for no. staining when none were tested
I=which(ViabStart %in% 0)
NumStain[I] = 0

# create empty matrices
s1=s2=s3=s4=s5=s6=g1=g2=g3=SSs1=SSs2=SSs3=SSs4=SSs5=SSs6=SSg1=SSg2=SSg3=NBs1=NBs2=NBs3=NBs4=NBs5=NBs6=PV1=PV2=PV3=matrix(0,nrow=numpops,ncol=3)

pv.list = list()
mat.list = list()

viabFinal<-viabFinal %>% dplyr::mutate(index = paste(site,round,age,bagNo))

mat.list = data.frame(matrix(NA,nrow=length(unique(viabFinal$index)),ncol=13))
colnames(mat.list) = c("site","round","age","bagNo","pv","PVtot","mult",'ns','vs',"gs","ng","PVmean","PVwts")
mat.list$site <- as.character(mat.list$site)
index = 1

pv1.list = data.frame(matrix(NA,nrow=60,ncol=6))
colnames(pv1.list) = c("site","round","age","PVmean","PVwts","PVratio")
index2 = 1

for(p in 1:numpops){
  
  PopNow = pops[p];
  
  ## Round 1
  # Age 1
  r = 1; a = 1; # estimate s1, g1, s2
  
   ## CALL GETVIAB
   source("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/getViab-troubleshoot.R")
  
  s2[p,r] = surv;
  SSs2[p,r] = SS;
  NBs2[p,r] = NB;
  PV1[p,r] = PVmean;
  
  #   I = which(Pop%in%PopNow & Round==r & Age==a & !is.na(JanIntact));
  #   NBs1[p,r] = length(I);
  #   SSs1[p,r] = 100*NBs1[p,r];
  #   SSg1[p,r] = sum( JanGerm[I] + (PV1[p,r]^(1/3))*JanIntact[I]);
  #   s1[p,r] = SSg1[p,r]/SSs1[p,r];
  #   g1[p,r] = sum(JanGerm[I])/SSg1[p,r];
  
  r = 2; a = 1; # estimate s1, g1, s2

  ## CALL GETVIAB
  source("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/getViab-troubleshoot.R")

  s2[p,r] = surv;
  SSs2[p,r] = SS;
  NBs2[p,r] = NB;
  PV1[p,r] = PVmean;
  
  r = 3; a = 1; # estimate s1, g1, s2

  ## CALL GETVIAB
  source("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/getViab-troubleshoot.R")

  s2[p,r] = surv;
  SSs2[p,r] = SS;
  NBs2[p,r] = NB;
  PV1[p,r] = PVmean;
  
} 

mat.list

saveRDS(pv1.list,"~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/pv1List.RDS")
saveRDS(mat.list,"~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/matList.RDS")
saveRDS(PV1,"~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/PV1.RDS")
# 
# vComp=viabFinal %>% dplyr::select(site,round,age,bagNo) %>%
#   tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
#   dplyr::filter(age==1)
# sComp=seedBagMaster %>% dplyr::select(site,round,age,bagNo) %>%
#   tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
#   dplyr::filter(age==1)
# 
# vComp$id[!(vComp$id %in% sComp$id)]
# 
# viabFinal %>% 
#   tidyr::unite(col='id', c(site,bagNo,round,age), sep="-", remove=FALSE) %>%
#   dplyr::filter(age==1) %>%
#   dplyr::filter(!(vComp$id %in% sComp$id)) %>% View

# 
#   s2[p,r] = surv;
#   SSs2[p,r] = SS;
#   NBs2[p,r] = NB;
#   PV1[p,r] = PVmean;

#   I = which(Pop%in%PopNow & Round==r & Age==a & !is.na(JanIntact));
#   NBs1[p,r] = length(I);
#   SSs1[p,r] = 100*NBs1[p,r];
#   SSg1[p,r] = sum( JanGerm[I] + (PV1[p,r]^(1/3))*JanIntact[I]);
#   s1[p,r] = SSg1[p,r]/SSs1[p,r];
#   g1[p,r] = sum(JanGerm[I])/SSg1[p,r];
#   
#   ## Round 1
#   # Age 2
#   r = 1; a = 2; # estimate s3, g2, s4
#   
#   ## CALL GETVIAB
#   source("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/getViab.R")
#   s4[p,r] = surv;
#   SSs4[p,r] = SS;
#   NBs4[p,r] = NB;
#   PV2[p,r] = PVmean;
#   
#   I = which(Pop%in%PopNow & Round==r & Age==a & !is.na(JanIntact));
#   NBs3[p,r] = length(I);
#   # estimate total number of ungerminated seeds still alive after 1 year in October
#   SSs3[p,r] = 100*s1[p,r]*(1-g1[p,r])*s2[p,r]*NBs3[p,r];
#   if( PV2[p,r] < PV1[p,r] ){
#     PVjan = (PV1[p,r]^(2/3))*(PV2[p,r]^(1/3))
#   } else {
#     PVjan = PV2[p,r]^(1/3)
#   }
#   ## SS for estimating G2 is number of VIABLE seeds in January
#   SSg2[p,r] = sum( JanGerm[I] + PVjan*JanIntact[I]);
#   s3[p,r] = min(SSg2[p,r]/SSs3[p,r], 1);
#   g2[p,r] = sum(JanGerm[I])/SSg2[p,r];
#   
#   ## Round 1
#   # Age 3
#   r = 1; a = 3; # estimate s5, g3, s6
#   
#   ## CALL GETVIAB
#   source("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/getViab.R")
#   s6[p,r] = surv;
#   SSs6[p,r] = SS;
#   NBs6[p,r] = NB;
#   PV3[p,r] = PVmean;
#   
#   I = which(Pop%in%PopNow & Round==r & Age==a & !is.na(JanIntact));
#   NBs5[p,r] = length(I);
#   # estimate total number of ungerminated seeds still alive after 2 year in October
#   SSs5[p,r] = 100*s1[p,r]*(1-g1[p,r])*s2[p,r]*s3[p,r]*(1-g2[p,r])*s4[p,r]*NBs5[p,r];
#   if( PV3[p,r] < PV2[p,r] ){
#     PVjan = (PV2[p,r]^(2/3))*(PV3[p,r]^(1/3))
#   } else {
#     PVjan = PV3[p,r]^(1/3)
#   }
#   ## SS for estimating G3 is number of VIABLE seeds in January
#   SSg3[p,r] = sum( JanGerm[I] + PVjan*JanIntact[I]);
#   s5[p,r] = min(SSg3[p,r]/SSs5[p,r], 1);
#   g3[p,r] = sum(JanGerm[I])/SSg3[p,r];
#   
#   ## Round 2
#   # Age 1
#   r = 2; a = 1; # estimate s1, g1, s2
#   
#   ## CALL GETVIAB
#   source("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/getViab.R")
#   s2[p,r] = surv;
#   SSs2[p,r] = SS;
#   NBs2[p,r] = NB;
#   PV1[p,r] = PVmean;
#   
#   I = which(Pop%in%PopNow & Round==r & Age==a & !is.na(JanIntact));
#   NBs1[p,r] = length(I);
#   SSs1[p,r] = 100*NBs1[p,r];
#   SSg1[p,r] = sum( JanGerm[I] + (PV1[p,r]^(1/3))*JanIntact[I]);
#   s1[p,r] = SSg1[p,r]/SSs1[p,r];
#   g1[p,r] = sum(JanGerm[I])/SSg1[p,r];
#   
#   ## Round 2
#   # Age 2
#   r = 2; a = 2; # estimate s3, g2, s4
#   
#   ## CALL GETVIAB
#   source("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/getViab.R")
#   s4[p,r] = surv;
#   SSs4[p,r] = SS;
#   NBs4[p,r] = NB;
#   PV2[p,r] = PVmean;
#   
#   I = which(Pop%in%PopNow & Round==r & Age==a & !is.na(JanIntact));
#   NBs3[p,r] = length(I);
#   # estimate total number of ungerminated seeds still alive after 1 year in October
#   SSs3[p,r] = 100*s1[p,r]*(1-g1[p,r])*s2[p,r]*NBs3[p,r];
#   if( PV2[p,r] < PV1[p,r] ){
#     PVjan = (PV1[p,r]^(2/3))*(PV2[p,r]^(1/3))
#   } else {
#     PVjan = PV2[p,r]^(1/3)
#   }
#   ## SS for estimating G2 is number of VIABLE seeds in January
#   SSg2[p,r] = sum( JanGerm[I] + PVjan*JanIntact[I]);
#   s3[p,r] = min(SSg2[p,r]/SSs3[p,r], 1);
#   g2[p,r] = sum(JanGerm[I])/SSg2[p,r];
#   
#   ## Round 3
#   # Age 1
#   r = 3; a = 1; # estimate s1, g1, s2
#   
#   ## CALL GETVIAB
#   source("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/getViab.R")
#   s2[p,r] = surv;
#   SSs2[p,r] = SS;
#   NBs2[p,r] = NB;
#   PV1[p,r] = PVmean;
#   
#   I = which(Pop%in%PopNow & Round==r & Age==a & !is.na(JanIntact));
#   NBs1[p,r] = length(I);
#   SSs1[p,r] = 100*NBs1[p,r];
#   SSg1[p,r] = sum( JanGerm[I] + (PV1[p,r]^(1/3))*JanIntact[I]);
#   s1[p,r] = SSg1[p,r]/SSs1[p,r];
#   g1[p,r] = sum(JanGerm[I])/SSg1[p,r];
#   
# }
# 
# 
# vitalRates = list(
#   g1=g1,g2=g2,g3=g3,s1=s1,s2=s2,s3=s3,s4=s4,s5=s5,s6=s6,
#   PV1=PV1
# )
# 
# saveRDS(vitalRates,file="~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/eckhartEstimates.RDS")
