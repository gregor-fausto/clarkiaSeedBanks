library(tidyverse)

df<-readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantTransects.RDS")
head(df)


df2<-readRDS("~/Dropbox/dataLibrary/postProcessingData/countUndamagedDamagedFruitsPerPlantTransects.RDS")
head(df2)

referenceTable <- df %>%
  dplyr::select(site,transect,position) %>%
  unique %>%
  dplyr::arrange(position) %>%
  dplyr::arrange(transect) %>%
  dplyr::arrange(site)

referenceTable2 <- df2 %>%
  dplyr::select(site,transect,position) %>%
  unique %>%
  dplyr::arrange(position) %>%
  dplyr::arrange(transect) %>%
  dplyr::arrange(site)

dim(referenceTable);dim(referenceTable2)
referenceTable<-union(referenceTable,referenceTable2)
dim(referenceTable)

years<-rep(2007:2018,each=dim(referenceTable)[1])

referenceTable=do.call("rbind", replicate(length(2007:2018), referenceTable, simplify = FALSE)) 

referenceTable = cbind(referenceTable,year=years)

plantCounts <- df %>% 
  dplyr::group_by(year,site,transect,position) %>%
  dplyr::summarise(count = n())

plantCounts2 <- df2 %>% 
  dplyr::group_by(year,site,transect,position) %>%
  dplyr::summarise(count = n())

plantCounts <- plantCounts %>%
  dplyr::bind_rows(plantCounts2)

fullDataset = referenceTable %>%
  dplyr::full_join(plantCounts,by=c("site","transect","position","year"))

presenceAbsence = fullDataset %>%
  dplyr::mutate(presentBinary = ifelse(is.na(count),0,1))

siteNames = unique(presenceAbsence$site)

pdf("~/Desktop/trajectories.pdf")
for( i in 1:20){
print(ggplot(presenceAbsence %>% dplyr::filter(site==siteNames[i])) +
  geom_point(aes(x=year,y=interaction(site,position,transect),color=as.factor(presentBinary))) +
  theme_bw() + facet_wrap(~site) +
  scale_color_manual(values=c("white", "black")) +
  theme(legend.position = "none"))
}
dev.off()  

presenceAbsenceAverage = fullDataset %>%
  dplyr::mutate(presentBinary = ifelse(is.na(count),0,1)) %>%
  dplyr::group_by(year,site) %>%
  dplyr::summarise(sum = sum(presentBinary)/n())

presenceAbsenceTransect = fullDataset %>%
  dplyr::mutate(presentBinary = ifelse(is.na(count),0,1)) %>%
  dplyr::group_by(year,site,transect) %>%
  dplyr::summarise(sum = sum(presentBinary)/n())

ggplot(presenceAbsenceAverage ) +
  geom_point(aes(x=year,y=sum)) +
  facet_wrap(~site)

ggplot( ) +
  geom_line(data=presenceAbsenceTransect,aes(x=year,y=sum,group=transect),alpha=0.25) +
  geom_line(data=presenceAbsenceAverage,aes(x=year,y=sum)) +
  facet_wrap(~site) +
  theme_bw()

## ALL PLOTS