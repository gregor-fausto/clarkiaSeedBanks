
data2006 <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantAllPlots.RDS")

data2013 <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countUndamagedDamagedFruitsPerPlantAllPlots.RDS")

df.sum <- df %>% 
  dplyr::group_by(site,year) %>%
  dplyr::summarise(p = sum(countDamagedFruitNumberPerPlant>0)/n()) 

ggplot(data=df.sum) +
  geom_point(aes(x=year,y=p)) +
  facet_wrap(~site) + theme_bw() + ylim(c(0,1))

dataAll <- data2013 %>%
  dplyr::mutate(countFruitNumberPerPlant = countUndamagedFruitNumberPerPlant + countDamagedFruitNumberPerPlant, damage = 0) %>%
  dplyr::select(-c(countUndamagedFruitNumberPerPlant,countDamagedFruitNumberPerPlant)) %>%
  dplyr::bind_rows(data2006) %>%
  dplyr::group_by(site,year) %>%
  dplyr::summarise(mu = mean(countFruitNumberPerPlant)) 

ggplot(data=dataAll) +
  geom_point(aes(x=year,y=mu)) +
  facet_wrap(~site) + theme_bw() 
