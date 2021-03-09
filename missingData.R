# countFruitsPerPlantTransects <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantTransects.RDS")
censusSeedlingsFruitingPlants <- readRDS("~/Dropbox/dataLibrary/postProcessingData/censusSeedlingsFruitingPlants.RDS")

tmp1 <- censusSeedlingsFruitingPlants %>%
  #dplyr::filter(!(seedlingNumber==0&fruitplNumber==0)) %>%
  dplyr::group_by(year,site) %>%
  dplyr::summarise(totalSeedlings = sum(seedlingNumber,na.rm=TRUE),
                   totalFruitingPlants = sum(fruitplNumber,na.rm=TRUE)) %>%
  dplyr::mutate(p.lx = totalFruitingPlants/totalSeedlings)

head(tmp1)
plot(tmp1$totalFruitingPlants/tmp1$totalSeedlings);abline(h=1)

View(tmp1)


countFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantAllPlots.RDS")
countUndamagedDamagedFruitsPerPlantAllPlots <- readRDS("~/Dropbox/dataLibrary/postProcessingData/countUndamagedDamagedFruitsPerPlantAllPlots.RDS")

tmp2 <- countFruitsPerPlantAllPlots %>%
  dplyr::group_by(year,site) %>%
  dplyr::summarise(total = sum(countFruitNumberPerPlant))

tmp3 <- countUndamagedDamagedFruitsPerPlantAllPlots %>%
  dplyr::mutate(total = sum(countUndamagedFruitNumberPerPlant,countDamagedFruitNumberPerPlant,na.rm=TRUE)) %>%
  dplyr::group_by(year,site) %>%
  dplyr::summarise(total = sum(total))

ref = tmp2 %>%
  dplyr::select(site,year)

ref3 = tmp3 %>%
  dplyr::select(site,year)

siteRef = as.factor(ref$site %>% unique)
yearRef = (as.numeric(ref$year) %>% unique)
yearRef2 = (as.numeric(ref3$year) %>% unique)
yearRef = (c(yearRef,yearRef2))

refMat=expand.grid(site=siteRef,year=yearRef)

tmpFruits = tmp2 %>%
  dplyr::bind_rows(tmp3)

fruits=refMat %>%
  dplyr::left_join(tmpFruits, by = c('site','year'))


# combine

combo=tmp1 %>%
  dplyr::select(site,year,p.lx) %>%
  dplyr::left_join(fruits,by=c('site','year'))

View(combo)

## SOME OF THE ISSUES TO ACCOUNT FOR

# 13 year*pop combinations had 0 survival and no plants outside of plots
# these are likely to be 'true zero' fitness
combo %>%
  dplyr::filter(is.na(total)&p.lx==0) 

# 7 year*pop combinations had no plants outside of the permanent plots
# these are likely to be 'true zero' fitness
combo %>%
  dplyr::filter(is.na(p.lx)&is.na(total))

# 7 year*pop combinations that had 0 seedlings in the permanent plots but
# had fruiting plants in the permanent plots and as well as outside of them
# missed, phenology-wise
combo %>%
  dplyr::filter(is.infinite(p.lx))

# 9 year*pop combinations had fruiting plants outside of the permanent plots
# some germination but missed in the plots (spatial miss)
combo %>%
  dplyr::filter(is.na(p.lx)&!is.na(total))

# 20 year*pop combinations had >0 survival but no plants outside the plots
# not sure what is up with these as there should be estimates on the fruits per plant?
# this is actually missing data from 2020
combo %>%
  dplyr::filter(is.na(total)&p.lx!=0)  %>% View


