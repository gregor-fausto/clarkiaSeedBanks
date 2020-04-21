### Seed rain and seedlings
\iffalse
```{r}
## need to turn this first part into RDS file
## Seedlings from transects 

### Census 1
dfSummary<-censusSeedlingsFruitingPlants %>%
  dplyr::select(site,transect,position,year,contains(c("seedling")))

## Fruits per plant from transects
countFruitsPerPlantTransects<-readRDS("~/Dropbox/dataLibrary/postProcessingData/countFruitsPerPlantTransects.RDS")

fruitSummary <- countFruitsPerPlantTransects %>%
  dplyr::group_by(site,year,transect,position) %>%
  dplyr::summarise(countFruitsPerPlot = sum(countFruitsPerPlant),
                   n = n())
```

```{r}
# add 1 to the year in the fruit summary dataset to match the seedling year
# fruiting plants in year t match seedlings in year t+1

years = 2009

df <- list()

for(i in 1:length(years)){
  tmp <-  dfSummary %>%
    dplyr::filter(year==years[i]) %>%
    dplyr::rename(n = seedlingNumber) %>%
    dplyr::mutate(year=2) %>%
    dplyr::filter(year==2 & n>0) %>%
    dplyr::mutate(position = as.double(position))
  
  tmp1 <- fruitSummary %>%
    dplyr::filter(year==years[[i]]-1) %>%
    dplyr::ungroup(year) %>%
    dplyr::mutate(year=1,
                  n = ifelse(is.na(n),0,n)) %>%
    dplyr::select(site,year,transect,position,n) 
  
  tmp2 <- fruitSummary %>%
    dplyr::filter(year==years[[i]]-2) %>%
    dplyr::ungroup(year) %>%
    dplyr::mutate(year=0,
                  n = ifelse(is.na(n),0,n)) %>%
    dplyr::select(site,year,transect,position,n)
  
  df[[i]] <- tmp %>% 
    dplyr::left_join(tmp1,by=c("site","transect","position")) %>% 
    dplyr::left_join(tmp2,by=c("site","transect","position")) 
  
}

# ggplot(df[[1]] ) + 
#   geom_line(aes(x=year,y=n,group=interaction(site,transect,position)),alpha=.5) +
#   facet_wrap(~site,scale='free') + theme_bw()

```

One option would be to run a model with a incorporation rate from one year old seeds as well as a parameter for a number of seeds that come from elsewhere. For sites where the 1 year old seeds make most of the contribution, the second parameter would be small. 

```{r}
fruitSummary2 = fruitSummary
fruitSummary2$year = fruitSummary2$year+2
fruitSummary2<-fruitSummary2 %>% dplyr::rename(countFruitsPerPlot2 = countFruitsPerPlot,
                                               n2=n)

fruitSummary$year=fruitSummary$year+1


# set NA = 0
joinedFruitSeedling<-dfSummary %>% mutate(position = as.double(position)) %>%
  dplyr::full_join(fruitSummary,by=c("site","year","transect","position")) %>%
  dplyr::full_join(fruitSummary2,by=c("site","year","transect","position")) %>% 
  #dplyr::filter(year>2008) %>%
  dplyr::mutate(countFruitsPerPlot2 = ifelse(is.na(countFruitsPerPlot2),0,countFruitsPerPlot2),
                countFruitsPerPlot = ifelse(is.na(countFruitsPerPlot),0,countFruitsPerPlot),
                seedlingCount = ifelse(is.na(seedlingNumber),0,seedlingNumber))

ggplot(data=dfSummary %>% dplyr::filter(year>2006&year<2009)) +
  geom_line(aes(x=as.factor(year),y=seedlingNumber,group=interaction(site,transect,position)),alpha=0.5) +
  theme_bw()
```

```{r}
# create a ratio of proportion of plots with seeds the year before
joinedFruitSeedlingSummary<-joinedFruitSeedling %>%
  dplyr::left_join(position,by="site") %>%
  dplyr::mutate(test = ifelse(countFruitsPerPlot==0&countFruitsPerPlot2&seedlingCount>0,1,0)) %>%
  group_by(site,year,easting) %>%
  dplyr::summarise(p=sum(test)/n()) 

ggplot(data=joinedFruitSeedlingSummary) +
  geom_histogram(aes(p)) +
  facet_wrap(~year,scales='free_y') + xlim(c(0,1))

ggplot() +
  geom_point(data=joinedFruitSeedlingSummary,aes(x=easting,y=p),size=.5,col='orange',alpha=0.5) +
  geom_point(data=joinedSummary,aes(x=easting,y=phat),size=1,col='black',alpha=0.5) +
  geom_smooth(data=joinedSummary,aes(x=easting,y=phat),method="lm") +
  facet_wrap(~year,scales='free',nrow=2) +
  ylab("Ratio of seedlings to seeds per plot") +
  theme_bw()

ggplot() +
  geom_point(data=joinedFruitSeedling,aes(x=year,y=p),size=.5,col='orange',alpha=0.5) +
  geom_point(data=joinedSummary,aes(x=year,y=phat),size=1,col='black',alpha=0.5) +
  facet_wrap(~site,scales='free',nrow=4) +
  ylab("Ratio of seedlings to seeds per plot") +
  theme_bw()
```
\fi