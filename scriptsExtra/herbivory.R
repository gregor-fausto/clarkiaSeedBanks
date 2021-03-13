# damaged undamaged fruits

data = readRDS("~/Dropbox/dataLibrary/postProcessingData/countSeedPerFruit.RDS")

data = data %>%
  dplyr::filter(year>2012) %>%
  dplyr::filter(demography==1)

data.sum = data %>%
  dplyr::group_by(site,year,damaged) %>%
  dplyr::summarise(n=n()) 

ggplot(data.sum) +geom_histogram(aes(n,fill=as.factor(damaged)),position='identity',alpha=.5)

library(tidyverse)

ggplot(data=data ) +
  geom_histogram(aes(sdno,fill=as.factor(damaged)),position='identity',alpha=.5) +
  facet_wrap(~year) +
  theme_bw()

delta = data %>%
  dplyr::group_by(site,year,damaged) %>%
  dplyr::summarise(mu = mean(sdno)) %>%
  tidyr::pivot_wider(names_from = damaged, values_from = mu) %>%
  dplyr::mutate(delta = `1`/`0`)


ggplot(data=delta ) +
  geom_histogram(aes(delta),position='identity',alpha=.5) +
  facet_wrap(~year) +
  geom_vline(xintercept=0) +
  theme_bw()

ggplot(data=delta ) +
  geom_point(aes(x=year,y=delta,color=site)) +
  geom_line(aes(x=year,y=delta,color=site)) +
 # facet_wrap(~site) +
  theme_bw()

ggplot(data=delta ) +
  geom_point(aes(x=year,y=delta,color=site)) +
  geom_line(aes(x=year,y=delta,color=site)) +
   facet_wrap(~site) +
  theme_bw()

ggplot(data=data %>% dplyr::filter(year==2013) ) +
  geom_histogram(aes(sdno,fill=as.factor(damaged)),position='identity',alpha=.5) +
  facet_wrap(~site,scales='free') +
  theme_bw()

ggplot(data=data %>% dplyr::filter(year==2014) ) +
  geom_histogram(aes(sdno,fill=as.factor(damaged)),position='identity',alpha=.5) +
  facet_wrap(~site,scales='free') +
  theme_bw()
