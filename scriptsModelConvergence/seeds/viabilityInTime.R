library(tidyverse)

viabFinal <- readRDS("~/Dropbox/clarkiaSeedBanks/scriptsModelChecking/seeds/viabilityFinal.RDS")

out<-viabFinal %>% dplyr::group_by(site,round,age) %>%
  dplyr::summarise(g = sum(germCount)/sum(germStart)) %>%
  dplyr::mutate(year=round+age+2005)

ggplot(out,aes(x=age,y=g,group=round)) +
  geom_line() +
  geom_point(aes(color=as.factor(year))) +
  facet_wrap(~site)

out<-viabFinal %>% dplyr::group_by(site,round,age) %>%
  dplyr::summarise(g = sum(viabStain,na.rm=TRUE)/sum(viabStart,na.rm=TRUE)) %>%
  dplyr::mutate(year=round+age+2005)

ggplot(out,aes(x=age,y=g,group=round)) +
  geom_line() +
  geom_point(aes(color=as.factor(year))) +
  facet_wrap(~site)

