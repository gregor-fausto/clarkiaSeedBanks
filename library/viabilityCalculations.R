library(tidyverse)
library(ggrepel)

# load data from the viability experiments
data <- readRDS("~/Dropbox/dataLibrary/workflow/data/viabilityExperiment.rds")

# parse data from the viability experiments
data <- data %>%
  dplyr::mutate(g = germCount/germStart, v = viabStain/viabStart) %>%
  dplyr::mutate(nu = g + v*(1-g)) %>%
  dplyr::group_by(siteViab,yearViab,ageViab) %>%
  dplyr::summarise(nu = mean(nu,na.rm=TRUE))

# for each site and year, calculate the v1 and v2

# plot the v2/v1 ratio
data_longer = data %>%
  tidyr::pivot_wider(names_from = ageViab, values_from = nu) %>%
  dplyr::mutate(ratio = `2`/`1`)

data_longer %>%
  dplyr::mutate(`2` = ifelse(is.na(`2`, mean())))

ggplot(data_longer,aes(x=`1`,y=`2`,color=yearViab,label=siteViab)) +
  geom_abline(intercept=0,slope=1,alpha=.5) +
  geom_point() +
  geom_text_repel(size=3,color="black") +
  theme_bw() + #xlim(c(0,1)) + ylim(c(0,1)) +
  xlab("Viability at end of year 1") +
  ylab("Viability at end of year 2")
