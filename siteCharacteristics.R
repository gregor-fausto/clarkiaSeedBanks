siteVariables <- read.csv("~/Dropbox/projects/clarkiaScripts/data/reshapeData/siteAbiotic.csv",header=TRUE)
head(siteVariables)

plot(siteVariables$easting,siteVariables$elevation)
plot(siteVariables$easting,siteVariables$raw.aspect)

p = ggplot(siteVariables, aes(x=raw.aspect, y=elevation, fill=elevation)) +
  geom_bar(binwidth=1, stat='identity') +theme_light() +
  scale_fill_gradient(low='red', high='white', limits=c(0,100)) +
  theme(axis.title.y=element_text(angle=0))
p + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1))
p + coord_polar()

siteVariables

df=data.frame(site=siteVariables$site,group=c("upper","upper","mid","north","mid","mid","lake","upper","lake","lake","lower","lower","lower","mid","upper","upper","lake","north","lake","lower"))
df2 = data.frame(group = c("lower","mid","upper","lake","north"),ord = c(1,2,3,4,5))
siteVariables = siteVariables %>% dplyr::left_join(df,by='site') %>% 
  dplyr::left_join(df2,by='group') %>%
  dplyr::mutate(group = fct_reorder(group, ord))
  
ggplot(siteVariables, aes(x = raw.aspect)) +
  geom_histogram(binwidth = 15, boundary = -7.5) +
  coord_polar() +
  scale_x_continuous(limits = c(0,360)) +
  theme_bw() +
  facet_wrap(~group,nrow=1)

ggplot(siteVariables, aes(x = raw.aspect)) +
  geom_histogram(binwidth = 15, boundary = -7.5) +
  coord_polar() +
  scale_x_continuous(limits = c(0,360)) +
  theme_bw()

plot(siteVariables$average_slope,siteVariables$linear_azimuth)


ggplot(siteVariables, aes(x = raw.aspect,y=average_slope)) +
  geom_point() +
  coord_polar() +
  scale_x_continuous(limits = c(0,360)) +
  theme_bw() +
  facet_wrap(~group,nrow=1)
