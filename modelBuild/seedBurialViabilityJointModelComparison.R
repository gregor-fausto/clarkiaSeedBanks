# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Model comparison
# -------------------------------------------------------------------
# -------------------------------------------------------------------

x<-list.files("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/")

x<-x[sapply(x,stringr::str_detect,pattern="Fit")]
x<-x[sapply(x,stringr::str_detect,pattern="Pool")]
x<-x[sapply(x,stringr::str_detect,pattern="seedBags")]

x<-paste0("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/",x)

lapply(x,load,.GlobalEnv)

library(MCMCvis)

load("/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/seedBagExperimentData.rds")


N = dim(seedBagExperiment)[1]
data <- seedBagExperiment

focalParameter="pi"

ss_pool = MCMCchains(zc_pool,params=focalParameter)
ss_nopool = MCMCchains(zc_nopool,params=focalParameter)
ss_partialpool = MCMCchains(zc_partialpool,params=focalParameter)
ss_partialpoollogit = MCMCchains(zc_partialpoollogit,params=focalParameter)
ss_partialpoolhyperpriors = MCMCchains(zc_partialpoolhyperpriors,params=focalParameter)

rm(zc_pool,zc_nopool,zc_partialpool,zc_partialpoollogit,zc_partialpoolhyperpriors)

#
p_pool <- t(apply(ss_pool,2,function(x) quantile(x,probs=c(.1,.5,.9)))) %>%
  data.frame %>% 
  dplyr::rename(theta_50 = X50.,theta_10 = X10., theta_90 = X90.)

#
p_nopool <- t(apply(ss_nopool,2,function(x) quantile(x,probs=c(.1,.5,.9)))) %>%
  data.frame %>% 
  dplyr::rename(theta_50 = X50.,theta_10 = X10., theta_90 = X90.)

#
p_partialpool <- t(apply(ss_partialpool,2,function(x) quantile(x,probs=c(.1,.5,.9)))) %>%
  data.frame %>% 
  dplyr::rename(theta_50 = X50.,theta_10 = X10., theta_90 = X90.)

# p_partialpool <- p_partialpool %>% 
#   bind_cols(bagNoUnique = 1:dim(p_partialpool)[1])
# 
# p_partialpool<- data %>%
#   dplyr::bind_cols(bagNoUnique = data$bagNo) %>%
#   dplyr::left_join(p_partialpool,by="bagNoUnique") 

#
p_partialpoollogit <- t(apply(ss_partialpoollogit,2,function(x) quantile(x,probs=c(.1,.5,.9)))) %>%
  data.frame %>% 
  dplyr::rename(theta_50 = X50.,theta_10 = X10., theta_90 = X90.)

# p_partialpoollogit <- p_partialpoollogit %>% 
#   bind_cols(bagNoUnique = 1:dim(p_partialpoollogit)[1])
# 
# p_partialpoollogit<- dat %>%
#   dplyr::bind_cols(bagNoUnique = data$bagNo) %>%
#   dplyr::left_join(p_partialpoollogit,by="bagNoUnique") 
#
p_partialpoolhyperpriors <- t(apply(ss_partialpoolhyperpriors,2,function(x) quantile(x,probs=c(.1,.5,.9)))) %>%
  data.frame %>% 
  dplyr::rename(theta_50 = X50.,theta_10 = X10., theta_90 = X90.)

# p_partialpoolhyperpriors <- p_partialpoolhyperpriors %>% 
#   bind_cols(bagNoUnique = 1:dim(p_partialpoolhyperpriors)[1])
# 
# p_partialpoolhyperpriors<- dat %>%
#   dplyr::bind_cols(bagNoUnique = data$bagNo) %>%
#   dplyr::left_join(p_partialpoolhyperpriors,by="bagNoUnique")

lapply(list(p_pool,p_nopool,p_partialpool,p_partialpoollogit,p_partialpoolhyperpriors),dim)

plot(p_partialpool$theta_50,p_partialpoollogit$theta_50);abline(a=0,b=1)
plot(p_partialpool$theta_50,p_partialpoolhyperpriors$theta_50);abline(a=0,b=1)

df_plot2<-data.frame(group = rep(data$bagNo, 5),
                     model = c(rep("complete pooling", N),
                               rep("no pooling", N),
                               rep("partial pooling", N),
                               rep('partial pooling, logit', N),
                               rep("partial pooling, hyperpriors", N)),
                     y = c(rep(p_pool$theta_50,N),
                           p_nopool$theta_50,
                           p_partialpool$theta_50,
                           p_partialpoollogit$theta_50,
                           p_partialpoolhyperpriors$theta_50))

dfPlot <- df_plot2 %>% 
  group_by(model) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = model, values_from = y) %>% 
  select(-row) %>%
  bind_cols(site=dat$site)

library(GGally)
pairsPlot <- ggpairs(dfPlot, columns = 3:6, 
                     ggplot2::aes(alpha = 0.4,color=as.factor(site)),
                     upper = list(continuous = wrap("cor", size = 2)))
pairsPlot

dfPlot %>% View

ggsave(plot=pairsPlot,
       filename="/Users/Gregor/Dropbox/dataLibrary/clarkiaSeedBanks/viabilityPairsPlot.png",
       width=12,height=12)

# ggplot() + 
#   geom_point(data = dfPlot,aes(x = `no pooling`,y = `partial pooling`)) + 
#   geom_errorbar(aes(x= dfPlot$`no pooling`, ymin = dfPlotLo$`partial pooling`,ymax = dfPlotHi$`partial pooling`)) +
#   geom_errorbarh(aes(y= dfPlot$`partial pooling`, xmin = dfPlotLo$`no pooling`,xmax = dfPlotHi$`no pooling`))
# 
# ggplot() + 
#   geom_point(data = dfPlot,aes(x = `partial pooling`,y = `partial pooling, logit`)) + 
#   geom_errorbar(aes(x= dfPlot$`partial pooling`, ymin = dfPlotLo$`partial pooling, logit`,ymax = dfPlotHi$`partial pooling, logit`)) +
#   geom_errorbarh(aes(y= dfPlot$`partial pooling, logit`, xmin = dfPlotLo$`partial pooling`,xmax = dfPlotHi$`partial pooling`))
# 
# ggplot() + 
#   geom_point(data = dfPlot,aes(x = `partial pooling`,y = `partial pooling, hyperpriors`)) + 
#   geom_errorbar(aes(x= dfPlot$`partial pooling`, ymin = dfPlotLo$`partial pooling, hyperpriors`,ymax = dfPlotHi$`partial pooling, hyperpriors`)) +
#   geom_errorbarh(aes(y= dfPlot$`partial pooling, hyperpriors`, xmin = dfPlotLo$`partial pooling`,xmax = dfPlotHi$`partial pooling`))

