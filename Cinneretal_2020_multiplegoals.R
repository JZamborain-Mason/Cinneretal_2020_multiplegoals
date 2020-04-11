#Code for Cinner et al.2020: Meeting fisheries, ecosystem function, and biodiversity goals in a human dominated world
#ARC Centre of Excellence for Coral Reef Studies
#R version 3.4.2 (2017-09-28)


##Remove everything from the environment
rm(list = ls())


##set working directory
#NOTE:Please add your working directory (i.e., location where data is stored)
setwd("") 

##load required libraries
library(ggplot2) ##H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016
library(MuMIn) #Kamil Barton (2019) #. MuMIn: Multi-Model Inference. R package version 1.43.6.
library(mice) #Stef van Buuren, Karin Groothuis-Oudshoorn (2011). mice: Multivariate Imputation by Chained Equations in.R. Journal of Statistical Software, 45(3), 1-67. URL https://www.jstatsoft.org/v45/i03/.
library(ggpubr) #Alboukadel Kassambara (2018). ggpubr: 'ggplot2' Based Publication Ready Plots. 
library(rstan) #Stan Development Team (2018). RStan: the R interface to Stan.
library(rstanarm) #Goodrich B, Gabry J, Ali I & Brilleman S. (2018). rstanarm: Bayesian applied regression modeling via
library(bayesplot) #Jonah Gabry and Tristan Mahr (2018). bayesplot: Plotting for Bayesian Models. 
library(tidyverse) #Hadley Wickham (2017). tidyverse: Easily Install and Load the 'Tidyverse'. 
library(broom) #David Robinson and Alex Hayes (2019). broom: Convert Statistical Analysis Objects into Tidy Tibbles.
library(coda) #Martyn Plummer, Nicky Best, Kate Cowles and Karen Vines (2006). CODA: Convergence Diagnosis and Output.Analysis for MCMC, R News, vol 6, 7-11
library(data.table) #Matt Dowle and Arun Srinivasan (2019). data.table: Extension of `data.frame`
library(lme4) #Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015) #. Fitting Linear Mixed-Effects Models.Using lme4. Journal of Statistical Software, 67(1), 1-48. doi:10.18637/jss.v067.i01.
library(brms) # Paul-Christian Bürkner (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. Journal of StatisticalSoftware, 80(1), 1-28. doi:10.18637/jss.v080.i01
library(DHARMa)#Florian Hartig (2019). DHARMa: Residual Diagnostics for Hierarchical (Multi-Level / Mixed) Regression Models. R package version 0.2.4. https://CRAN.R-project.org/package=DHARMa
library(plyr)#Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical Software, 40(1), 1-29. URL http://www.jstatsoft.org/v40/i01/.
library(networkD3) # J.J. Allaire, Christopher Gandrud, Kenton Russell and CJ Yetman (2017). networkD3: D3 JavaScript Network Graphs from R. R package version 0.4. https://CRAN.R-project.org/package=networkD3


#Functions
#to draw histogram on pairs plot
panel_hist = function(x, ...)
{
    usr = par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h = hist(x, plot = FALSE)
    breaks = h$breaks; nB = length(breaks)
    y = h$counts; y = y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}
#Standardize function for continuous variables
standardize = function(x){(x-mean(x, na.rm=T))/(2*sd(x, na.rm=T))} 

#variance inflation factors  function for mixed models
vif.mer = function (fit) {
  ## adapted from rms::vif
  v <-vcov(fit)
  nam = names(fixef(fit))
  ## exclude intercepts
  ns = sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v = v[-(1:ns), -(1:ns), drop = FALSE]
    nam = nam[-(1:ns)]
  }
  d = diag(v)^0.5
  v = diag(solve(v/(d %o% d)))
  names(v) = nam
  v
}

#correlation function
panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr = par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = abs(cor(x, y, method = "pearson",use = "complete.obs"))
  txt = format(c(r, 0.123456789), digits=digits)[1]
  txt = paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor = 0.9/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r*2)
}
#function to extract the probabilities from the prediced distribution
quantInv = function(distr, value) ecdf(distr)(value)

#upload data
alldata=read.csv("Cinneretal2020_multiplegoals_data.csv",header = T)

#mpa size and age restrictions
alldata$mpacondition=ifelse(alldata$Protection=="UnfishedHigh" &(alldata$MPAage>4|alldata$MPAage==4) & (alldata$NTZarea>2 |alldata$NTZarea==0),1,0)
alldata$keep=ifelse(!alldata$Protection=="UnfishedHigh",1,ifelse(alldata$Protection=="UnfishedHigh" &alldata$mpacondition==1,1,0))
#alldata=alldata[alldata$keep==1,] #Uncomment if you want to use more strick guidelines for the marine reserves used


##Explanatory variables: relevel categorical variables and standardize continuous variables

#reef-scale covariates
alldata$DepthCategory=relevel(alldata$DepthCategory,ref="4-10m")
alldata$CleanHabitat=relevel(alldata$CleanHabitat,ref="Slope")
alldata$Protection=relevel(alldata$Protection,ref="Fished")
alldata$CensusMethod=relevel(alldata$CensusMethod,ref="Standard belt transect")
alldata$sTotal_sampling_area=standardize(log(alldata$Total_sampling_area))
#gravities with different exponents
alldata$sgrav_tot1=standardize(log(alldata$gravtot5001+min(alldata$gravtot5001[alldata$gravtot5001>0])))
alldata$sgrav_tot2=standardize(log(alldata$gravtot5002+min(alldata$gravtot5002[alldata$gravtot5002>0])))
alldata$sgrav_tot3=standardize(log(alldata$gravtot5003+min(alldata$gravtot5003[alldata$gravtot5003>0])))

#reef cluster-scale covariates
alldata$sOcean_prod=standardize(log(alldata$Ocean_prod))
alldata$sClimate_stress=standardize(alldata$Climate_stress)
alldata$sRegional_population_growth=standardize(alldata$Regional_population_growth)

#nation/state-scale covariates
alldata$sReef_fish_landings_per_km2=standardize(log(alldata$Reef_fish_landings_per_km2+1))
alldata$sLarger_pop_size=standardize(log(alldata$Larger_pop_size+1))
alldata$sHDI=standardize(alldata$HDI)

##check that multicolinearity is not  a concern
#pairs correlations among covariates

#smaller scale
pairs(~sgrav_tot2+
        Protection+
        sRegional_population_growth+
        sOcean_prod+
        sClimate_stress+
        DepthCategory+
        CleanHabitat+
        sTotal_sampling_area+
        CensusMethod, data= alldata,lower.panel=panel.cor )
#larger scale
pairs(~sHDI+
  sLarger_pop_size+
  sReef_fish_landings_per_km2, data= alldata,lower.panel=panel.cor )

#variance inflation factor
VIF.table=as.data.frame(vif.mer(lmer(log(Biomass_above20cm+1)~
                                       DepthCategory+ CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot2+
                                       sRegional_population_growth+sOcean_prod+sClimate_stress+
                                       sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+(1|Larger/ReefCluster), data=alldata)))
colnames(VIF.table)="VIF"
print(VIF.table)

##exploring and transforming response variables 

#Parrotfish scraping: continuous variable with almost 31% of the data containing 0's
#hurdle model: two part model. Part 1 estimates the probability of observing the function (binomial), and part 2 estimates the standardized effect size of the covaraites and the mean herbivore function given that we observed the function (normally distributed when log transformed: gaussian)
length(alldata$Scraping_potential[alldata$Scraping_potential==0])/length(alldata$Scraping_potential[!is.na(alldata$Scraping_potential)])
hist(log(alldata$Scraping_potential+1))

#Trait diversity: continuous variable with no zeros. 
hist(log(alldata$Trait_diversity))
alldata$tTrait_diversity=log(alldata$Trait_diversity)

#Reef fish biomass> 20cm: only 3 % of the data are 0's. Consequently, we log+1 trasform it so it has a normal distribution
length(alldata$Biomass_above20cm[alldata$Biomass_above20cm==0])/length(alldata$Biomass_above20cm)
hist(log(alldata$Biomass_above20cm+1))
alldata$tBiomass_above20cm=log(alldata$Biomass_above20cm+1)

#correlations among response variables
pairs(~tBiomass_above20cm+tTrait_diversity+log(Scraping_potential+1),data=alldata,lower.panel=panel.cor, 
      pch = 21, bg = "darkgrey",labels=c("log(Biomass >20cm+1)","log(Trait diversity)", "log(Parrotfish scraping+1)"),cex.labels=1.5,font.labels=2,diag.panel =panel_hist, hist.col="grey")

#correlation of non-zero data
nonzero=alldata[alldata$Scraping_potential>0 & alldata$Biomass_above20cm>0,]
pairs(~tBiomass_above20cm+tTrait_diversity+log(Scraping_potential),data=nonzero,lower.panel=panel.cor, 
      pch = 21, bg = "darkgrey",labels=c("log(Biomass >20cm+1)","log(Trait diversity)", "log(Parrotfish scraping)"),cex.labels=1.5,font.labels=2,diag.panel =panel_hist, hist.col="grey",main = "Non-zero data")
#correlations among explanatory variables
pairs(~DepthCategory+Geographic_Basin+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot2+
        sRegional_population_growth+sOcean_prod+sClimate_stress+
        sHDI+sLarger_pop_size+sReef_fish_landings_per_km2,data=alldata,lower.panel=panel.cor, 
      pch = 21, bg = "darkgrey",labels=c("Depth Category","Basin","Habitat type", "Management
protection", "Sampling
method","Sampling area","Total gravity","Local
population
growth","Ocean
productivity", "Climate stress
index","HDI","Population size",  "Reef fish
landings"),cex.labels=0.9,font.labels=2)

#countries included
countries=as.data.frame(unique(alldata$Larger))
#write.csv(countries, "countries.included.csv",row.names=F)


##Models###########################################################################################################################
#For each response variable we run 4 models:
#(i) null;(ii) using gravity exponent1;(iii)using gravity exponent2 ; and (iv) using gravity exponent3
#For the paper we use gravity exponent 2

#Parrotfish scraping#################################################
her_model_NULL=brm(Scraping_potential~
                      (1|Larger/ReefCluster),
                    data=alldata[!is.na(alldata$Scraping_potential),],family=hurdle_lognormal(link = "identity", link_sigma = "log",
                    link_hu = "logit"),
                    iter=20000,  warmup=19000,
                    chains=4) 

her_model_totexp1=brm(Scraping_potential~
                                DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot1+
                                sRegional_population_growth+sOcean_prod+sClimate_stress+
                                sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                                (1|Larger/ReefCluster),
                              data=alldata[!is.na(alldata$Scraping_potential),],family=hurdle_lognormal(link = "identity", link_sigma = "log",
                              link_hu = "logit"),
                              iter=10000,  warmup=9000,
                              chains=4) 
her_model_totexp2=brm(Scraping_potential~
                                DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot2+
                                sRegional_population_growth+sOcean_prod+sClimate_stress+
                                sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                                (1|Larger/ReefCluster),
                              data=alldata[!is.na(alldata$Scraping_potential),],family=hurdle_lognormal(link = "identity", link_sigma = "log",
                              link_hu = "logit"),
                              iter=10000,  warmup=9000,
                              chains=4) 
her_model_totexp3=brm(Scraping_potential~
                                DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot3+
                                sRegional_population_growth+sOcean_prod+sClimate_stress+
                                sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                                (1|Larger/ReefCluster),
                              data=alldata[!is.na(alldata$Scraping_potential),],family=hurdle_lognormal(link = "identity", link_sigma = "log",
                              link_hu = "logit"),
                              iter=10000,  warmup=9000,
                              chains=4) 
model_comparison2=loo(her_model_NULL,her_model_totexp1,her_model_totexp2,her_model_totexp3)
model_comparison2=as.data.frame(model_comparison2$diffs)

#check diagnostics and model fit of used model
pp_check(her_model_totexp2,nsamples=4000)+xlim(c(0,300))
yrep_hf=as.data.frame(posterior_predict(her_model_totexp2))
x = createDHARMa(simulatedResponse=t(yrep_hf), observedResponse=alldata$Scraping_potential[!is.na(alldata$Scraping_potential)])
windows()
plotResiduals(x)

#model posterior
posteriorestimates_hf=as.data.frame(posterior_summary(her_model_totexp2))

#Trait diversity##################################################

#correlation with previous FDis
fundiv_model_NULL=brm(tTrait_diversity~
                                (1|Larger/ReefCluster),
                              data=alldata[!is.na(alldata$Trait_diversity),],family=gaussian,
                              iter=10000,  warmup=9000,
                              chains=4, 
                              thin=10) 
fundiv_model_totexp1=brm(tTrait_diversity~
                                   DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot1+
                                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                                   (1|Larger/ReefCluster),
                                 data=alldata[!is.na(alldata$Trait_diversity),],family=gaussian,
                                 iter=10000,  warmup=9000,
                                 chains=4) 
fundiv_model_totexp2=brm(tTrait_diversity~
                                   DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot2+
                                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                                   (1|Larger/ReefCluster),
                                 data=alldata[!is.na(alldata$Trait_diversity),],family=gaussian,
                                 iter=10000,  warmup=9000,
                                 chains=4) 
fundiv_model_totexp3=brm(tTrait_diversity~
                                   DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot3+
                                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                                   (1|Larger/ReefCluster),
                                 data=alldata[!is.na(alldata$Trait_diversity),],family=gaussian,
                                 iter=10000,  warmup=9000,
                                 chains=4) 
#model comparison
ftot2=loo(fundiv_model_NULL,fundiv_model_totexp1,fundiv_model_totexp2,fundiv_model_totexp3)
model_comparison2=rbind(model_comparison2,ftot2$diffs)

#check model fit
pp_check(fundiv_model_totexp2, nsamples=4000)
hist(resid(fundiv_model_totexp2)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(fundiv_model_totexp2), x=fitted(fundiv_model_totexp2)))

#model posterior
posteriorestimates_fd=as.data.frame(posterior_summary(fundiv_model_totexp2))

#Biomass of reef fish >20cm ##################################################
B20_model_NULL=brm(tBiomass_above20cm~
                             (1|Larger/ReefCluster),
                           data=alldata,family=gaussian,
                           iter=10000,  warmup=9000,
                           chains=4, 
                           thin=10) 

B20_model_totexp1=brm(tBiomass_above20cm~
                                DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot1+
                                sRegional_population_growth+sOcean_prod+sClimate_stress+
                                sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                                (1|Larger/ReefCluster),
                              data=alldata,family=gaussian,
                              iter=10000,  warmup=9000,
                              chains=4) 
B20_model_totexp2=brm(tBiomass_above20cm~
                                DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot2+
                                sRegional_population_growth+sOcean_prod+sClimate_stress+
                                sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                                (1|Larger/ReefCluster),
                              data=alldata,family=gaussian,
                              iter=10000,  warmup=9000,
                              chains=4) 
B20_model_totexp3=brm(tBiomass_above20cm~
                                DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot3+
                                sRegional_population_growth+sOcean_prod+sClimate_stress+
                                sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                                (1|Larger/ReefCluster),
                              data=alldata,family=gaussian,
                              iter=10000,  warmup=9000,
                              chains=4) 
#model comparison
btot2=loo(B20_model_totexp1,B20_model_totexp2,B20_model_totexp3)
model_comparison2=rbind(model_comparison2,btot2$diffs)
#write.csv(model_comparison2, "totgrav_sensitiity_comparison.csv", row.names=T)

#check model fit
pp_check(B20_model_totexp2,nsamples=NULL)
hist(resid(B20_model_totexp2)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(B20_model_totexp2), x=fitted(B20_model_totexp2)))

#model posterior
posteriorestimates_b20=as.data.frame(posterior_summary(B20_model_totexp2))


#coefficient plots
#Parrotfish scraping
coefplot_hf=posteriorestimates_hf[2:18,]
coefplot_hf$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Total gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_hf$sign=ifelse(coefplot_hf$Q2.5<0 & coefplot_hf$Q97.5<0, "negative",ifelse(coefplot_hf$Q2.5>0 & coefplot_hf$Q97.5>0, "positive", "no effect"))
coefplot_hf$strength=ifelse(coefplot_hf$sign=="no effect", "open", "closed")
coefplot_hf$order=c(4,5,6,7,8,11,12,3,2,1,17,16,10,9,13,14,15)
coefplot_hf[order(coefplot_hf$order),]
coefplot_hf$variable = factor(coefplot_hf$variable, levels = coefplot_hf$variable[order(coefplot_hf$order)])

coefplot_ps=
  ggplot(coefplot_hf,aes(x=variable,y=Estimate,ymin=Q2.5,ymax=Q97.5))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  scale_y_continuous("Standardized effect size",limits=c(-3,3.6))+xlab("")+ggtitle("Parrotfish scraping")+
  theme(axis.text.y=element_blank())


#BIOMASS REEF FISH >20 CM
coefplot_B20=posteriorestimates_b20[2:18,]
coefplot_B20$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Total gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_B20$sign=ifelse(coefplot_B20$Q2.5<0 & coefplot_B20$Q97.5<0, "negative",ifelse(coefplot_B20$Q2.5>0 & coefplot_B20$Q97.5>0, "positive", "no effect"))
coefplot_B20$strength=ifelse(coefplot_B20$sign=="no effect", "open", "closed")
coefplot_B20$order=c(4,5,6,7,8,11,12,3,2,1,17,16,10,9,13,14,15)
coefplot_B20[order(coefplot_B20$order),]
coefplot_B20$variable = factor(coefplot_B20$variable, levels = coefplot_B20$variable[order(coefplot_B20$order)])

coefplot_b20=
  ggplot(coefplot_B20,aes(x=variable,y=Estimate,ymin=Q2.5,ymax=Q97.5))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  scale_y_continuous("",limits=c(-3,3.6))+xlab("")+ggtitle("Biomass > 20cm")+theme_classic()


#Trait diversity
coefplot_fd=posteriorestimates_fd[2:18,]
coefplot_fd$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Total gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_fd$sign=ifelse(coefplot_fd$Q2.5<0 & coefplot_fd$Q97.5<0, "negative",ifelse(coefplot_fd$Q2.5>0 & coefplot_fd$Q97.5>0, "positive", "no effect"))
coefplot_fd$strength=ifelse(coefplot_fd$sign=="no effect", "open", "closed")
coefplot_fd$order=c(4,5,6,7,8,11,12,3,2,1,17,16,10,9,13,14,15)
coefplot_fd[order(coefplot_fd$order),]
coefplot_fd$variable = factor(coefplot_fd$variable, levels = coefplot_fd$variable[order(coefplot_fd$order)])

coefplot_ts=
  ggplot(coefplot_fd,aes(x=variable,y=Estimate,ymin=Q2.5,ymax=Q97.5))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  scale_y_continuous("",limits=c(-3,3.2))+xlab("")+ggtitle("Trait diversity")+
  theme(axis.text.y=element_blank())

#combined coefplot figure
ggarrange(coefplot_b20,coefplot_ps,coefplot_ts,nrow=1,ncol=3, widths=c(1.9,1,1))

##################################################################################

# probabilities of passing different thresholds

#first we marginalize the response variables (i.e., take the effect of the things we dont want off)
#we are marginalizing 3 different ways
#1: taking the effect of only sampling area and sampling methodology as if everyone had been consistent sampling (for plotting the maps)
#2: taking the effect of sampling and habitat/depth sampled so everything is consistent (to estimate the number of sites passing different thresholds given that they were all slopes and at 4-10 m depth)
#3: taking off the effect of the above and also the estimated random effect (i.e., at average global conditions) (to estimate the benchmarks and thresholds)

alldata$marg_B20=(alldata$Biomass_above20cm)/exp(fixef(B20_model_totexp2)["sTotal_sampling_area",1]*alldata$sTotal_sampling_area+
                                                   ifelse(alldata$CensusMethod=="Distance sampling",fixef(B20_model_totexp2)["CensusMethodDistancesampling",1],
                                                          ifelse(alldata$CensusMethod=="Point intercept",fixef(B20_model_totexp2)["CensusMethodPointintercept",1],0))+
                                                   ifelse(alldata$DepthCategory==">10m",fixef(B20_model_totexp2)["DepthCategory>10m",1],
                                                          ifelse(alldata$DepthCategory=="0-4m",fixef(B20_model_totexp2)["DepthCategory0M4m",1],0))+
                                                   ifelse(alldata$CleanHabitat=="Crest",fixef(B20_model_totexp2)["CleanHabitatCrest",1],
                                                          ifelse(alldata$CleanHabitat=="Flat",fixef(B20_model_totexp2)["CleanHabitatFlat",1],
                                                                 ifelse(alldata$CleanHabitat=="Lagoon_Back reef",fixef(B20_model_totexp2)["CleanHabitatLagoon_Backreef",1],0))))

alldata$marg_B20_habitat=(alldata$Biomass_above20cm)/exp(fixef(B20_model_totexp2)["sTotal_sampling_area",1]*alldata$sTotal_sampling_area+
                                                           ifelse(alldata$CensusMethod=="Distance sampling",fixef(B20_model_totexp2)["CensusMethodDistancesampling",1],
                                                                  ifelse(alldata$CensusMethod=="Point intercept",fixef(B20_model_totexp2)["CensusMethodPointintercept",1],0)))


#match random effects with data
ranefB20_largerreefcluster=as.data.frame(ranef(B20_model_totexp2)$`Larger:ReefCluster`[,1,1])
colnames(ranefB20_largerreefcluster)="ranef_larger_reefcluster_B20"
ranefB20_largerreefcluster$Larger_ReefCluster=rownames(ranefB20_largerreefcluster)
ranefB20_larger=as.data.frame(ranef(B20_model_totexp2)$Larger[,1,1])
colnames(ranefB20_larger)="ranef_larger_B20"
ranefB20_larger$Larger=rownames(ranefB20_larger)
alldata$Larger_ReefCluster=paste(alldata$Larger, alldata$ReefCluster, sep="_")
alldata=merge(alldata, ranefB20_largerreefcluster, by="Larger_ReefCluster",all.x=T)
alldata=merge(alldata, ranefB20_larger, by="Larger",all.x=T)


alldata$marg_B20_ranef=(alldata$Biomass_above20cm)/exp(alldata$ranef_larger_reefcluster_B20+alldata$ranef_larger_B20+
                                                         fixef(B20_model_totexp2)["sTotal_sampling_area",1]*alldata$sTotal_sampling_area+
                                                         ifelse(alldata$CensusMethod=="Distance sampling",fixef(B20_model_totexp2)["CensusMethodDistancesampling",1],
                                                                ifelse(alldata$CensusMethod=="Point intercept",fixef(B20_model_totexp2)["CensusMethodPointintercept",1],0))+
                                                         ifelse(alldata$DepthCategory==">10m",fixef(B20_model_totexp2)["DepthCategory>10m",1],
                                                                ifelse(alldata$DepthCategory=="0-4m",fixef(B20_model_totexp2)["DepthCategory0M4m",1],0))+
                                                         ifelse(alldata$CleanHabitat=="Crest",fixef(B20_model_totexp2)["CleanHabitatCrest",1],
                                                                ifelse(alldata$CleanHabitat=="Flat",fixef(B20_model_totexp2)["CleanHabitatFlat",1],
                                                                       ifelse(alldata$CleanHabitat=="Lagoon_Back reef",fixef(B20_model_totexp2)["CleanHabitatLagoon_Backreef",1],0))))


raneffundiv_largerreefcluster=as.data.frame(ranef(fundiv_model_totexp2)$`Larger:ReefCluster`[,1,1])
colnames(raneffundiv_largerreefcluster)="ranef_larger_reefcluster_td"
raneffundiv_largerreefcluster$Larger_ReefCluster=rownames(raneffundiv_largerreefcluster)
raneffundiv_larger=as.data.frame(ranef(fundiv_model_totexp2)$Larger[,1,1])
colnames(raneffundiv_larger)="ranef_larger_td"
raneffundiv_larger$Larger=rownames(raneffundiv_larger)
alldata=merge(alldata, raneffundiv_largerreefcluster, by="Larger_ReefCluster",all.x=T)
alldata=merge(alldata, raneffundiv_larger, by="Larger",all.x=T)

alldata$marg_td_ranef=alldata$Trait_diversity/exp((alldata$ranef_larger_reefcluster_td+alldata$ranef_larger_td+
                                                 fixef(fundiv_model_totexp2)["sTotal_sampling_area",1]*alldata$sTotal_sampling_area+
                                                 ifelse(alldata$CensusMethod=="Distance sampling",fixef(fundiv_model_totexp2)["CensusMethodDistancesampling",1],
                                                        ifelse(alldata$CensusMethod=="Point intercept",fixef(fundiv_model_totexp2)["CensusMethodPointintercept",1],0))+
                                                 ifelse(alldata$DepthCategory==">10m",fixef(fundiv_model_totexp2)["DepthCategory>10m",1],
                                                        ifelse(alldata$DepthCategory=="0-4m",fixef(fundiv_model_totexp2)["DepthCategory0M4m",1],0))+
                                                 ifelse(alldata$CleanHabitat=="Crest",fixef(fundiv_model_totexp2)["CleanHabitatCrest",1],
                                                        ifelse(alldata$CleanHabitat=="Flat",fixef(fundiv_model_totexp2)["CleanHabitatFlat",1],
                                                               ifelse(alldata$CleanHabitat=="Lagoon_Back reef",fixef(fundiv_model_totexp2)["CleanHabitatLagoon_Backreef",1],0)))))

alldata$marg_td=alldata$Trait_diversity/exp((fixef(fundiv_model_totexp2)["sTotal_sampling_area",1]*alldata$sTotal_sampling_area+
                                           ifelse(alldata$CensusMethod=="Distance sampling",fixef(fundiv_model_totexp2)["CensusMethodDistancesampling",1],
                                                  ifelse(alldata$CensusMethod=="Point intercept",fixef(fundiv_model_totexp2)["CensusMethodPointintercept",1],0))+
                                           ifelse(alldata$DepthCategory==">10m",fixef(fundiv_model_totexp2)["DepthCategory>10m",1],
                                                  ifelse(alldata$DepthCategory=="0-4m",fixef(fundiv_model_totexp2)["DepthCategory0M4m",1],0))+
                                           ifelse(alldata$CleanHabitat=="Crest",fixef(fundiv_model_totexp2)["CleanHabitatCrest",1],
                                                  ifelse(alldata$CleanHabitat=="Flat",fixef(fundiv_model_totexp2)["CleanHabitatFlat",1],
                                                         ifelse(alldata$CleanHabitat=="Lagoon_Back reef",fixef(fundiv_model_totexp2)["CleanHabitatLagoon_Backreef",1],0)))))

alldata$marg_td_habitat=alldata$Trait_diversity/exp((fixef(fundiv_model_totexp2)["sTotal_sampling_area",1]*alldata$sTotal_sampling_area+
                                                   ifelse(alldata$CensusMethod=="Distance sampling",fixef(fundiv_model_totexp2)["CensusMethodDistancesampling",1],
                                                          ifelse(alldata$CensusMethod=="Point intercept",fixef(fundiv_model_totexp2)["CensusMethodPointintercept",1],0))))

ranefher_largerreefcluster=as.data.frame(ranef(her_model_totexp2)$`Larger:ReefCluster`[,1,1])
colnames(ranefher_largerreefcluster)="ranef_larger_reefcluster_ps"
ranefher_largerreefcluster$Larger_ReefCluster=rownames(ranefher_largerreefcluster)
ranefher_larger=as.data.frame(ranef(her_model_totexp2)$Larger[,1,1])
colnames(ranefher_larger)="ranef_larger_ps"
ranefher_larger$Larger=rownames(ranefher_larger)
alldata=merge(alldata, ranefher_largerreefcluster, by="Larger_ReefCluster",all.x=T)
alldata=merge(alldata, ranefher_larger, by="Larger",all.x=T)
alldata$marg_ps_ranef=alldata$Scraping_potential/exp(alldata$ranef_larger_reefcluster_ps+alldata$ranef_larger_ps+
                                                        fixef(her_model_totexp2)["sTotal_sampling_area",1]*alldata$sTotal_sampling_area+
                                                        ifelse(alldata$CensusMethod=="Distance sampling",fixef(her_model_totexp2)["CensusMethodDistancesampling",1],
                                                               ifelse(alldata$CensusMethod=="Point intercept",fixef(her_model_totexp2)["CensusMethodPointintercept",1],0))+
                                                        ifelse(alldata$DepthCategory==">10m",fixef(her_model_totexp2)["DepthCategory>10m",1],
                                                               ifelse(alldata$DepthCategory=="0-4m",fixef(her_model_totexp2)["DepthCategory0M4m",1],0))+
                                                        ifelse(alldata$CleanHabitat=="Crest",fixef(her_model_totexp2)["CleanHabitatCrest",1],
                                                               ifelse(alldata$CleanHabitat=="Flat",fixef(her_model_totexp2)["CleanHabitatFlat",1],
                                                                      ifelse(alldata$CleanHabitat=="Lagoon_Back reef",fixef(her_model_totexp2)["CleanHabitatLagoon_Backreef",1],0))))

alldata$marg_ps=alldata$Scraping_potential/exp(fixef(her_model_totexp2)["sTotal_sampling_area",1]*alldata$sTotal_sampling_area+
                                                  ifelse(alldata$CensusMethod=="Distance sampling",fixef(her_model_totexp2)["CensusMethodDistancesampling",1],
                                                         ifelse(alldata$CensusMethod=="Point intercept",fixef(her_model_totexp2)["CensusMethodPointintercept",1],0))+
                                                  ifelse(alldata$DepthCategory==">10m",fixef(her_model_totexp2)["DepthCategory>10m",1],
                                                         ifelse(alldata$DepthCategory=="0-4m",fixef(her_model_totexp2)["DepthCategory0M4m",1],0))+
                                                  ifelse(alldata$CleanHabitat=="Crest",fixef(her_model_totexp2)["CleanHabitatCrest",1],
                                                         ifelse(alldata$CleanHabitat=="Flat",fixef(her_model_totexp2)["CleanHabitatFlat",1],
                                                                ifelse(alldata$CleanHabitat=="Lagoon_Back reef",fixef(her_model_totexp2)["CleanHabitatLagoon_Backreef",1],0))))

alldata$marg_ps_habitat=alldata$Scraping_potential/exp(fixef(her_model_totexp2)["sTotal_sampling_area",1]*alldata$sTotal_sampling_area+
                                                          ifelse(alldata$CensusMethod=="Distance sampling",fixef(her_model_totexp2)["CensusMethodDistancesampling",1],
                                                                 ifelse(alldata$CensusMethod=="Point intercept",fixef(her_model_totexp2)["CensusMethodPointintercept",1],0)))

#Probabilities of passing different thresholds

#non-NA data
nonNAHFdata=alldata[!is.na(alldata$Scraping_potential),]
nonNAtd=alldata[!is.na(alldata$Trait_diversity),]

#Benchmarks: change p to 0.8 or 0.95 for different benckmarks used
rherb=exp(quantile(log(alldata$marg_ps_ranef+1),p=0.9, na.rm=T))-1 
rfunc=exp(quantile(log(alldata$marg_td_ranef),p=0.9, na.rm=T))
rpred=exp(quantile(log(alldata$marg_B20_ranef+1),p=0.9, na.rm=T))-1

#fitted values
biomass_fitted=fitted(B20_model_totexp2)[,1]
alldata$pred_biomass=exp(biomass_fitted-1)
herb_fitted=fitted(her_model_totexp2)[,1]
nonNAHFdata$pred_herb=herb_fitted
func_fitted=fitted(fundiv_model_totexp2)[,1]
nonNAtd$pred_fundiv=exp(func_fitted)

#jitter lat and long points
alldata$Site_Lat2=alldata$Site_Lat+runif(length(alldata$Site_Lat), min=0, max=4)
alldata$Site_Long2=alldata$Site_Long+runif(length(alldata$Site_Long), min=0, max=3)
alldata$Site_Lat2=ifelse(alldata$Site_Lat2>23.5, alldata$Site_Lat,alldata$Site_Lat2)
alldata$lon2 <- ifelse(alldata$Site_Long2 < -25, alldata$Site_Long2 + 360, alldata$Site_Long2) 


#merge the predicted data
data_predvalues=alldata[,c("UniqueSite","pred_biomass","Protection","Scraping_potential","Site_Lat","Site_Long")]
data_predvalues=merge(data_predvalues,nonNAHFdata[,c("UniqueSite","pred_herb")], by="UniqueSite",all.x=T)
data_predvalues=merge(data_predvalues,nonNAtd[,c("UniqueSite","pred_fundiv")], by="UniqueSite",all.x=T)

#reefs passing the thresholds for the marginalized data
sites_b20_25=length(alldata$marg_B20[alldata$marg_B20>rpred*0.25])/length(alldata$marg_B20)
sites_b20_50=length(alldata$marg_B20[alldata$marg_B20>rpred*0.50])/length(alldata$marg_B20)
sites_b20_75=length(alldata$marg_B20[alldata$marg_B20>rpred*0.75])/length(alldata$marg_B20)
sites_her_25=length(alldata$marg_ps[!is.na(alldata$marg_ps)&alldata$marg_ps>rherb*0.25])/length(alldata$marg_ps[!is.na(alldata$marg_ps)])
sites_her_50=length(alldata$marg_ps[!is.na(alldata$marg_ps)&alldata$marg_ps>rherb*0.50])/length(alldata$marg_ps[!is.na(alldata$marg_ps)])
sites_her_75=length(alldata$marg_ps[!is.na(alldata$marg_ps)&alldata$marg_ps>rherb*0.75])/length(alldata$marg_ps[!is.na(alldata$marg_ps)])
sites_fun_25=length(alldata$marg_td[!is.na(alldata$marg_td)&alldata$marg_td>rfunc*0.25])/length(alldata$marg_td[!is.na(alldata$marg_td)])
sites_fun_50=length(alldata$marg_td[!is.na(alldata$marg_td)&alldata$marg_td>rfunc*0.50])/length(alldata$marg_td[!is.na(alldata$marg_td)])
sites_fun_75=length(alldata$marg_td[!is.na(alldata$marg_td)&alldata$marg_td>rfunc*0.75])/length(alldata$marg_td[!is.na(alldata$marg_td)])


#data with all three metrics
data_predvalues2=alldata[!is.na(alldata$marg_td)&!is.na(alldata$marg_ps),]
length(data_predvalues2$UniqueSite)

#sites in our data passing the different thresholds
data_predvalues2$tot_25=ifelse(data_predvalues2$marg_td>rfunc*0.25 & data_predvalues2$marg_ps>rherb*0.25& data_predvalues2$marg_B20>rpred*0.25,1,0)
data_predvalues2$tot_50=ifelse(data_predvalues2$marg_td>rfunc*0.50 & data_predvalues2$marg_ps>rherb*0.50& data_predvalues2$marg_B20>rpred*0.50,1,0)
data_predvalues2$tot_75=ifelse(data_predvalues2$marg_td>rfunc*0.75 & data_predvalues2$marg_ps>rherb*0.75& data_predvalues2$marg_B20>rpred*0.75,1,0)


#from data with all metrics available: number of sites pasisng different thresholds
sites2_b20_25=length(data_predvalues2$marg_B20[data_predvalues2$marg_B20>rpred*0.25])/length(data_predvalues2$marg_B20)
sites2_b20_50=length(data_predvalues2$marg_B20[data_predvalues2$marg_B20>rpred*0.50])/length(data_predvalues2$marg_B20)
sites2_b20_75=length(data_predvalues2$marg_B20[data_predvalues2$marg_B20>rpred*0.75])/length(data_predvalues2$marg_B20)
sites2_her_25=length(data_predvalues2$marg_ps[!is.na(data_predvalues2$marg_ps)&data_predvalues2$marg_ps>rherb*0.25])/length(data_predvalues2$marg_ps[!is.na(data_predvalues2$marg_ps)])
sites2_her_50=length(data_predvalues2$marg_ps[!is.na(data_predvalues2$marg_ps)&data_predvalues2$marg_ps>rherb*0.50])/length(data_predvalues2$marg_ps[!is.na(data_predvalues2$marg_ps)])
sites2_her_75=length(data_predvalues2$marg_ps[!is.na(data_predvalues2$marg_ps)&data_predvalues2$marg_ps>rherb*0.75])/length(data_predvalues2$marg_ps[!is.na(data_predvalues2$marg_ps)])
sites2_fun_25=length(data_predvalues2$marg_td[!is.na(data_predvalues2$marg_td)&data_predvalues2$marg_td>rfunc*0.25])/length(data_predvalues2$marg_td[!is.na(data_predvalues2$marg_td)])
sites2_fun_50=length(data_predvalues2$marg_td[!is.na(data_predvalues2$marg_td)&data_predvalues2$marg_td>rfunc*0.50])/length(data_predvalues2$marg_td[!is.na(data_predvalues2$marg_td)])
sites2_fun_75=length(data_predvalues2$marg_td[!is.na(data_predvalues2$marg_td)&data_predvalues2$marg_td>rfunc*0.75])/length(data_predvalues2$marg_td[!is.na(data_predvalues2$marg_td)])
sites_tot_25=length(data_predvalues2$UniqueSite[data_predvalues2$marg_td>rfunc*0.25 & data_predvalues2$marg_ps>rherb*0.25& data_predvalues2$marg_B20>rpred*0.25])/length(data_predvalues2$UniqueSite)
sites_tot_50=length(data_predvalues2$UniqueSite[data_predvalues2$marg_td>rfunc*0.50 & data_predvalues2$marg_ps>rherb*0.50& data_predvalues2$marg_B20>rpred*0.50])/length(data_predvalues2$UniqueSite)
sites_tot_75=length(data_predvalues2$UniqueSite[data_predvalues2$marg_td>rfunc*0.75 & data_predvalues2$marg_ps>rherb*0.75& data_predvalues2$marg_B20>rpred*0.75])/length(data_predvalues2$UniqueSite)


#openly fished sites:proportion
fished=data_predvalues2[data_predvalues2$Protection=="Fished",]
fsites_tot_25=length(fished$UniqueSite[fished$marg_td>rfunc*0.25 & fished$marg_ps>rherb*0.25& fished$marg_B20>rpred*0.25])/length(fished$UniqueSite)
fsites_tot_50=length(fished$UniqueSite[fished$marg_td>rfunc*0.50 & fished$marg_ps>rherb*0.50& fished$marg_B20>rpred*0.50])/length(fished$UniqueSite)
fsites_tot_75=length(fished$UniqueSite[fished$marg_td>rfunc*0.75 & fished$marg_ps>rherb*0.75& fished$marg_B20>rpred*0.75])/length(fished$UniqueSite)


#create dataset to store values of sites passing different threshols
sites_passing=matrix(NA, ncol=4,nrow=36)
sites_passing[,1]=c(rep("alldata", length.out=21), rep("all_three (n=1662)", length.out=12),rep("fished_all_three(n=1092)", length.out=3))
sites_passing[,2]=c(rep("benchmarks", length.out=12), rep("b20 (n=1798)", length.out=3),rep("ps (n=1662)", length.out=3),rep("td (n=1662)", length.out=3),rep("b20", length.out=3),rep("ps", length.out=3),rep("td", length.out=3),rep("3 metrics", length.out=3),rep("3 metrics", length.out=3))
sites_passing[,3]=c("rpred","rpred25","rpred50","rpred75","rherb","rherb25","rherb50","rherb75","rfunc", "rfunc25","rfunc50","rfunc75",rep(c("propsites>25","propsites>50","propsites>75"), length.out=24))
sites_passing[,4]=c(rpred,rpred*0.25,rpred*0.5,rpred*0.75,rherb,rherb*0.25,rherb*0.50,rherb*0.75,rfunc,rfunc*0.25,rfunc*0.50,rfunc*0.75,sites_b20_25,sites_b20_50,sites_b20_75,sites_her_25,sites_her_50,sites_her_75,sites_fun_25,sites_fun_50,sites_fun_75,sites2_b20_25,sites2_b20_50,sites2_b20_75,sites2_her_25,sites2_her_50,sites2_her_75,sites2_fun_25,sites2_fun_50,sites2_fun_75,sites_tot_25,sites_tot_50,sites_tot_75,fsites_tot_25,fsites_tot_50,fsites_tot_75)
sites_passing=as.data.frame(sites_passing)
colnames(sites_passing)=c("data","variable","variable2","totexp2")
#write.csv(sites_passing, "numbermanuscript_totexp2_90margranef.csv", row.names=F)

#Biomass >20cm
#for each row, get the probability of the column
prob25 = function(x){length(x[(exp(x)-1)>rpred*0.25])/length(x)} 
prob50 = function(x){length(x[(exp(x)-1)>rpred*0.50])/length(x)} 
prob75 = function(x){length(x[(exp(x)-1)>rpred*0.75])/length(x)} 

#Create modelled data
MyData_b20=expand.grid(sgrav_tot2=seq(range(alldata$sgrav_tot2)[[1]],range(alldata$sgrav_tot2)[[2]],length=101),
                       Protection=c("Fished","Restricted","UnfishedHigh"),
                       DepthCategory=c(">10m"," 0-4m","4-10m"),
                       CleanHabitat=c("Crest","Flat","Lagoon_Back reef","Slope"),
                       sTotal_sampling_area=0,
                       CensusMethod=c("Distance sampling","Point intercept","Standard belt transect"),
                       sRegional_population_growth=0,
                       sOcean_prod=0,
                       sHDI=0,
                       sReef_fish_landings_per_km2=0,
                       sClimate_stress=0,
                       sLarger_pop_size=0)

X = model.matrix(~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot2+
                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2,data=MyData_b20)

coefs_b20=as.matrix(B20_model_totexp2)[, paste("b",row.names(fixef(B20_model_totexp2)),sep="_")] #alternative

fit_b20=coefs_b20 %*% t(X)
#for each row, get the probability of the column
try25=as.data.frame(apply(fit_b20, 2, prob25) )
colnames(try25)="estimate25"
try50=as.data.frame(apply(fit_b20, 2, prob50) )
colnames(try50)="estimate50"
try75=as.data.frame(apply(fit_b20, 2, prob75) )
colnames(try75)="estimate75"

MyData_b20=MyData_b20 %>% cbind(tidyMCMC(as.mcmc(fit_b20),
                                         conf.int=T, conf.method='HPDinterval'))
MyData_b20=cbind(MyData_b20,try25,try50,try75)

# Create Delta
dat_Fished_b20 = MyData_b20[which(MyData_b20$CleanHabitat=="Slope" & MyData_b20$Protection=="Fished" & MyData_b20$CensusMethod=="Standard belt transect" &
                                    MyData_b20$DepthCategory=="4-10m"),]
dat_Reserve_b20 = MyData_b20[which(MyData_b20$CleanHabitat=="Slope" & MyData_b20$Protection=="UnfishedHigh" & MyData_b20$CensusMethod=="Standard belt transect" &
                                     MyData_b20$DepthCategory=="4-10m"),]
dat_Rest_b20 = MyData_b20[which(MyData_b20$CleanHabitat=="Slope" & MyData_b20$Protection=="Restricted" & MyData_b20$CensusMethod=="Standard belt transect" &
                                  MyData_b20$DepthCategory=="4-10m"),]

D1 = (dat_Reserve_b20[,"estimate25"] - dat_Fished_b20[,"estimate25"])
D2 = (dat_Reserve_b20[,"estimate50"] - dat_Fished_b20[,"estimate50"])
D3 = (dat_Reserve_b20[,"estimate75"] - dat_Fished_b20[,"estimate75"])

delta_FR_b20 = cbind(dat_Fished_b20,D1,D2,D3)

D1 = (dat_Rest_b20[,"estimate25"] - dat_Fished_b20[,"estimate25"])
D2 = (dat_Rest_b20[,"estimate50"] - dat_Fished_b20[,"estimate50"])
D3 = (dat_Rest_b20[,"estimate75"] - dat_Fished_b20[,"estimate75"])

delta_FRest_b20 = cbind(dat_Fished_b20,D1,D2,D3)

D1 = (dat_Reserve_b20[,"estimate25"] - dat_Rest_b20[,"estimate25"])
D2 = (dat_Reserve_b20[,"estimate50"] - dat_Rest_b20[,"estimate50"])
D3 = (dat_Reserve_b20[,"estimate75"] - dat_Rest_b20[,"estimate75"])

delta_RestR_b20 = cbind(dat_Fished_b20,D1,D2,D3)

#now do the same for herbivore function
#for each row, get the probability of the column
prob25 = function(x){length(x[(exp(x))>rherb*0.25])/length(x)} 
prob50 = function(x){length(x[(exp(x))>rherb*0.50])/length(x)} 
prob75 = function(x){length(x[(exp(x))>rherb*0.75])/length(x)} 

MyData_hf=expand.grid(sgrav_tot2=seq(range(alldata$sgrav_tot2)[[1]],range(alldata$sgrav_tot2)[[2]],length=101),
                      Protection=c("Fished","Restricted","UnfishedHigh"),
                      DepthCategory=c(">10m"," 0-4m","4-10m"),
                      CleanHabitat=c("Crest","Flat","Lagoon_Back reef","Slope"),
                      sTotal_sampling_area=0,
                      CensusMethod=c("Distance sampling","Point intercept","Standard belt transect"),
                      sRegional_population_growth=0,
                      sOcean_prod=0,
                      
                      sHDI=0,
                      sReef_fish_landings_per_km2=0,
                      sClimate_stress=0,
                      sLarger_pop_size=0)

X = model.matrix(~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot2+
                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2,data=MyData_hf)

coefs_hf=as.matrix(her_model_totexp2)[, paste("b",row.names(fixef(her_model_totexp2)),sep="_")] #alternative

fit_hf=coefs_hf %*% t(X)
dim(fit_hf)

#for each row, get the probability of the column
try25=as.data.frame(apply(fit_hf, 2, prob25) )
colnames(try25)="estimate25"
try50=as.data.frame(apply(fit_hf, 2, prob50) )
colnames(try50)="estimate50"
try75=as.data.frame(apply(fit_hf, 2, prob75) )
colnames(try75)="estimate75"

MyData_hf=MyData_hf %>% cbind(tidyMCMC(as.mcmc(fit_hf),
                                       conf.int=T, conf.method='HPDinterval'))
MyData_hf=cbind(MyData_hf,try25,try50,try75)

# Create Delta
dat_Fished_hf = MyData_hf[which(MyData_hf$CleanHabitat=="Slope" & MyData_hf$Protection=="Fished" & MyData_hf$CensusMethod=="Standard belt transect" &
                                  MyData_hf$DepthCategory=="4-10m"),]
dat_Reserve_hf = MyData_hf[which(MyData_hf$CleanHabitat=="Slope" & MyData_hf$Protection=="UnfishedHigh" & MyData_hf$CensusMethod=="Standard belt transect" &
                                   MyData_hf$DepthCategory=="4-10m"),]
dat_Rest_hf = MyData_hf[which(MyData_hf$CleanHabitat=="Slope" & MyData_hf$Protection=="Restricted" & MyData_hf$CensusMethod=="Standard belt transect" &
                                MyData_hf$DepthCategory=="4-10m"),]

D1 = (dat_Reserve_hf[,"estimate25"] - dat_Fished_hf[,"estimate25"])
D2 = (dat_Reserve_hf[,"estimate50"] - dat_Fished_hf[,"estimate50"])
D3 = (dat_Reserve_hf[,"estimate75"] - dat_Fished_hf[,"estimate75"])

delta_FR_hf = cbind(dat_Fished_hf,D1,D2,D3)

D1 = (dat_Rest_hf[,"estimate25"] - dat_Fished_hf[,"estimate25"])
D2 = (dat_Rest_hf[,"estimate50"] - dat_Fished_hf[,"estimate50"])
D3 = (dat_Rest_hf[,"estimate75"] - dat_Fished_hf[,"estimate75"])

delta_FRest_hf = cbind(dat_Fished_hf,D1,D2,D3)

D1 = (dat_Reserve_hf[,"estimate25"] - dat_Rest_hf[,"estimate25"])
D2 = (dat_Reserve_hf[,"estimate50"] - dat_Rest_hf[,"estimate50"])
D3 = (dat_Reserve_hf[,"estimate75"] - dat_Rest_hf[,"estimate75"])

delta_RestR_hf = cbind(dat_Fished_hf,D1,D2,D3)

#trait diversity
prob25 = function(x){length(x[exp(x)>rfunc*0.25])/length(x)} 
prob50 = function(x){length(x[exp(x)>rfunc*0.50])/length(x)} 
prob75 = function(x){length(x[exp(x)>rfunc*0.75])/length(x)} 

MyData_td=expand.grid(sgrav_tot2=seq(range(alldata$sgrav_tot2)[[1]],range(alldata$sgrav_tot2)[[2]],length=101),
                      Protection=c("Fished","Restricted","UnfishedHigh"),
                      DepthCategory=c(">10m"," 0-4m","4-10m"),
                      CleanHabitat=c("Crest","Flat","Lagoon_Back reef","Slope"),
                      sTotal_sampling_area=0,
                      CensusMethod=c("Distance sampling","Point intercept","Standard belt transect"),
                      sRegional_population_growth=0,
                      sOcean_prod=0,
                      sHDI=0,
                      sReef_fish_landings_per_km2=0,
                      sClimate_stress=0,
                      sLarger_pop_size=0)

X = model.matrix(~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot2+
                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2,data=MyData_td)

coefs_td=as.matrix(fundiv_model_totexp2)[, paste("b",row.names(fixef(fundiv_model_totexp2)),sep="_")] #alternative

fit_td=coefs_td %*% t(X)

#for each row, get the probability of the column
try25=as.data.frame(apply(fit_td, 2, prob25) )
colnames(try25)="estimate25"
try50=as.data.frame(apply(fit_td, 2, prob50) )
colnames(try50)="estimate50"
try75=as.data.frame(apply(fit_td, 2, prob75) )
colnames(try75)="estimate75"

MyData_td=MyData_td %>% cbind(tidyMCMC(as.mcmc(fit_td),
                                       conf.int=T, conf.method='HPDinterval'))
MyData_td=cbind(MyData_td,try25,try50,try75)

# Create Delta
dat_Fished_td = MyData_td[which(MyData_td$CleanHabitat=="Slope" & MyData_td$Protection=="Fished" & MyData_td$CensusMethod=="Standard belt transect" &
                                  MyData_td$DepthCategory=="4-10m"),]
dat_Reserve_td = MyData_td[which(MyData_td$CleanHabitat=="Slope" & MyData_td$Protection=="UnfishedHigh" & MyData_td$CensusMethod=="Standard belt transect" &
                                   MyData_td$DepthCategory=="4-10m"),]
dat_Rest_td = MyData_td[which(MyData_td$CleanHabitat=="Slope" & MyData_td$Protection=="Restricted" & MyData_td$CensusMethod=="Standard belt transect" &
                                MyData_td$DepthCategory=="4-10m"),]

D1 = (dat_Reserve_td[,"estimate25"] - dat_Fished_td[,"estimate25"])
D2 = (dat_Reserve_td[,"estimate50"] - dat_Fished_td[,"estimate50"])
D3 = (dat_Reserve_td[,"estimate75"] - dat_Fished_td[,"estimate75"])

delta_FR_td = cbind(dat_Fished_td,D1,D2,D3)

D1 = (dat_Rest_td[,"estimate25"] - dat_Fished_td[,"estimate25"])
D2 = (dat_Rest_td[,"estimate50"] - dat_Fished_td[,"estimate50"])
D3 = (dat_Rest_td[,"estimate75"] - dat_Fished_td[,"estimate75"])

delta_FRest_td = cbind(dat_Fished_td,D1,D2,D3)

D1 = (dat_Reserve_td[,"estimate25"] - dat_Rest_td[,"estimate25"])
D2 = (dat_Reserve_td[,"estimate50"] - dat_Rest_td[,"estimate50"])
D3 = (dat_Reserve_td[,"estimate75"] - dat_Rest_td[,"estimate75"])

delta_RestR_td = cbind(dat_Fished_td,D1,D2,D3)


#all three goals (assuming independence (multiply outcomes))
# Create Delta
dat_Fished_3 = as.data.frame(cbind((dat_Fished_b20$estimate25*dat_Fished_hf$estimate25*dat_Fished_td$estimate25),(dat_Fished_b20$estimate50*dat_Fished_hf$estimate50*dat_Fished_td$estimate50),(dat_Fished_b20$estimate75*dat_Fished_hf$estimate75*dat_Fished_td$estimate75),dat_Fished_b20$sgrav_tot2))
colnames(dat_Fished_3 )  =c("estimate25","estimate50","estimate75","sgrav_tot2")
dat_Reserve_3 =  as.data.frame(cbind((dat_Reserve_b20$estimate25*dat_Reserve_hf$estimate25*dat_Reserve_td$estimate25),(dat_Reserve_b20$estimate50*dat_Reserve_hf$estimate50*dat_Reserve_td$estimate50),(dat_Reserve_b20$estimate75*dat_Reserve_hf$estimate75*dat_Reserve_td$estimate75),dat_Reserve_b20$sgrav_tot2))
colnames(dat_Reserve_3 )  =c("estimate25","estimate50","estimate75","sgrav_tot2")
dat_Rest_3 = as.data.frame(cbind((dat_Rest_b20$estimate25*dat_Rest_hf$estimate25*dat_Rest_td$estimate25),(dat_Rest_b20$estimate50*dat_Rest_hf$estimate50*dat_Rest_td$estimate50),(dat_Rest_b20$estimate75*dat_Rest_hf$estimate75*dat_Rest_td$estimate75),dat_Rest_b20$sgrav_tot2))
colnames(dat_Rest_3 )  =c("estimate25","estimate50","estimate75","sgrav_tot2")

D1 = (dat_Reserve_3[,"estimate25"] - dat_Fished_3[,"estimate25"])
D2 = (dat_Reserve_3[,"estimate50"] - dat_Fished_3[,"estimate50"])
D3 = (dat_Reserve_3[,"estimate75"] - dat_Fished_3[,"estimate75"])

delta_FR_3 = cbind(dat_Fished_3,D1,D2,D3)

D1 = (dat_Rest_3[,"estimate25"] - dat_Fished_3[,"estimate25"])
D2 = (dat_Rest_3[,"estimate50"] - dat_Fished_3[,"estimate50"])
D3 = (dat_Rest_3[,"estimate75"] - dat_Fished_3[,"estimate75"])

delta_FRest_3 = cbind(dat_Fished_3,D1,D2,D3)

D1 = (dat_Reserve_3[,"estimate25"] - dat_Rest_3[,"estimate25"])
D2 = (dat_Reserve_3[,"estimate50"] - dat_Rest_3[,"estimate50"])
D3 = (dat_Reserve_3[,"estimate75"] - dat_Rest_3[,"estimate75"])

delta_RestR_3 = cbind(dat_Fished_3,D1,D2,D3)


#probability figure

p5_bmrf = ggplot(delta_FR_b20)+
  geom_line(aes(sgrav_tot2,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.05,0.8),labels=NULL) + theme_classic() 

p6_brf = ggplot(delta_FRest_b20)+
  geom_line(aes(sgrav_tot2,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.justification=c(1,0), legend.position=c(0,1), line= element_blank(), axis.ticks = element_line()) +
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.05,0.8),labels=NULL) + theme_classic()

p7_brr = ggplot(delta_RestR_b20)+
  geom_line(aes(sgrav_tot2,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.justification=c(1,0), legend.position=c(0,1), line= element_blank(), axis.ticks = element_line()) +
  scale_x_continuous("")+
  scale_y_continuous("",limits=c(-0.05,0.8),labels=NULL) + theme_classic()

p1_bfci = ggplot(MyData_b20[which(MyData_b20$CleanHabitat=="Slope" & MyData_b20$Protection=="Fished" & MyData_b20$CensusMethod=="Standard belt transect" &
                                    MyData_b20$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_tot2,estimate25),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,estimate50),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + theme_classic()

p2_bmrci = ggplot(MyData_b20[which(MyData_b20$CleanHabitat=="Slope" & MyData_b20$Protection=="UnfishedHigh" & MyData_b20$CensusMethod=="Standard belt transect" &
                                     MyData_b20$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_tot2,estimate25),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,estimate50),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + theme_classic()


p3_brci = ggplot(MyData_b20[which(MyData_b20$CleanHabitat=="Slope" & MyData_b20$Protection=="Restricted" & MyData_b20$CensusMethod=="Standard belt transect" &
                                    MyData_b20$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_tot2,estimate25),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,estimate50),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + theme_classic()


# figures for manuscript
p5_hmrf = ggplot(delta_FR_hf)+
  geom_line(aes(sgrav_tot2,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.05,0.8),labels=NULL) + theme_classic() 

p6_hrf = ggplot(delta_FRest_hf)+
  geom_line(aes(sgrav_tot2,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.justification=c(1,0), legend.position=c(0,1), line= element_blank(), axis.ticks = element_line()) +
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.05,0.8),labels=NULL) + theme_classic()

p7_hrr = ggplot(delta_RestR_hf)+
  geom_line(aes(sgrav_tot2,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.justification=c(1,0), legend.position=c(0,1), line= element_blank(), axis.ticks = element_line()) +
  scale_x_continuous("")+
  scale_y_continuous("",limits=c(-0.05,0.8),labels=NULL) + theme_classic()

p1_hfci = ggplot(MyData_hf[which(MyData_hf$CleanHabitat=="Slope" & MyData_hf$Protection=="Fished" & MyData_hf$CensusMethod=="Standard belt transect" &
                                   MyData_hf$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_tot2,estimate25),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,estimate50),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + theme_classic()

p2_hmrci = ggplot(MyData_hf[which(MyData_hf$CleanHabitat=="Slope" & MyData_hf$Protection=="UnfishedHigh" & MyData_hf$CensusMethod=="Standard belt transect" &
                                    MyData_hf$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_tot2,estimate25),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,estimate50),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + theme_classic()


p3_hrci = ggplot(MyData_hf[which(MyData_hf$CleanHabitat=="Slope" & MyData_hf$Protection=="Restricted" & MyData_hf$CensusMethod=="Standard belt transect" &
                                   MyData_hf$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_tot2,estimate25),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,estimate50),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + theme_classic()


# figures for manuscript
p5_rmrf = ggplot(delta_FR_td)+
  geom_line(aes(sgrav_tot2,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.05,0.8),labels=NULL) + theme_classic() 

p6_rrf = ggplot(delta_FRest_td)+
  geom_line(aes(sgrav_tot2,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.justification=c(1,0), legend.position=c(0,1), line= element_blank(), axis.ticks = element_line()) +
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.05,0.8),labels=NULL) + theme_classic()

p7_rrr = ggplot(delta_RestR_td)+
  geom_line(aes(sgrav_tot2,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.justification=c(1,0), legend.position=c(0,1), line= element_blank(), axis.ticks = element_line()) +
  scale_x_continuous("")+
  scale_y_continuous("",limits=c(-0.05,0.8),labels=NULL) + theme_classic()

p1_rfci = ggplot(MyData_td[which(MyData_td$CleanHabitat=="Slope" & MyData_td$Protection=="Fished" & MyData_td$CensusMethod=="Standard belt transect" &
                                   MyData_td$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_tot2,estimate25),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,estimate50),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + theme_classic()

p2_rmrci = ggplot(MyData_td[which(MyData_td$CleanHabitat=="Slope" & MyData_td$Protection=="UnfishedHigh" & MyData_td$CensusMethod=="Standard belt transect" &
                                    MyData_td$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_tot2,estimate25),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,estimate50),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + theme_classic()


p3_rrci = ggplot(MyData_td[which(MyData_td$CleanHabitat=="Slope" & MyData_td$Protection=="Restricted" & MyData_td$CensusMethod=="Standard belt transect" &
                                   MyData_td$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_tot2,estimate25),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,estimate50),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + theme_classic()


p5_3mrf = ggplot(delta_FR_3)+
  geom_line(aes(sgrav_tot2,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.05,0.8)) + theme_classic() 

p6_3rf = ggplot(delta_FRest_3)+
  geom_line(aes(sgrav_tot2,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.justification=c(1,0), legend.position=c(0,1), line= element_blank(), axis.ticks = element_line()) +
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.05,0.8)) + theme_classic()

p7_3rr = ggplot(delta_RestR_3)+
  geom_line(aes(sgrav_tot2,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.justification=c(1,0), legend.position=c(0,1), line= element_blank(), axis.ticks = element_line()) +
  scale_x_continuous("")+
  scale_y_continuous("",limits=c(-0.05,0.8)) + theme_classic()

p1_3fci = ggplot(dat_Fished_3)+
  geom_line(aes(sgrav_tot2,estimate25),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,estimate50),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1)) + theme_classic()

p2_3mrci = ggplot(dat_Reserve_3)+
  geom_line(aes(sgrav_tot2,estimate25),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,estimate50),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1)) + theme_classic()


p3_3rci = ggplot(dat_Rest_3)+
  geom_line(aes(sgrav_tot2,estimate25),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_tot2,estimate50),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_tot2,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1)) + theme_classic()

#manuscript
windows()
Fig.grouped=ggarrange(p1_3fci, p1_bfci,p1_hfci,p1_rfci,p2_3mrci,p2_bmrci,p2_hmrci,p2_rmrci,p3_3rci,p3_brci,p3_hrci,p3_rrci,p5_3mrf,p5_bmrf,p5_hmrf,p5_rmrf,p6_3rf,p6_brf,p6_hrf,p6_rrf, p7_3rr,p7_brr,p7_hrr,p7_rrr,nrow=6, ncol=4, widths = c(1.2,1,1,1), heights=c(1,1,1,1,1,1.2), labels=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X"),font.label = list(size = 10, color = "black", face = "bold", family = "Helvetica"))
annotate_figure(Fig.grouped,left = text_grob("Difference in probability                                                     Probability", rot = 90),
                bottom=text_grob("Std Gravity (log +min transformed)"))

#######################################################################################################
##maps of marginalized response variables for manuscript
#correlations
pdf("FigureS1.pdf", width=7.25,family="Helvetica")
pairs(~log(marg_B20_habitat+1)+log(marg_td_habitat) +log(marg_ps_habitat+1),data=alldata,lower.panel=panel.cor, 
      pch = 21, bg = "darkgrey",labels=c("log(Biomass >20cm+1)","log(Trait diversity)", "log(Parrotfish
scraping+1)"),cex.labels=1.5,font.labels=2,diag.panel =panel_hist, hist.col="grey")
dev.off()


#extract map and match latitude
mapWorld <- map_data('world', wrap=c(-25,335), ylim=c(-55,75))
alldata$lon2 <- ifelse(alldata$Site_Long2 < -25, alldata$Site_Long2 + 360, alldata$Site_Long2) 
#specify gradient colours
mycolours = c("dodgerblue2", "lightskyblue1", "white", "darkorange", "orangered3") 

Fig.1a=ggplot() + 
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group),fill = "grey", color = "grey") +
  coord_fixed(xlim = c(30, 320),  ylim = c(30, -30), ratio = 1.3)+
  geom_point(data=alldata[!is.na(alldata$marg_ps),], aes(x=lon2, y=Site_Lat2, fill=log(alldata$marg_ps_habitat[!is.na(alldata$marg_ps_habitat)]+1)),colour="black", pch=21, size=3)+
  scale_fill_gradientn(colors = mycolours, name="log(Scraping
potential +1)",guide=F)+  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+scale_x_continuous("",breaks=c(80,180,280),labels=c(-100,0,100))+
  scale_y_continuous("",breaks=c(-20,-10,0,10,20))+ theme_classic()

Fig.1b=ggplot() + 
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group),fill = "grey", color = "grey") +
  coord_fixed(xlim = c(30, 320),  ylim = c(30, -30), ratio = 1.3)+
  geom_point(data=alldata, aes(x=lon2, y=Site_Lat2, fill=log(alldata$marg_B20_habitat+1)),colour="black", pch=21, size=3)+
  scale_fill_gradientn(colors = mycolours, name="log (Biomass
> 20 cm +1)",guide=F)+  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+scale_x_continuous("",breaks=c(80,180,280),labels=c(-100,0,100))+
  scale_y_continuous("",breaks=c(-20,-10,0,10,20))+ theme_classic()

Fig.1c=ggplot() + 
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group),fill = "grey", color = "grey") +
  coord_fixed(xlim = c(30, 320),  ylim = c(30, -30), ratio = 1.3)+
  geom_point(data=alldata[!is.na(alldata$marg_td),], aes(x=lon2, y=Site_Lat2, fill=log(alldata$marg_td_habitat[!is.na(alldata$marg_td_habitat)] )),colour="black", pch=21, size=3)+
  scale_fill_gradientn(colors = mycolours, name="log(Trait
diversity)",guide=F)+  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+scale_x_continuous("",breaks=c(80,180,280),labels=c(-100,0,100))+
  scale_y_continuous("",breaks=c(-20,-10,0,10,20))+ theme_classic()

data_predvalues2$category3goals=ifelse(data_predvalues2$tot_25==0,"< 25", ifelse(data_predvalues2$tot_25==1 &data_predvalues2$tot_50==0&data_predvalues2$tot_75==0,"25-50",ifelse(data_predvalues2$tot_50==1 &data_predvalues2$tot_75==0,"50-75","> 75")))
FigureX=ggplot() + 
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group),fill = "grey", color = "grey") +
  coord_fixed(xlim = c(30, 320),  ylim = c(30, -30), ratio = 1.3)+
  geom_point(data=data_predvalues2, aes(x=lon2, y=Site_Lat2, fill=data_predvalues2$category3goals, col=data_predvalues2$category3goals), pch=21, size=3)+
  scale_fill_manual(values=c("> 75"="#7a0177","50-75"="#c51b8a","25-50"="#f768a1","< 25"="black"),guide=F)+
  scale_colour_manual(values=c("> 75"="#7a0177","50-75"="#c51b8a","25-50"="#f768a1","< 25"="black"),guide=F)+
  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+scale_x_continuous("",breaks=c(80,180,280),labels=c(-100,0,100))+
  scale_y_continuous("",breaks=c(-20,-10,0,10,20))+ theme_classic()

Fig.1=ggarrange(Fig.1b,Fig.1a,Fig.1c,FigureX,nrow=4, ncol=1)
windows()
Fig.1
#####################################################################
#conservation gains for our sites given their socio-ecological context

#we simmulate new data from the posterior maintaining sampling constant
test=alldata
test[1:dim(test)[1],"CensusMethod"]="Standard belt transect"
test[1:dim(test)[1],"CleanHabitat"]="Slope"
test[1:dim(test)[1],"DepthCategory"]="4-10m"
test[1:dim(test)[1],"sTotal_sampling_area"]=0
test=droplevels(test)
B20test_posterior=posterior_predict(B20_model_totexp2,newdata=test,re_formula = NULL)
alldata$B20_newdata=exp(colMeans(B20test_posterior))-1

#we simmulate new data from the posterior maintaining sampling constant but making all sites reserves
reserves=test
reserves[1:dim(reserves)[1],"Protection"]="UnfishedHigh"
reserves=droplevels(reserves)
B20reserves_posterior=posterior_predict(B20_model_totexp2,newdata=reserves,re_formula = NULL)
alldata$B20_newdata_reserves=exp(colMeans(B20reserves_posterior))-1

#we simmulate new data from the posterior maintaining sampling constant but making all sites restricted
restricted=test
restricted[1:dim(restricted)[1],"Protection"]="Restricted"
restricted=droplevels(restricted)
B20restricted_posterior=posterior_predict(B20_model_totexp2,newdata=restricted,re_formula = NULL)
alldata$B20_newdata_restricted=exp(colMeans(B20restricted_posterior))-1

#to calculate the probabilities of pasisng differnet threshols
prob25 <- function(x){length(x[(exp(x)-1)>rpred*0.25])/length(x)} 
prob50 <- function(x){length(x[(exp(x)-1)>rpred*0.50])/length(x)} 
prob75 <- function(x){length(x[(exp(x)-1)>rpred*0.75])/length(x)} 

#for each row, get the probability of the column
newdata25=as.data.frame(apply(B20test_posterior, 2, prob25) )
colnames(newdata25)="B20newdata25"
newdata50=as.data.frame(apply(B20test_posterior, 2, prob50) )
colnames(newdata50)="B20newdata50"
newdata75=as.data.frame(apply(B20test_posterior, 2, prob75) )
colnames(newdata75)="B20newdata75"
newdata=cbind(alldata,newdata25,newdata50,newdata75)

reserves25=as.data.frame(apply(B20reserves_posterior, 2, prob25) )
colnames(reserves25)="B20reserves25"
reserves50=as.data.frame(apply(B20reserves_posterior, 2, prob50) )
colnames(reserves50)="B20reserves50"
reserves75=as.data.frame(apply(B20reserves_posterior, 2, prob75) )
colnames(reserves75)="B20reserves75"
newdata=cbind(newdata,reserves25,reserves50,reserves75)

restricted25=as.data.frame(apply(B20restricted_posterior, 2, prob25) )
colnames(restricted25)="B20restricted25"
restricted50=as.data.frame(apply(B20restricted_posterior, 2, prob50) )
colnames(restricted50)="B20restricted50"
restricted75=as.data.frame(apply(B20restricted_posterior, 2, prob75) )
colnames(restricted75)="B20restricted75"
newdata=cbind(newdata,restricted25,restricted50,restricted75)
colnames(newdata)

# we do the same for trait diversity
fundivtest_posterior=posterior_predict(fundiv_model_totexp2,newdata=test[!is.na(test$Trait_diversity),],re_formula = NULL)
newdata[!is.na(newdata$Trait_diversity),"fundiv_newdata"]=exp(colMeans(fundivtest_posterior))
fundivreserves_posterior=posterior_predict(fundiv_model_totexp2,newdata=reserves[!is.na(reserves$Trait_diversity),],re_formula = NULL)
newdata[!is.na(newdata$Trait_diversity),"fundiv_newdata_reserves"]=exp(colMeans(fundivreserves_posterior))
fundivrestricted_posterior=posterior_predict(fundiv_model_totexp2,newdata=restricted[!is.na(restricted$Trait_diversity),],re_formula = NULL)
newdata[!is.na(newdata$Trait_diversity),"fundiv_newdata_restricted"]=exp(colMeans(fundivrestricted_posterior))

prob25 <- function(x){length(x[exp(x)>rfunc*0.25])/length(x)} 
prob50 <- function(x){length(x[exp(x)>rfunc*0.50])/length(x)} 
prob75 <- function(x){length(x[exp(x)>rfunc*0.75])/length(x)} 

newdata[!is.na(newdata$Trait_diversity),"fundivnewdata25"]=as.data.frame(apply(fundivtest_posterior, 2, prob25) )
newdata[!is.na(newdata$Trait_diversity),"fundivnewdata50"]=as.data.frame(apply(fundivtest_posterior, 2, prob50) )
newdata[!is.na(newdata$Trait_diversity),"fundivnewdata75"]=as.data.frame(apply(fundivtest_posterior, 2, prob75) )


newdata[!is.na(newdata$Trait_diversity),"fundivreserves25"]=as.data.frame(apply(fundivreserves_posterior, 2, prob25) )
newdata[!is.na(newdata$Trait_diversity),"fundivreserves50"]=as.data.frame(apply(fundivreserves_posterior, 2, prob50) )
newdata[!is.na(newdata$Trait_diversity),"fundivreserves75"]=as.data.frame(apply(fundivreserves_posterior, 2, prob75) )

newdata[!is.na(newdata$Trait_diversity),"fundivrestricted25"]=as.data.frame(apply(fundivrestricted_posterior, 2, prob25) )
newdata[!is.na(newdata$Trait_diversity),"fundivrestricted50"]=as.data.frame(apply(fundivrestricted_posterior, 2, prob50) )
newdata[!is.na(newdata$Trait_diversity),"fundivrestricted75"]=as.data.frame(apply(fundivrestricted_posterior, 2, prob75) )

# we do the same for parrotfish sraping
hertest_posterior=posterior_predict(her_model_totexp2,newdata=test[!is.na(test$Scraping_potential),],re_formula = NULL)
newdata[!is.na(newdata$Scraping_potential),"her_newdata"]=colMeans(hertest_posterior)
herreserves_posterior=posterior_predict(her_model_totexp2,newdata=reserves[!is.na(reserves$Scraping_potential),],re_formula = NULL)
newdata[!is.na(newdata$Scraping_potential),"her_newdata_reserves"]=colMeans(herreserves_posterior)
herrestricted_posterior=posterior_predict(her_model_totexp2,newdata=restricted[!is.na(restricted$Scraping_potential),],re_formula = NULL)
newdata[!is.na(newdata$Scraping_potential),"her_newdata_restricted"]=colMeans(herrestricted_posterior)

prob25 <- function(x){length(x[x>rherb*0.25])/length(x)} 
prob50 <- function(x){length(x[x>rherb*0.50])/length(x)} 
prob75 <- function(x){length(x[x>rherb*0.75])/length(x)} 

newdata[!is.na(newdata$Scraping_potential),"hernewdata25"]=as.data.frame(apply(hertest_posterior, 2, prob25) )
newdata[!is.na(newdata$Scraping_potential),"hernewdata50"]=as.data.frame(apply(hertest_posterior, 2, prob50) )
newdata[!is.na(newdata$Scraping_potential),"hernewdata75"]=as.data.frame(apply(hertest_posterior, 2, prob75) )


newdata[!is.na(newdata$Scraping_potential),"herreserves25"]=as.data.frame(apply(herreserves_posterior, 2, prob25) )
newdata[!is.na(newdata$Scraping_potential),"herreserves50"]=as.data.frame(apply(herreserves_posterior, 2, prob50) )
newdata[!is.na(newdata$Scraping_potential),"herreserves75"]=as.data.frame(apply(herreserves_posterior, 2, prob75) )

newdata[!is.na(newdata$Scraping_potential),"herrestricted25"]=as.data.frame(apply(herrestricted_posterior, 2, prob25) )
newdata[!is.na(newdata$Scraping_potential),"herrestricted50"]=as.data.frame(apply(herrestricted_posterior, 2, prob50) )
newdata[!is.na(newdata$Scraping_potential),"herrestricted75"]=as.data.frame(apply(herrestricted_posterior, 2, prob75) )

#all three together
newdata$totnewdata25=newdata$B20newdata25*newdata$hernewdata25*newdata$fundivnewdata25
newdata$totnewdata50=newdata$B20newdata50*newdata$hernewdata50*newdata$fundivnewdata50
newdata$totnewdata75=newdata$B20newdata75*newdata$hernewdata75*newdata$fundivnewdata75
newdata$totreserves25=newdata$B20reserves25*newdata$herreserves25*newdata$fundivreserves25
newdata$totreserves50=newdata$B20reserves50*newdata$herreserves50*newdata$fundivreserves50
newdata$totreserves75=newdata$B20reserves75*newdata$herreserves75*newdata$fundivreserves75
newdata$totrestricted25=newdata$B20restricted25*newdata$herrestricted25*newdata$fundivrestricted25
newdata$totrestricted50=newdata$B20restricted50*newdata$herrestricted50*newdata$fundivrestricted50
newdata$totrestricted75=newdata$B20restricted75*newdata$herrestricted75*newdata$fundivrestricted75

#newdata with all three metrics present
newdata2=newdata[!is.na(newdata$totnewdata75),]

#originally fished sites only
fished_newdata=newdata[newdata$Protection=="Fished",]

#B20
fished_newdata$predicted_pred25=ifelse(fished_newdata$B20_newdata<rpred*0.25,0,1)
fished_newdata$predicted_pred50=ifelse(fished_newdata$B20_newdata<rpred*0.50,0,1)
fished_newdata$predicted_pred75=ifelse(fished_newdata$B20_newdata<rpred*0.75,0,1)
fished_newdata$unfished_pred25=ifelse(fished_newdata$B20_newdata_reserves<rpred*0.25,0,1)
fished_newdata$unfished_pred50=ifelse(fished_newdata$B20_newdata_reserves<rpred*0.50,0,1)
fished_newdata$unfished_pred75=ifelse(fished_newdata$B20_newdata_reserves<rpred*0.75,0,1)
fished_newdata$restricted_pred25=ifelse(fished_newdata$B20_newdata_restricted<rpred*0.25,0,1)
fished_newdata$restricted_pred50=ifelse(fished_newdata$B20_newdata_restricted<rpred*0.50,0,1)
fished_newdata$restricted_pred75=ifelse(fished_newdata$B20_newdata_restricted<rpred*0.75,0,1)


fished_newdata$unfished_changeB20=ifelse(fished_newdata$predicted_pred25==0 &fished_newdata$unfished_pred75==1, "0to75", 
                                         ifelse(fished_newdata$predicted_pred25==0 & fished_newdata$unfished_pred50==1 & fished_newdata$unfished_pred75==0, "0to50",
                                                ifelse(fished_newdata$predicted_pred25==0 & fished_newdata$unfished_pred25==1 & fished_newdata$unfished_pred50==0, "0to25",
                                                       ifelse(fished_newdata$predicted_pred25==0 & fished_newdata$unfished_pred25==0, "0to0",
                                                              ifelse(fished_newdata$predicted_pred25==1 &fished_newdata$predicted_pred50==0 &fished_newdata$unfished_pred75==1, "25to75",
                                                                     ifelse(fished_newdata$predicted_pred25==1 &fished_newdata$predicted_pred50==1&fished_newdata$predicted_pred75==0 &fished_newdata$unfished_pred75==1, "50to75",
                                                                            ifelse(fished_newdata$predicted_pred25==1 &fished_newdata$predicted_pred50==1&fished_newdata$predicted_pred75==1 &fished_newdata$unfished_pred75==1, "75to75",
                                                                                   ifelse(fished_newdata$predicted_pred25==1 &fished_newdata$predicted_pred50==0 &fished_newdata$unfished_pred75==0&fished_newdata$unfished_pred50==1, "25to50",
                                                                                          ifelse(fished_newdata$predicted_pred25==1 &fished_newdata$predicted_pred50==1 &fished_newdata$unfished_pred75==0&fished_newdata$unfished_pred50==1, "50to50",
                                                                                                 ifelse(fished_newdata$predicted_pred25==1 &fished_newdata$predicted_pred50==0 &fished_newdata$unfished_pred25==1 &fished_newdata$unfished_pred75==0&fished_newdata$unfished_pred50==0, "25to25",-999))))))))))
fished_newdata$restricted_changeB20=ifelse(fished_newdata$predicted_pred25==0 &fished_newdata$restricted_pred75==1, "0to75", 
                                           ifelse(fished_newdata$predicted_pred25==0 & fished_newdata$restricted_pred50==1 & fished_newdata$restricted_pred75==0, "0to50",
                                                  ifelse(fished_newdata$predicted_pred25==0 & fished_newdata$restricted_pred25==1 & fished_newdata$restricted_pred50==0, "0to25",
                                                         ifelse(fished_newdata$predicted_pred25==0 & fished_newdata$restricted_pred25==0, "0to0",
                                                                ifelse(fished_newdata$predicted_pred25==1 &fished_newdata$predicted_pred50==0 &fished_newdata$restricted_pred75==1, "25to75",
                                                                       ifelse(fished_newdata$predicted_pred25==1 &fished_newdata$predicted_pred50==1&fished_newdata$predicted_pred75==0 &fished_newdata$restricted_pred75==1, "50to75",
                                                                              ifelse(fished_newdata$predicted_pred25==1 &fished_newdata$predicted_pred50==1&fished_newdata$predicted_pred75==1 &fished_newdata$restricted_pred75==1, "75to75",
                                                                                     ifelse(fished_newdata$predicted_pred25==1 &fished_newdata$predicted_pred50==0 &fished_newdata$restricted_pred75==0&fished_newdata$restricted_pred50==1, "25to50",
                                                                                            ifelse(fished_newdata$predicted_pred25==1 &fished_newdata$predicted_pred50==1 &fished_newdata$restricted_pred75==0&fished_newdata$restricted_pred50==1, "50to50",
                                                                                                   ifelse(fished_newdata$predicted_pred25==1 &fished_newdata$predicted_pred50==0 &fished_newdata$restricted_pred25==1 &fished_newdata$restricted_pred75==0&fished_newdata$restricted_pred50==0, "25to25",-999))))))))))


#crete matric of conservation gains to use for the alluvial plots
matrixpred=matrix(NA,nrow=24,ncol=5)
matrixpred[1:4,1]=rep(0,4)
matrixpred[1:4,2]=c(1:4)
matrixpred[1:4,3]=c(length(fished_newdata$B20_newdata[fished_newdata$B20_newdata>rpred*0.75]),length(fished_newdata$B20_newdata[fished_newdata$B20_newdata>rpred*0.50&fished_newdata$B20_newdata<rpred*0.75]),length(fished_newdata$B20_newdata[fished_newdata$B20_newdata>rpred*0.25&fished_newdata$B20_newdata<rpred*0.50]),length(fished_newdata$B20_newdata[fished_newdata$B20_newdata<rpred*0.25]))
matrixpred[1:4,4]=rep("Initial", 4)
matrixpred[1:4,5]=c("75-100%","50-75%","25-50%","0-25%")
matrixpred[5:8,1]=rep(4,4)
matrixpred[5:8,2]=c(8,7,6,5)
matrixpred[5:8,4]=rep("reserve",4)
matrixpred[5:8,5]=c("0to0","0to25","0to50","0to75")
matrixpred[5:8,3]=c(length(fished_newdata$unfished_changeB20[fished_newdata$unfished_changeB20=="0to0"]),length(fished_newdata$unfished_changeB20[fished_newdata$unfished_changeB20=="0to25"]),length(fished_newdata$unfished_changeB20[fished_newdata$unfished_changeB20=="0to50"]),length(fished_newdata$unfished_changeB20[fished_newdata$unfished_changeB20=="0to75"]))
matrixpred[9:11,1]=rep(3,3)
matrixpred[9:11,2]=c(7,6,5)
matrixpred[9:11,4]=rep("reserve",3)
matrixpred[9:11,5]=c("25to25","25to50","25to75")
matrixpred[9:11,3]=c(length(fished_newdata$unfished_changeB20[fished_newdata$unfished_changeB20=="25to25"]),length(fished_newdata$unfished_changeB20[fished_newdata$unfished_changeB20=="25to50"]),length(fished_newdata$unfished_changeB20[fished_newdata$unfished_changeB20=="25to75"]))
matrixpred[12:13,1]=rep(2,2)
matrixpred[12:13,2]=c(6,5)
matrixpred[12:13,4]=rep("reserve",2)
matrixpred[12:13,5]=c("50to50","50to75")
matrixpred[12:13,3]=c(length(fished_newdata$unfished_changeB20[fished_newdata$unfished_changeB20=="50to50"]),length(fished_newdata$unfished_changeB20[fished_newdata$unfished_changeB20=="50to75"]))
matrixpred[14,1]=1
matrixpred[14,2]=5
matrixpred[14,4]="reserve"
matrixpred[14,5]="75to75"
matrixpred[14,3]=length(fished_newdata$unfished_changeB20[fished_newdata$unfished_changeB20=="75to75"])
matrixpred[15:18,1]=rep(4,4)
matrixpred[15:18,2]=c(8,7,6,5)
matrixpred[15:18,4]=rep("restricted",4)
matrixpred[15:18,5]=c("0to0","0to25","0to50","0to75")
matrixpred[15:18,3]=c(length(fished_newdata$restricted_changeB20[fished_newdata$restricted_changeB20=="0to0"]),length(fished_newdata$restricted_changeB20[fished_newdata$restricted_changeB20=="0to25"]),length(fished_newdata$restricted_changeB20[fished_newdata$restricted_changeB20=="0to50"]),length(fished_newdata$restricted_changeB20[fished_newdata$restricted_changeB20=="0to75"]))
matrixpred[19:21,1]=rep(3,3)
matrixpred[19:21,2]=c(7,6,5)
matrixpred[19:21,4]=rep("restricted",3)
matrixpred[19:21,5]=c("25to25","25to50","25to75")
matrixpred[19:21,3]=c(length(fished_newdata$restricted_changeB20[fished_newdata$restricted_changeB20=="25to25"]),length(fished_newdata$restricted_changeB20[fished_newdata$restricted_changeB20=="25to50"]),length(fished_newdata$restricted_changeB20[fished_newdata$restricted_changeB20=="25to75"]))
matrixpred[22:23,1]=rep(2,2)
matrixpred[22:23,2]=c(6,5)
matrixpred[22:23,4]=rep("restricted",2)
matrixpred[22:23,5]=c("50to50","50to75")
matrixpred[22:23,3]=c(length(fished_newdata$restricted_changeB20[fished_newdata$restricted_changeB20=="50to50"]),length(fished_newdata$restricted_changeB20[fished_newdata$restricted_changeB20=="50to75"]))
matrixpred[24,1]=1
matrixpred[24,2]=5
matrixpred[24,4]="restricted"
matrixpred[24,5]="75to75"
matrixpred[24,3]=length(fished_newdata$restricted_changeB20[fished_newdata$restricted_changeB20=="75to75"])


#alluvial charts for B20
pred=as.data.frame(matrixpred)
colnames(pred)=c("from","to","value","scenario","clarified")
pred$from=as.numeric(as.character(pred$from))
pred$to=as.numeric(as.character(pred$to))
pred$value=as.numeric(as.character(pred$value))
pred=pred[order(pred$from),]
nodes <- data.frame(c("75-100%","50-75%","25-50%","0-25%","75-100%","50-75%","25-50%","0-25%"))
names(nodes) <- "name"

#Keeping only open-fished sites that are transformed into restricted fishing
Rest <- pred[which(pred$scenario != "reserve"),]
link <- Rest[,c(1:3)] #delete scenarii
link <- link[-c(1:4),] #initial pools not required
link$from <- link$from-1 # this package requires a starting point = 0
link$to <- link$to-1
link$zeroout=ifelse(link$value==0,1,0)
link=link[link$zeroout==0,]
link$zeroout=NULL
link$group=ifelse(link$from==0,"group1", ifelse(link$from==1,"group2", ifelse(link$from==2,"group3", "group4")))

# prepare color scale
my_color <- 'd3.scaleOrdinal() .domain(["group1", "group2","group3", "group4", "group5", "group6", "group7", "group8"]) 
.range(["#7a0177", "#c51b8a","#f768a1","black", "#7a0177", "#c51b8a","#f768a1","black"])'


nodes$group = as.factor(c("group1", "group2","group3", "group4", "group5", "group6", "group7", "group8"))
networkR_pred  <- sankeyNetwork(Links = link, Nodes = nodes, Source = "from",
                                Target = "to", Value = "value", NodeID = "name",
                                units = "", fontSize = 13, nodeWidth = 50, colourScale=my_color, NodeGroup = "group", LinkGroup = "group", width=300)
networkR_pred              

#Keeping only open-fished sites that are transformed into marine reserves
MR <- pred[which(pred$scenario != "restricted"),]
link <- MR[,c(1:3)]
link <- link[-c(1:4),]
link$from <- link$from-1
link$to <- link$to-1
link$zeroout=ifelse(link$value==0,1,0)
link=link[link$zeroout==0,]
link$zeroout=NULL
link$group=ifelse(link$from==0,"group1", ifelse(link$from==1,"group2", ifelse(link$from==2,"group3", "group4")))
nodes$sites=c(round(MR[1:4,3],digits=3),round(sum(fished_newdata$unfished_pred75),digits=3),round(sum(fished_newdata$unfished_pred50[fished_newdata$unfished_pred75==0]),digits=3),round(sum(fished_newdata$unfished_pred25[fished_newdata$unfished_pred50==0]),digits=3),round(length(fished_newdata$unfished_pred25[fished_newdata$unfished_pred25==0]),digits=3))
nodes$sites=c(round(MR[1:4,3]/length(fished_newdata$unfished_pred75),digits=3)*100,round(sum(fished_newdata$unfished_pred75)/length(fished_newdata$unfished_pred75),digits=3)*100,round(sum(fished_newdata$unfished_pred50[fished_newdata$unfished_pred75==0])/length(fished_newdata$unfished_pred75),digits=3)*100,round(sum(fished_newdata$unfished_pred25[fished_newdata$unfished_pred50==0])/length(fished_newdata$unfished_pred75),digits=3)*100,round(length(fished_newdata$unfished_pred25[fished_newdata$unfished_pred25==0])/length(fished_newdata$unfished_pred75),digits=3)*100)
nodes$sitesnolabel=rep("",length(nodes$sites))
networkMR_pred <- sankeyNetwork(Links = link, Nodes = nodes, Source = "from",
                                Target = "to", Value = "value", NodeID = "sites",
                                units = , fontSize = 14, nodeWidth = 50,colourScale=my_color, NodeGroup = "group", LinkGroup = "group", width=300)
networkMR_pred 



#do the same for trait diversity
fished_newdata2=fished_newdata[!is.na(fished_newdata$fundiv_newdata),]
fished_newdata2$predicted_func25=ifelse(fished_newdata2$fundiv_newdata<rfunc*0.25,0,1)
fished_newdata2$predicted_func50=ifelse(fished_newdata2$fundiv_newdata<rfunc*0.50,0,1)
fished_newdata2$predicted_func75=ifelse(fished_newdata2$fundiv_newdata<rfunc*0.75,0,1)
fished_newdata2$unfished_func25=ifelse(fished_newdata2$fundiv_newdata_reserves<rfunc*0.25,0,1)
fished_newdata2$unfished_func50=ifelse(fished_newdata2$fundiv_newdata_reserves<rfunc*0.50,0,1)
fished_newdata2$unfished_func75=ifelse(fished_newdata2$fundiv_newdata_reserves<rfunc*0.75,0,1)
fished_newdata2$restricted_func25=ifelse(fished_newdata2$fundiv_newdata_restricted<rfunc*0.25,0,1)
fished_newdata2$restricted_func50=ifelse(fished_newdata2$fundiv_newdata_restricted<rfunc*0.50,0,1)
fished_newdata2$restricted_func75=ifelse(fished_newdata2$fundiv_newdata_restricted<rfunc*0.75,0,1)


fished_newdata2$unfished_changefundiv=ifelse(fished_newdata2$predicted_func25==0 &fished_newdata2$unfished_func75==1, "0to75", 
                                             ifelse(fished_newdata2$predicted_func25==0 & fished_newdata2$unfished_func50==1 & fished_newdata2$unfished_func75==0, "0to50",
                                                    ifelse(fished_newdata2$predicted_func25==0 & fished_newdata2$unfished_func25==1 & fished_newdata2$unfished_func50==0, "0to25",
                                                           ifelse(fished_newdata2$predicted_func25==0 & fished_newdata2$unfished_func25==0, "0to0",
                                                                  ifelse(fished_newdata2$predicted_func25==1 &fished_newdata2$predicted_func50==0 &fished_newdata2$unfished_func75==1, "25to75",
                                                                         ifelse(fished_newdata2$predicted_func25==1 &fished_newdata2$predicted_func50==1&fished_newdata2$predicted_func75==0 &fished_newdata2$unfished_func75==1, "50to75",
                                                                                ifelse(fished_newdata2$predicted_func25==1 &fished_newdata2$predicted_func50==1&fished_newdata2$predicted_func75==1 &fished_newdata2$unfished_func75==1, "75to75",
                                                                                       ifelse(fished_newdata2$predicted_func25==1 &fished_newdata2$predicted_func50==0 &fished_newdata2$unfished_func75==0&fished_newdata2$unfished_func50==1, "25to50",
                                                                                              ifelse(fished_newdata2$predicted_func25==1 &fished_newdata2$predicted_func50==1 &fished_newdata2$unfished_func75==0&fished_newdata2$unfished_func50==1, "50to50",
                                                                                                     ifelse(fished_newdata2$predicted_func25==1 &fished_newdata2$predicted_func50==0 &fished_newdata2$unfished_func25==1 &fished_newdata2$unfished_func75==0&fished_newdata2$unfished_func50==0, "25to25",-999))))))))))
fished_newdata2$restricted_changefundiv=ifelse(fished_newdata2$predicted_func25==0 &fished_newdata2$restricted_func75==1, "0to75", 
                                               ifelse(fished_newdata2$predicted_func25==0 & fished_newdata2$restricted_func50==1 & fished_newdata2$restricted_func75==0, "0to50",
                                                      ifelse(fished_newdata2$predicted_func25==0 & fished_newdata2$restricted_func25==1 & fished_newdata2$restricted_func50==0, "0to25",
                                                             ifelse(fished_newdata2$predicted_func25==0 & fished_newdata2$restricted_func25==0, "0to0",
                                                                    ifelse(fished_newdata2$predicted_func25==1 &fished_newdata2$predicted_func50==0 &fished_newdata2$restricted_func75==1, "25to75",
                                                                           ifelse(fished_newdata2$predicted_func25==1 &fished_newdata2$predicted_func50==1&fished_newdata2$predicted_func75==0 &fished_newdata2$restricted_func75==1, "50to75",
                                                                                  ifelse(fished_newdata2$predicted_func25==1 &fished_newdata2$predicted_func50==1&fished_newdata2$predicted_func75==1 &fished_newdata2$restricted_func75==1, "75to75",
                                                                                         ifelse(fished_newdata2$predicted_func25==1 &fished_newdata2$predicted_func50==0 &fished_newdata2$restricted_func75==0&fished_newdata2$restricted_func50==1, "25to50",
                                                                                                ifelse(fished_newdata2$predicted_func25==1 &fished_newdata2$predicted_func50==1 &fished_newdata2$restricted_func75==0&fished_newdata2$restricted_func50==1, "50to50",
                                                                                                       ifelse(fished_newdata2$predicted_func25==1 &fished_newdata2$predicted_func50==0 &fished_newdata2$restricted_func25==1 &fished_newdata2$restricted_func75==0&fished_newdata2$restricted_func50==0, "25to25",-999))))))))))



matrixfunc=matrix(NA,nrow=24,ncol=5)
matrixfunc[1:4,1]=rep(0,4)
matrixfunc[1:4,2]=c(1:4)
sum(c(length(fished_newdata2$fundiv_newdata[fished_newdata2$fundiv_newdata>rfunc*0.75]),length(fished_newdata2$fundiv_newdata[fished_newdata2$fundiv_newdata>rfunc*0.50&fished_newdata2$fundiv_newdata<rfunc*0.75]),length(fished_newdata2$fundiv_newdata[fished_newdata2$fundiv_newdata>rfunc*0.25&fished_newdata2$fundiv_newdata<rfunc*0.50]),length(fished_newdata2$fundiv_newdata[fished_newdata2$fundiv_newdata<rfunc*0.25])))
matrixfunc[1:4,3]=c(length(fished_newdata2$fundiv_newdata[fished_newdata2$fundiv_newdata>rfunc*0.75]),length(fished_newdata2$fundiv_newdata[fished_newdata2$fundiv_newdata>rfunc*0.50&fished_newdata2$fundiv_newdata<rfunc*0.75]),length(fished_newdata2$fundiv_newdata[fished_newdata2$fundiv_newdata>rfunc*0.25&fished_newdata2$fundiv_newdata<rfunc*0.50]),length(fished_newdata2$fundiv_newdata[fished_newdata2$fundiv_newdata<rfunc*0.25]))
matrixfunc[1:4,4]=rep("Initial", 4)
matrixfunc[1:4,5]=c("75-100%","50-75%","25-50%","0-25%")
matrixfunc[5:8,1]=rep(4,4)
matrixfunc[5:8,2]=c(8,7,6,5)
matrixfunc[5:8,4]=rep("reserve",4)
matrixfunc[5:8,5]=c("0to0","0to25","0to50","0to75")
matrixfunc[5:8,3]=c(length(fished_newdata2$unfished_changefundiv[fished_newdata2$unfished_changefundiv=="0to0"]),length(fished_newdata2$unfished_changefundiv[fished_newdata2$unfished_changefundiv=="0to25"]),length(fished_newdata2$unfished_changefundiv[fished_newdata2$unfished_changefundiv=="0to50"]),length(fished_newdata2$unfished_changefundiv[fished_newdata2$unfished_changefundiv=="0to75"]))
matrixfunc[9:11,1]=rep(3,3)
matrixfunc[9:11,2]=c(7,6,5)
matrixfunc[9:11,4]=rep("reserve",3)
matrixfunc[9:11,5]=c("25to25","25to50","25to75")
matrixfunc[9:11,3]=c(length(fished_newdata2$unfished_changefundiv[fished_newdata2$unfished_changefundiv=="25to25"]),length(fished_newdata2$unfished_changefundiv[fished_newdata2$unfished_changefundiv=="25to50"]),length(fished_newdata2$unfished_changefundiv[fished_newdata2$unfished_changefundiv=="25to75"]))
matrixfunc[12:13,1]=rep(2,2)
matrixfunc[12:13,2]=c(6,5)
matrixfunc[12:13,4]=rep("reserve",2)
matrixfunc[12:13,5]=c("50to50","50to75")
matrixfunc[12:13,3]=c(length(fished_newdata2$unfished_changefundiv[fished_newdata2$unfished_changefundiv=="50to50"]),length(fished_newdata2$unfished_changefundiv[fished_newdata2$unfished_changefundiv=="50to75"]))
matrixfunc[14,1]=1
matrixfunc[14,2]=5
matrixfunc[14,4]="reserve"
matrixfunc[14,5]="75to75"
matrixfunc[14,3]=length(fished_newdata2$unfished_changefundiv[fished_newdata2$unfished_changefundiv=="75to75"])
matrixfunc[15:18,1]=rep(4,4)
matrixfunc[15:18,2]=c(8,7,6,5)
matrixfunc[15:18,4]=rep("restricted",4)
matrixfunc[15:18,5]=c("0to0","0to25","0to50","0to75")
matrixfunc[15:18,3]=c(length(fished_newdata2$restricted_changefundiv[fished_newdata2$restricted_changefundiv=="0to0"]),length(fished_newdata2$restricted_changefundiv[fished_newdata2$restricted_changefundiv=="0to25"]),length(fished_newdata2$restricted_changefundiv[fished_newdata2$restricted_changefundiv=="0to50"]),length(fished_newdata2$restricted_changefundiv[fished_newdata2$restricted_changefundiv=="0to75"]))
matrixfunc[19:21,1]=rep(3,3)
matrixfunc[19:21,2]=c(7,6,5)
matrixfunc[19:21,4]=rep("restricted",3)
matrixfunc[19:21,5]=c("25to25","25to50","25to75")
matrixfunc[19:21,3]=c(length(fished_newdata2$restricted_changefundiv[fished_newdata2$restricted_changefundiv=="25to25"]),length(fished_newdata2$restricted_changefundiv[fished_newdata2$restricted_changefundiv=="25to50"]),length(fished_newdata2$restricted_changefundiv[fished_newdata2$restricted_changefundiv=="25to75"]))
matrixfunc[22:23,1]=rep(2,2)
matrixfunc[22:23,2]=c(6,5)
matrixfunc[22:23,4]=rep("restricted",2)
matrixfunc[22:23,5]=c("50to50","50to75")
matrixfunc[22:23,3]=c(length(fished_newdata2$restricted_changefundiv[fished_newdata2$restricted_changefundiv=="50to50"]),length(fished_newdata2$restricted_changefundiv[fished_newdata2$restricted_changefundiv=="50to75"]))
matrixfunc[24,1]=1
matrixfunc[24,2]=5
matrixfunc[24,4]="restricted"
matrixfunc[24,5]="75to75"
matrixfunc[24,3]=length(fished_newdata2$restricted_changefundiv[fished_newdata2$restricted_changefundiv=="75to75"])

#alluvial charts
fundiv=as.data.frame(matrixfunc)
colnames(fundiv)=c("from","to","value","scenario","clarified")
fundiv$from=as.numeric(as.character(fundiv$from))
fundiv$to=as.numeric(as.character(fundiv$to))
fundiv$value=as.numeric(as.character(fundiv$value))
fundiv=fundiv[order(fundiv$from),]
nodes <- data.frame(c("75-100%","50-75%","25-50%","0-25%","75-100%","50-75%","25-50%","0-25%"))
names(nodes) <- "name"

#Keeping only open-fished sites that are transformed into restricted fishing
Rest <- fundiv[which(fundiv$scenario != "reserve"),]
link <- Rest[,c(1:3)] #delete scenarii
link <- link[-c(1:4),] #initial pools not required
link$from <- link$from-1 # this package requires a starting point = 0
link$to <- link$to-1
link$zeroout=ifelse(link$value==0,1,0)
link=link[link$zeroout==0,]
link$zeroout=NULL
link$group=ifelse(link$from==0,"group1", ifelse(link$from==1,"group2", ifelse(link$from==2,"group3", "group4")))
nodes$group = as.factor(c("group1", "group2","group3", "group4", "group5", "group6", "group7", "group8"))

networkR_func  <- sankeyNetwork(Links = link, Nodes = nodes, Source = "from",
                                Target = "to", Value = "value", NodeID = "name",
                                units = "", fontSize = 13, nodeWidth = 50, colourScale=my_color, NodeGroup = "group", LinkGroup = "group", width=300)

networkR_func              

#Keeping only open-fished sites that are transformed into marine reserves
MR <- fundiv[which( fundiv$scenario != "restricted"),]
link <- MR[,c(1:3)]
link <- link[-c(1:4),]
link$from <- link$from-1
link$to <- link$to-1
link$zeroout=ifelse(link$value==0,1,0)
link=link[link$zeroout==0,]
link$zeroout=NULL
link$group=ifelse(link$from==0,"group1", ifelse(link$from==1,"group2", ifelse(link$from==2,"group3", "group4")))
nodes$sites=c(round(MR[1:4,3],digits=3),round(sum(fished_newdata2$unfished_func75),digits=3),round(sum(fished_newdata2$unfished_func50[fished_newdata2$unfished_func75==0]),digits=3),round(sum(fished_newdata2$unfished_func25[fished_newdata2$unfished_func50==0]),digits=3),round(length(fished_newdata2$unfished_func25[fished_newdata2$unfished_func25==0]),digits=3))
nodes$sites=c(round(MR[1:4,3]/length(fished_newdata2$unfished_func75),digits=3)*100,round(sum(fished_newdata2$unfished_func75)/length(fished_newdata2$unfished_func75),digits=3)*100,round(sum(fished_newdata2$unfished_func50[fished_newdata2$unfished_func75==0])/length(fished_newdata2$unfished_func75),digits=3)*100,round(sum(fished_newdata2$unfished_func25[fished_newdata2$unfished_func50==0])/length(fished_newdata2$unfished_func75),digits=3)*100,round(length(fished_newdata2$unfished_func25[fished_newdata2$unfished_func25==0])/length(fished_newdata2$unfished_func75),digits=3)*100)
nodes$sitesnolabel=rep("",length(nodes$sites))
networkMR_func <- sankeyNetwork(Links = link, Nodes = nodes, Source = "from",
                                Target = "to", Value = "value", NodeID = "sites",
                                units = "", fontSize = 13, nodeWidth = 50,colourScale=my_color, NodeGroup = "group", LinkGroup = "group", width=300)


networkMR_func 


#parrotfsh scraping
fished_newdata3=fished_newdata2[!is.na(fished_newdata2$her_newdata),]
fished_newdata3$predicted_herb25=ifelse(fished_newdata3$her_newdata<rherb*0.25,0,1)
fished_newdata3$predicted_herb50=ifelse(fished_newdata3$her_newdata<rherb*0.50,0,1)
fished_newdata3$predicted_herb75=ifelse(fished_newdata3$her_newdata<rherb*0.75,0,1)
fished_newdata3$unfished_herb25=ifelse(fished_newdata3$her_newdata_reserves<rherb*0.25,0,1)
fished_newdata3$unfished_herb50=ifelse(fished_newdata3$her_newdata_reserves<rherb*0.50,0,1)
fished_newdata3$unfished_herb75=ifelse(fished_newdata3$her_newdata_reserves<rherb*0.75,0,1)
fished_newdata3$restricted_herb25=ifelse(fished_newdata3$her_newdata_restricted<rherb*0.25,0,1)
fished_newdata3$restricted_herb50=ifelse(fished_newdata3$her_newdata_restricted<rherb*0.50,0,1)
fished_newdata3$restricted_herb75=ifelse(fished_newdata3$her_newdata_restricted<rherb*0.75,0,1)

fished_newdata3$unfished_changeherb=ifelse(fished_newdata3$predicted_herb25==0 &fished_newdata3$unfished_herb75==1, "0to75", 
                                           ifelse(fished_newdata3$predicted_herb25==0 & fished_newdata3$unfished_herb50==1 & fished_newdata3$unfished_herb75==0, "0to50",
                                                  ifelse(fished_newdata3$predicted_herb25==0 & fished_newdata3$unfished_herb25==1 & fished_newdata3$unfished_herb50==0, "0to25",
                                                         ifelse(fished_newdata3$predicted_herb25==0 & fished_newdata3$unfished_herb25==0, "0to0",
                                                                ifelse(fished_newdata3$predicted_herb25==1 &fished_newdata3$predicted_herb50==0 &fished_newdata3$unfished_herb75==1, "25to75",
                                                                       ifelse(fished_newdata3$predicted_herb25==1 &fished_newdata3$predicted_herb50==1&fished_newdata3$predicted_herb75==0 &fished_newdata3$unfished_herb75==1, "50to75",
                                                                              ifelse(fished_newdata3$predicted_herb25==1 &fished_newdata3$predicted_herb50==1&fished_newdata3$predicted_herb75==1 &fished_newdata3$unfished_herb75==1, "75to75",
                                                                                     ifelse(fished_newdata3$predicted_herb25==1 &fished_newdata3$predicted_herb50==0 &fished_newdata3$unfished_herb75==0&fished_newdata3$unfished_herb50==1, "25to50",
                                                                                            ifelse(fished_newdata3$predicted_herb25==1 &fished_newdata3$predicted_herb50==1 &fished_newdata3$unfished_herb75==0&fished_newdata3$unfished_herb50==1, "50to50",
                                                                                                   ifelse(fished_newdata3$predicted_herb25==1 &fished_newdata3$predicted_herb50==0 &fished_newdata3$unfished_herb25==1 &fished_newdata3$unfished_herb75==0&fished_newdata3$unfished_herb50==0, "25to25",-999))))))))))
fished_newdata3$restricted_changeherb=ifelse(fished_newdata3$predicted_herb25==0 &fished_newdata3$restricted_herb75==1, "0to75", 
                                             ifelse(fished_newdata3$predicted_herb25==0 & fished_newdata3$restricted_herb50==1 & fished_newdata3$restricted_herb75==0, "0to50",
                                                    ifelse(fished_newdata3$predicted_herb25==0 & fished_newdata3$restricted_herb25==1 & fished_newdata3$restricted_herb50==0, "0to25",
                                                           ifelse(fished_newdata3$predicted_herb25==0 & fished_newdata3$restricted_herb25==0, "0to0",
                                                                  ifelse(fished_newdata3$predicted_herb25==1 &fished_newdata3$predicted_herb50==0 &fished_newdata3$restricted_herb75==1, "25to75",
                                                                         ifelse(fished_newdata3$predicted_herb25==1 &fished_newdata3$predicted_herb50==1&fished_newdata3$predicted_herb75==0 &fished_newdata3$restricted_herb75==1, "50to75",
                                                                                ifelse(fished_newdata3$predicted_herb25==1 &fished_newdata3$predicted_herb50==1&fished_newdata3$predicted_herb75==1 &fished_newdata3$restricted_herb75==1, "75to75",
                                                                                       ifelse(fished_newdata3$predicted_herb25==1 &fished_newdata3$predicted_herb50==0 &fished_newdata3$restricted_herb75==0&fished_newdata3$restricted_herb50==1, "25to50",
                                                                                              ifelse(fished_newdata3$predicted_herb25==1 &fished_newdata3$predicted_herb50==1 &fished_newdata3$restricted_herb75==0&fished_newdata3$restricted_herb50==1, "50to50",
                                                                                                     ifelse(fished_newdata3$predicted_herb25==1 &fished_newdata3$predicted_herb50==0 &fished_newdata3$restricted_herb25==1 &fished_newdata3$restricted_herb75==0&fished_newdata3$restricted_herb50==0, "25to25",-999))))))))))



matrixherb=matrix(NA,nrow=24,ncol=5)
matrixherb[1:4,1]=rep(0,4)
matrixherb[1:4,2]=c(1:4)
matrixherb[1:4,3]=c(length(fished_newdata3$her_newdata[fished_newdata3$her_newdata>rherb*0.75]),length(fished_newdata3$her_newdata[fished_newdata3$her_newdata>rherb*0.50&fished_newdata3$her_newdata<rherb*0.75]),length(fished_newdata3$her_newdata[fished_newdata3$her_newdata>rherb*0.25&fished_newdata3$her_newdata<rherb*0.50]),length(fished_newdata3$her_newdata[fished_newdata3$her_newdata<rherb*0.25]))
matrixherb[1:4,4]=rep("Initial", 4)
matrixherb[1:4,5]=c("75-100%","50-75%","25-50%","0-25%")
matrixherb[5:8,1]=rep(4,4)
matrixherb[5:8,2]=c(8,7,6,5)
matrixherb[5:8,4]=rep("reserve",4)
matrixherb[5:8,5]=c("0to0","0to25","0to50","0to75")
matrixherb[5:8,3]=c(length(fished_newdata3$unfished_changeherb[fished_newdata3$unfished_changeherb=="0to0"]),length(fished_newdata3$unfished_changeherb[fished_newdata3$unfished_changeherb=="0to25"]),length(fished_newdata3$unfished_changeherb[fished_newdata3$unfished_changeherb=="0to50"]),length(fished_newdata3$unfished_changeherb[fished_newdata3$unfished_changeherb=="0to75"]))
matrixherb[9:11,1]=rep(3,3)
matrixherb[9:11,2]=c(7,6,5)
matrixherb[9:11,4]=rep("reserve",3)
matrixherb[9:11,5]=c("25to25","25to50","25to75")
matrixherb[9:11,3]=c(length(fished_newdata3$unfished_changeherb[fished_newdata3$unfished_changeherb=="25to25"]),length(fished_newdata3$unfished_changeherb[fished_newdata3$unfished_changeherb=="25to50"]),length(fished_newdata3$unfished_changeherb[fished_newdata3$unfished_changeherb=="25to75"]))
matrixherb[12:13,1]=rep(2,2)
matrixherb[12:13,2]=c(6,5)
matrixherb[12:13,4]=rep("reserve",2)
matrixherb[12:13,5]=c("50to50","50to75")
matrixherb[12:13,3]=c(length(fished_newdata3$unfished_changeherb[fished_newdata3$unfished_changeherb=="50to50"]),length(fished_newdata3$unfished_changeherb[fished_newdata3$unfished_changeherb=="50to75"]))
matrixherb[14,1]=1
matrixherb[14,2]=5
matrixherb[14,4]="reserve"
matrixherb[14,5]="75to75"
matrixherb[14,3]=length(fished_newdata3$unfished_changeherb[fished_newdata3$unfished_changeherb=="75to75"])
matrixherb[15:18,1]=rep(4,4)
matrixherb[15:18,2]=c(8,7,6,5)
matrixherb[15:18,4]=rep("restricted",4)
matrixherb[15:18,5]=c("0to0","0to25","0to50","0to75")
matrixherb[15:18,3]=c(length(fished_newdata3$restricted_changeherb[fished_newdata3$restricted_changeherb=="0to0"]),length(fished_newdata3$restricted_changeherb[fished_newdata3$restricted_changeherb=="0to25"]),length(fished_newdata3$restricted_changeherb[fished_newdata3$restricted_changeherb=="0to50"]),length(fished_newdata3$restricted_changeherb[fished_newdata3$restricted_changeherb=="0to75"]))
matrixherb[19:21,1]=rep(3,3)
matrixherb[19:21,2]=c(7,6,5)
matrixherb[19:21,4]=rep("restricted",3)
matrixherb[19:21,5]=c("25to25","25to50","25to75")
matrixherb[19:21,3]=c(length(fished_newdata3$restricted_changeherb[fished_newdata3$restricted_changeherb=="25to25"]),length(fished_newdata3$restricted_changeherb[fished_newdata3$restricted_changeherb=="25to50"]),length(fished_newdata3$restricted_changeherb[fished_newdata3$restricted_changeherb=="25to75"]))
matrixherb[22:23,1]=rep(2,2)
matrixherb[22:23,2]=c(6,5)
matrixherb[22:23,4]=rep("restricted",2)
matrixherb[22:23,5]=c("50to50","50to75")
matrixherb[22:23,3]=c(length(fished_newdata3$restricted_changeherb[fished_newdata3$restricted_changeherb=="50to50"]),length(fished_newdata3$restricted_changeherb[fished_newdata3$restricted_changeherb=="50to75"]))
matrixherb[24,1]=1
matrixherb[24,2]=5
matrixherb[24,4]="restricted"
matrixherb[24,5]="75to75"
matrixherb[24,3]=length(fished_newdata3$restricted_changeherb[fished_newdata3$restricted_changeherb=="75to75"])

#alluvial charts
herb=as.data.frame(matrixherb)
colnames(herb)=c("from","to","value","scenario","clarified")
herb$from=as.numeric(as.character(herb$from))
herb$to=as.numeric(as.character(herb$to))
herb$value=as.numeric(as.character(herb$value))
herb=herb[order(herb$from),]
nodes <- data.frame(c("75-100%","50-75%","25-50%","0-25%","75-100%","50-75%","25-50%","0-25%"))
names(nodes) <- "name"

#Keeping only open-fished sites that are transformed into restricted fishing
Rest <- herb[which(herb$scenario != "reserve"),]
link <- Rest[,c(1:3)] #delete scenarii
link <- link[-c(1:4),] #initial pools not required
link$from <- link$from-1 # this package requires a starting point = 0
link$to <- link$to-1
link$zeroout=ifelse(link$value==0,1,0)
link=link[link$zeroout==0,]
link$zeroout=NULL
link$group=ifelse(link$from==0,"group1", ifelse(link$from==1,"group2", ifelse(link$from==2,"group3", "group4")))
nodes$group = as.factor(c("group1", "group2","group3", "group4", "group5", "group6", "group7", "group8"))

networkR_herb  <- sankeyNetwork(Links = link, Nodes = nodes, Source = "from",
                                Target = "to", Value = "value", NodeID = "name",
                                units = "", fontSize = 13, nodeWidth = 50, colourScale=my_color, NodeGroup = "group", LinkGroup = "group", width=300)

networkR_herb              

#Keeping only open-fished sites that are transformed into marine reserves
MR <- herb[which( herb$scenario != "restricted"),]
link <- MR[,c(1:3)]
link <- link[-c(1:4),]
link$from <- link$from-1
link$to <- link$to-1
link$zeroout=ifelse(link$value==0,1,0)
link=link[link$zeroout==0,]
link$zeroout=NULL
link$group=ifelse(link$from==0,"group1", ifelse(link$from==1,"group2", ifelse(link$from==2,"group3", "group4")))
nodes$sites=c(round(MR[1:4,3],digits=3),round(sum(fished_newdata3$unfished_herb75),digits=3),round(sum(fished_newdata3$unfished_herb50[fished_newdata3$unfished_herb75==0]),digits=3),round(sum(fished_newdata3$unfished_herb25[fished_newdata3$unfished_herb50==0]),digits=3),round(length(fished_newdata3$unfished_herb25[fished_newdata3$unfished_herb25==0]),digits=3))
nodes$sites=c(round(MR[1:4,3]/length(fished_newdata3$unfished_herb75),digits=3)*100,round(sum(fished_newdata3$unfished_herb75)/length(fished_newdata3$unfished_herb75),digits=3)*100,round(sum(fished_newdata3$unfished_herb50[fished_newdata3$unfished_herb75==0])/length(fished_newdata3$unfished_herb75),digits=3)*100,round(sum(fished_newdata3$unfished_herb25[fished_newdata3$unfished_herb50==0])/length(fished_newdata3$unfished_herb75),digits=3)*100,round(length(fished_newdata3$unfished_herb25[fished_newdata3$unfished_herb25==0])/length(fished_newdata3$unfished_herb75),digits=3)*100)
nodes$sitesnolabel=rep("",length(nodes$sites))
networkMR_herb <- sankeyNetwork(Links = link, Nodes = nodes, Source = "from",
                                Target = "to", Value = "value", NodeID = "sites",
                                units = "", fontSize = 13, nodeWidth = 50,colourScale=my_color, NodeGroup = "group", LinkGroup = "group", width=300)


networkMR_herb 



#all three goals
fished_newdata3$predicted_tot25=ifelse(fished_newdata3$predicted_pred25==1&fished_newdata3$predicted_func25==1&fished_newdata3$predicted_herb25==1,1,0)
fished_newdata3$predicted_tot50=ifelse(fished_newdata3$predicted_pred50==1&fished_newdata3$predicted_func50==1&fished_newdata3$predicted_herb50==1,1,0)
fished_newdata3$predicted_tot75=ifelse(fished_newdata3$predicted_pred75==1&fished_newdata3$predicted_func75==1&fished_newdata3$predicted_herb75==1,1,0)
fished_newdata3$unfished_tot25=ifelse(fished_newdata3$unfished_pred25==1&fished_newdata3$unfished_func25==1&fished_newdata3$unfished_herb25==1,1,0)
fished_newdata3$unfished_tot50=ifelse(fished_newdata3$unfished_pred50==1&fished_newdata3$unfished_func50==1&fished_newdata3$unfished_herb50==1,1,0)
fished_newdata3$unfished_tot75=ifelse(fished_newdata3$unfished_pred75==1&fished_newdata3$unfished_func75==1&fished_newdata3$unfished_herb75==1,1,0)
fished_newdata3$restricted_tot25=ifelse(fished_newdata3$restricted_pred25==1&fished_newdata3$restricted_func25==1&fished_newdata3$restricted_herb25==1,1,0)
fished_newdata3$restricted_tot50=ifelse(fished_newdata3$restricted_pred50==1&fished_newdata3$restricted_func50==1&fished_newdata3$restricted_herb50==1,1,0)
fished_newdata3$restricted_tot75=ifelse(fished_newdata3$restricted_pred75==1&fished_newdata3$restricted_func75==1&fished_newdata3$restricted_herb75==1,1,0)

fished_newdata3$unfished_changetot=ifelse(fished_newdata3$predicted_tot25==0 &fished_newdata3$unfished_tot75==1, "0to75", 
                                          ifelse(fished_newdata3$predicted_tot25==0 & fished_newdata3$unfished_tot50==1 & fished_newdata3$unfished_tot75==0, "0to50",
                                                 ifelse(fished_newdata3$predicted_tot25==0 & fished_newdata3$unfished_tot25==1 & fished_newdata3$unfished_tot50==0, "0to25",
                                                        ifelse(fished_newdata3$predicted_tot25==0 & fished_newdata3$unfished_tot25==0, "0to0",
                                                               ifelse(fished_newdata3$predicted_tot25==1 &fished_newdata3$predicted_tot50==0 &fished_newdata3$unfished_tot75==1, "25to75",
                                                                      ifelse(fished_newdata3$predicted_tot25==1 &fished_newdata3$predicted_tot50==1&fished_newdata3$predicted_tot75==0 &fished_newdata3$unfished_tot75==1, "50to75",
                                                                             ifelse(fished_newdata3$predicted_tot25==1 &fished_newdata3$predicted_tot50==1&fished_newdata3$predicted_tot75==1 &fished_newdata3$unfished_tot75==1, "75to75",
                                                                                    ifelse(fished_newdata3$predicted_tot25==1 &fished_newdata3$predicted_tot50==0 &fished_newdata3$unfished_tot75==0&fished_newdata3$unfished_tot50==1, "25to50",
                                                                                           ifelse(fished_newdata3$predicted_tot25==1 &fished_newdata3$predicted_tot50==1 &fished_newdata3$unfished_tot75==0&fished_newdata3$unfished_tot50==1, "50to50",
                                                                                                  ifelse(fished_newdata3$predicted_tot25==1 &fished_newdata3$predicted_tot50==0 &fished_newdata3$unfished_tot25==1 &fished_newdata3$unfished_tot75==0&fished_newdata3$unfished_tot50==0, "25to25",-999))))))))))
fished_newdata3$restricted_changetot=ifelse(fished_newdata3$predicted_tot25==0 &fished_newdata3$restricted_tot75==1, "0to75", 
                                            ifelse(fished_newdata3$predicted_tot25==0 & fished_newdata3$restricted_tot50==1 & fished_newdata3$restricted_tot75==0, "0to50",
                                                   ifelse(fished_newdata3$predicted_tot25==0 & fished_newdata3$restricted_tot25==1 & fished_newdata3$restricted_tot50==0, "0to25",
                                                          ifelse(fished_newdata3$predicted_tot25==0 & fished_newdata3$restricted_tot25==0, "0to0",
                                                                 ifelse(fished_newdata3$predicted_tot25==1 &fished_newdata3$predicted_tot50==0 &fished_newdata3$restricted_tot75==1, "25to75",
                                                                        ifelse(fished_newdata3$predicted_tot25==1 &fished_newdata3$predicted_tot50==1&fished_newdata3$predicted_tot75==0 &fished_newdata3$restricted_tot75==1, "50to75",
                                                                               ifelse(fished_newdata3$predicted_tot25==1 &fished_newdata3$predicted_tot50==1&fished_newdata3$predicted_tot75==1 &fished_newdata3$restricted_tot75==1, "75to75",
                                                                                      ifelse(fished_newdata3$predicted_tot25==1 &fished_newdata3$predicted_tot50==0 &fished_newdata3$restricted_tot75==0&fished_newdata3$restricted_tot50==1, "25to50",
                                                                                             ifelse(fished_newdata3$predicted_tot25==1 &fished_newdata3$predicted_tot50==1 &fished_newdata3$restricted_tot75==0&fished_newdata3$restricted_tot50==1, "50to50",
                                                                                                    ifelse(fished_newdata3$predicted_tot25==1 &fished_newdata3$predicted_tot50==0 &fished_newdata3$restricted_tot25==1 &fished_newdata3$restricted_tot75==0&fished_newdata3$restricted_tot50==0, "25to25",-999))))))))))


matrixtot=matrix(NA,nrow=24,ncol=5)
matrixtot[1:4,1]=rep(0,4)
matrixtot[1:4,2]=c(1:4)
matrixtot[1:4,3]=c(length(fished_newdata3$predicted_tot25[fished_newdata3$predicted_tot25==0]),sum(fished_newdata3$predicted_tot25)-sum(fished_newdata3$predicted_tot50),sum(fished_newdata3$predicted_tot50)-sum(fished_newdata3$predicted_tot75),sum(fished_newdata3$predicted_tot75))
matrixtot[1:4,4]=rep("Initial", 4)
matrixtot[1:4,5]=c("75-100%","50-75%","25-50%","0-25%")
matrixtot[5:8,1]=rep(4,4)
matrixtot[5:8,2]=c(8,7,6,5)
matrixtot[5:8,4]=rep("reserve",4)
matrixtot[5:8,5]=c("0to0","0to25","0to50","0to75")
matrixtot[5:8,3]=c(length(fished_newdata3$unfished_changetot[fished_newdata3$unfished_changetot=="0to0"]),length(fished_newdata3$unfished_changetot[fished_newdata3$unfished_changetot=="0to25"]),length(fished_newdata3$unfished_changetot[fished_newdata3$unfished_changetot=="0to50"]),length(fished_newdata3$unfished_changetot[fished_newdata3$unfished_changetot=="0to75"]))
matrixtot[9:11,1]=rep(3,3)
matrixtot[9:11,2]=c(7,6,5)
matrixtot[9:11,4]=rep("reserve",3)
matrixtot[9:11,5]=c("25to25","25to50","25to75")
matrixtot[9:11,3]=c(length(fished_newdata3$unfished_changetot[fished_newdata3$unfished_changetot=="25to25"]),length(fished_newdata3$unfished_changetot[fished_newdata3$unfished_changetot=="25to50"]),length(fished_newdata3$unfished_changetot[fished_newdata3$unfished_changetot=="25to75"]))
matrixtot[12:13,1]=rep(2,2)
matrixtot[12:13,2]=c(6,5)
matrixtot[12:13,4]=rep("reserve",2)
matrixtot[12:13,5]=c("50to50","50to75")
matrixtot[12:13,3]=c(length(fished_newdata3$unfished_changetot[fished_newdata3$unfished_changetot=="50to50"]),length(fished_newdata3$unfished_changetot[fished_newdata3$unfished_changetot=="50to75"]))
matrixtot[14,1]=1
matrixtot[14,2]=5
matrixtot[14,4]="reserve"
matrixtot[14,5]="75to75"
matrixtot[14,3]=length(fished_newdata3$unfished_changetot[fished_newdata3$unfished_changetot=="75to75"])
matrixtot[15:18,1]=rep(4,4)
matrixtot[15:18,2]=c(8,7,6,5)
matrixtot[15:18,4]=rep("restricted",4)
matrixtot[15:18,5]=c("0to0","0to25","0to50","0to75")
matrixtot[15:18,3]=c(length(fished_newdata3$restricted_changetot[fished_newdata3$restricted_changetot=="0to0"]),length(fished_newdata3$restricted_changetot[fished_newdata3$restricted_changetot=="0to25"]),length(fished_newdata3$restricted_changetot[fished_newdata3$restricted_changetot=="0to50"]),length(fished_newdata3$restricted_changetot[fished_newdata3$restricted_changetot=="0to75"]))
matrixtot[19:21,1]=rep(3,3)
matrixtot[19:21,2]=c(7,6,5)
matrixtot[19:21,4]=rep("restricted",3)
matrixtot[19:21,5]=c("25to25","25to50","25to75")
matrixtot[19:21,3]=c(length(fished_newdata3$restricted_changetot[fished_newdata3$restricted_changetot=="25to25"]),length(fished_newdata3$restricted_changetot[fished_newdata3$restricted_changetot=="25to50"]),length(fished_newdata3$restricted_changetot[fished_newdata3$restricted_changetot=="25to75"]))
matrixtot[22:23,1]=rep(2,2)
matrixtot[22:23,2]=c(6,5)
matrixtot[22:23,4]=rep("restricted",2)
matrixtot[22:23,5]=c("50to50","50to75")
matrixtot[22:23,3]=c(length(fished_newdata3$restricted_changetot[fished_newdata3$restricted_changetot=="50to50"]),length(fished_newdata3$restricted_changetot[fished_newdata3$restricted_changetot=="50to75"]))
matrixtot[24,1]=1
matrixtot[24,2]=5
matrixtot[24,4]="restricted"
matrixtot[24,5]="75to75"
matrixtot[24,3]=length(fished_newdata3$restricted_changetot[fished_newdata3$restricted_changetot=="75to75"])

#alluvial charts
tot=as.data.frame(matrixtot)
colnames(tot)=c("from","to","value","scenario","clarified")
tot$from=as.numeric(as.character(tot$from))
tot$to=as.numeric(as.character(tot$to))
tot$value=as.numeric(as.character(tot$value))
tot=tot[order(tot$from),]
nodes <- data.frame(c("75-100%","50-75%","25-50%","0-25%","75-100%","50-75%","25-50%","0-25%"))
names(nodes) <- "name"

#Keeping only open-fished sites that are transformed into restricted fishing
Rest <- tot[which(tot$scenario != "reserve"),]
link <- Rest[,c(1:3)] #delete scenarii
link <- link[-c(1:4),] #initial pools not required
link$from <- link$from-1 # this package requires a starting point = 0
link$to <- link$to-1
link$zeroout=ifelse(link$value==0,1,0)
link=link[link$zeroout==0,]
link$zeroout=NULL
link$group=ifelse(link$from==0,"group1", ifelse(link$from==1,"group2", ifelse(link$from==2,"group3", "group4")))
nodes$group = as.factor(c("group1", "group2","group3", "group4", "group5", "group6", "group7", "group8"))

networkR_tot  <- sankeyNetwork(Links = link, Nodes = nodes, Source = "from",
                               Target = "to", Value = "value", NodeID = "name",
                               units = "", fontSize = 13, nodeWidth = 50, colourScale=my_color, NodeGroup = "group", LinkGroup = "group", width=300)

networkR_tot              

#Keeping only open-fished sites that are transformed into marine reserves
MR <- tot[which( tot$scenario != "restricted"),]
link <- MR[,c(1:3)]
link <- link[-c(1:4),]
link$from <- link$from-1
link$to <- link$to-1
link$zeroout=ifelse(link$value==0,1,0)
link=link[link$zeroout==0,]
link$zeroout=NULL
link$group=ifelse(link$from==0,"group1", ifelse(link$from==1,"group2", ifelse(link$from==2,"group3", "group4")))
nodes$sites=c(round(MR[1:4,3],digits=3),round(sum(fished_newdata3$unfished_tot75),digits=3),round(sum(fished_newdata3$unfished_tot50[fished_newdata3$unfished_tot75==0]),digits=3),round(sum(fished_newdata3$unfished_tot25[fished_newdata3$unfished_tot50==0]),digits=3),round(length(fished_newdata3$unfished_tot25[fished_newdata3$unfished_tot25==0]),digits=3))
nodes$sites=c(rev(round(MR[1:4,3]/length(fished_newdata3$unfished_tot75),digits=3)*100),round(sum(fished_newdata3$unfished_tot75)/length(fished_newdata3$unfished_tot75),digits=3)*100,round(sum(fished_newdata3$unfished_tot50[fished_newdata3$unfished_tot75==0])/length(fished_newdata3$unfished_tot75),digits=3)*100,round(sum(fished_newdata3$unfished_tot25[fished_newdata3$unfished_tot50==0])/length(fished_newdata3$unfished_tot75),digits=3)*100,round(length(fished_newdata3$unfished_tot25[fished_newdata3$unfished_tot25==0])/length(fished_newdata3$unfished_tot75),digits=3)*100)
nodes$sitesnolabel=rep("",length(nodes$sites))


networkMR_tot <- sankeyNetwork(Links = link, Nodes = nodes, Source = "from",
                               Target = "to", Value = "value", NodeID = "sites",
                               units = "", fontSize = 13, nodeWidth = 50,colourScale=my_color, NodeGroup = "group", LinkGroup = "group", width=300)


networkMR_tot 

##################################################################################
#run analyses with more stringent MPA definition
mpa_restricted=alldata[alldata$keep==1,]

her_mpa_restricted<-brm(Scraping_potential~
                          DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot2+
                          sRegional_population_growth+sOcean_prod+sClimate_stress+
                          sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                          (1|Larger/ReefCluster),
                        data=mpa_restricted[!is.na(mpa_restricted$Scraping_potential),],family=hurdle_lognormal(link = "identity", link_sigma = "log",
                                                                                                                 link_hu = "logit"),
                        iter=10000,  warmup=9000,
                        chains=4) 

fundiv_mpa_restricted<-brm(tTrait_diversity~
                             DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot2+
                             sRegional_population_growth+sOcean_prod+sClimate_stress+
                             sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                             (1|Larger/ReefCluster),
                           data=mpa_restricted[!is.na(mpa_restricted$tTrait_diversity),],family=gaussian,
                           iter=10000,  warmup=9000,
                           chains=4) 

B20_mpa_restricted<-brm(tBiomass_above20cm~
                          DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot2+
                          sRegional_population_growth+sOcean_prod+sClimate_stress+
                          sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                          (1|Larger/ReefCluster),
                        data=mpa_restricted,family=gaussian,
                        iter=10000,  warmup=9000,
                        chains=4) 

#for each row, get the probability of the column
prob25 <- function(x){length(x[(exp(x)-1)>rpred*0.25])/length(x)} 
prob50 <- function(x){length(x[(exp(x)-1)>rpred*0.50])/length(x)} 
prob75 <- function(x){length(x[(exp(x)-1)>rpred*0.75])/length(x)} 

#Create modelled data
MyData_b20<-expand.grid(sgrav_tot2=seq(range(alldata$sgrav_tot2)[[1]],range(alldata$sgrav_tot2)[[2]],length=101),
                        Protection=c("Fished","Restricted","UnfishedHigh"),
                        DepthCategory=c(">10m"," 0-4m","4-10m"),
                        CleanHabitat=c("Crest","Flat","Lagoon_Back reef","Slope"),
                        sTotal_sampling_area=0,
                        CensusMethod=c("Distance sampling","Point intercept","Standard belt transect"),
                        sRegional_population_growth=0,
                        sOcean_prod=0,
                        sHDI=0,
                        sReef_fish_landings_per_km2=0,
                        sClimate_stress=0,
                        sLarger_pop_size=0)

X <- model.matrix(~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot2+
                    sRegional_population_growth+sOcean_prod+sClimate_stress+
                    sHDI+sLarger_pop_size+sReef_fish_landings_per_km2,data=MyData_b20)

coefs_b20=as.matrix(B20_mpa_restricted)[, paste("b",row.names(fixef(B20_mpa_restricted)),sep="_")] #alternative

fit_b20=coefs_b20 %*% t(X)
dim(fit_b20)

#for each row, get the probability of the column
try25=as.data.frame(apply(fit_b20, 2, prob25) )
colnames(try25)="estimate25"
try50=as.data.frame(apply(fit_b20, 2, prob50) )
colnames(try50)="estimate50"
try75=as.data.frame(apply(fit_b20, 2, prob75) )
colnames(try75)="estimate75"

MyData_b20=MyData_b20 %>% cbind(tidyMCMC(as.mcmc(fit_b20),
                                         conf.int=T, conf.method='HPDinterval'))
MyData_b20=cbind(MyData_b20,try25,try50,try75)

# Create Delta
Reserves_conditioned_b20<- MyData_b20[which(MyData_b20$CleanHabitat=="Slope" & MyData_b20$Protection=="UnfishedHigh" & MyData_b20$CensusMethod=="Standard belt transect" &
                                              MyData_b20$DepthCategory=="4-10m"),]

#for each row, get the probability of the column
prob25 <- function(x){length(x[(exp(x))>rherb*0.25])/length(x)} 
prob50 <- function(x){length(x[(exp(x))>rherb*0.50])/length(x)} 
prob75 <- function(x){length(x[(exp(x))>rherb*0.75])/length(x)} 

MyData_hf<-expand.grid(sgrav_tot2=seq(range(alldata$sgrav_tot2)[[1]],range(alldata$sgrav_tot2)[[2]],length=101),
                       Protection=c("Fished","Restricted","UnfishedHigh"),
                       DepthCategory=c(">10m"," 0-4m","4-10m"),
                       CleanHabitat=c("Crest","Flat","Lagoon_Back reef","Slope"),
                       sTotal_sampling_area=0,
                       CensusMethod=c("Distance sampling","Point intercept","Standard belt transect"),
                       sRegional_population_growth=0,
                       sOcean_prod=0,
                       sHDI=0,
                       sReef_fish_landings_per_km2=0,
                       sClimate_stress=0,
                       sLarger_pop_size=0)

X <- model.matrix(~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot2+
                    sRegional_population_growth+sOcean_prod+sClimate_stress+
                    sHDI+sLarger_pop_size+sReef_fish_landings_per_km2,data=MyData_hf)

coefs_hf=as.matrix(her_mpa_restricted)[, paste("b",row.names(fixef(her_mpa_restricted)),sep="_")] #alternative

fit_hf=coefs_hf %*% t(X)
dim(fit_hf)
summary(fit_hf)
#for each row, get the probability of the column
try25=as.data.frame(apply(fit_hf, 2, prob25) )
colnames(try25)="estimate25"
try50=as.data.frame(apply(fit_hf, 2, prob50) )
colnames(try50)="estimate50"
try75=as.data.frame(apply(fit_hf, 2, prob75) )
colnames(try75)="estimate75"

MyData_hf=MyData_hf %>% cbind(tidyMCMC(as.mcmc(fit_hf),
                                       conf.int=T, conf.method='HPDinterval'))
MyData_hf=cbind(MyData_hf,try25,try50,try75)

Reserves_conditioned_hf<- MyData_hf[which(MyData_hf$CleanHabitat=="Slope" & MyData_hf$Protection=="UnfishedHigh" & MyData_hf$CensusMethod=="Standard belt transect" &
                                            MyData_hf$DepthCategory=="4-10m"),]
#trait diversity
prob25 <- function(x){length(x[exp(x)>rfunc*0.25])/length(x)} 
prob50 <- function(x){length(x[exp(x)>rfunc*0.50])/length(x)} 
prob75 <- function(x){length(x[exp(x)>rfunc*0.75])/length(x)} 

MyData_td<-expand.grid(sgrav_tot2=seq(range(alldata$sgrav_tot2)[[1]],range(alldata$sgrav_tot2)[[2]],length=101),
                       Protection=c("Fished","Restricted","UnfishedHigh"),
                       DepthCategory=c(">10m"," 0-4m","4-10m"),
                       CleanHabitat=c("Crest","Flat","Lagoon_Back reef","Slope"),
                       sTotal_sampling_area=0,
                       CensusMethod=c("Distance sampling","Point intercept","Standard belt transect"),
                       sRegional_population_growth=0,
                       sOcean_prod=0,
                       sHDI=0,
                       sReef_fish_landings_per_km2=0,
                       sClimate_stress=0,
                       sLarger_pop_size=0)

X <- model.matrix(~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot2+
                    sRegional_population_growth+sOcean_prod+sClimate_stress+
                    sHDI+sLarger_pop_size+sReef_fish_landings_per_km2,data=MyData_td)

coefs_td=as.matrix(fundiv_mpa_restricted)[, paste("b",row.names(fixef(fundiv_mpa_restricted)),sep="_")] #alternative

fit_td=coefs_td %*% t(X)
dim(fit_td)

#for each row, get the probability of the column
try25=as.data.frame(apply(fit_td, 2, prob25) )
colnames(try25)="estimate25"
try50=as.data.frame(apply(fit_td, 2, prob50) )
colnames(try50)="estimate50"
try75=as.data.frame(apply(fit_td, 2, prob75) )
colnames(try75)="estimate75"

MyData_td=MyData_td %>% cbind(tidyMCMC(as.mcmc(fit_td),
                                       conf.int=T, conf.method='HPDinterval'))
MyData_td=cbind(MyData_td,try25,try50,try75)

Reserves_conditioned_td<- MyData_td[which(MyData_td$CleanHabitat=="Slope" & MyData_td$Protection=="UnfishedHigh" & MyData_td$CensusMethod=="Standard belt transect" &
                                            MyData_td$DepthCategory=="4-10m"),]


#all three goals (assuming independence (multiply outcomes))
Reserves_conditioned_3 <-  as.data.frame(cbind((Reserves_conditioned_b20$estimate25*Reserves_conditioned_hf$estimate25*Reserves_conditioned_td$estimate25),(Reserves_conditioned_b20$estimate50*Reserves_conditioned_hf$estimate50*Reserves_conditioned_td$estimate50),(Reserves_conditioned_b20$estimate75*Reserves_conditioned_hf$estimate75*Reserves_conditioned_td$estimate75),Reserves_conditioned_b20$sgrav_tot2))
colnames(Reserves_conditioned_3)  =c("estimate25","estimate50","estimate75","sgrav_tot2")

#differences in probability by using the more stringent mpa categorization
mpa_cond_td <- ggplot(NULL)+
  geom_line(aes(x=Reserves_conditioned_td$sgrav_tot2,y=(Reserves_conditioned_td$estimate25-dat_Reserve_td$estimate25)),colour="#fbb4b9",size=1)+
  geom_line(aes(x=Reserves_conditioned_td$sgrav_tot2,y=(Reserves_conditioned_td$estimate50-dat_Reserve_td$estimate50)),colour="#f768a1",size=1)+
  geom_line(aes(x=Reserves_conditioned_td$sgrav_tot2,y=(Reserves_conditioned_td$estimate75-dat_Reserve_td$estimate75)),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous("",limits=c(-0.3,0.3),labels=NULL) + theme_classic()+xlab("")

mpa_cond_hf <- ggplot(NULL)+
  geom_line(aes(x=Reserves_conditioned_hf$sgrav_tot2,y=(Reserves_conditioned_hf$estimate25-dat_Reserve_hf$estimate25)),colour="#fbb4b9",size=1)+
  geom_line(aes(x=Reserves_conditioned_hf$sgrav_tot2,y=(Reserves_conditioned_hf$estimate50-dat_Reserve_hf$estimate50)),colour="#f768a1",size=1)+
  geom_line(aes(x=Reserves_conditioned_hf$sgrav_tot2,y=(Reserves_conditioned_hf$estimate75-dat_Reserve_hf$estimate75)),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous("",limits=c(-0.3,0.3),labels=NULL) + theme_classic()+xlab("")

mpa_cond_b20 <- ggplot(NULL)+
  geom_line(aes(x=Reserves_conditioned_b20$sgrav_tot2,y=(Reserves_conditioned_b20$estimate25-dat_Reserve_b20$estimate25)),colour="#fbb4b9",size=1)+
  geom_line(aes(x=Reserves_conditioned_b20$sgrav_tot2,y=(Reserves_conditioned_b20$estimate50-dat_Reserve_b20$estimate50)),colour="#f768a1",size=1)+
  geom_line(aes(x=Reserves_conditioned_b20$sgrav_tot2,y=(Reserves_conditioned_b20$estimate75-dat_Reserve_b20$estimate75)),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous("",limits=c(-0.3,0.3),labels=NULL) + theme_classic()+xlab("")


mpa_cond_3<- ggplot(NULL)+
  geom_line(aes(x=Reserves_conditioned_3$sgrav_tot2,y=(Reserves_conditioned_3$estimate25-dat_Reserve_3$estimate25)),colour="#fbb4b9",size=1)+
  geom_line(aes(x=Reserves_conditioned_3$sgrav_tot2,y=(Reserves_conditioned_3$estimate50-dat_Reserve_3$estimate50)),colour="#f768a1",size=1)+
  geom_line(aes(x=Reserves_conditioned_3$sgrav_tot2,y=(Reserves_conditioned_3$estimate75-dat_Reserve_3$estimate75)),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous("",limits=c(-0.3,0.3)) + theme_classic() +xlab("")

ggarrange(mpa_cond_3 ,mpa_cond_b20,mpa_cond_hf,mpa_cond_td ,nrow=1,ncol=4, widths=c(1.1,1,1,1))

#do the alluvial plots for the more stringent reserves
test_mpa=mpa_restricted
test_mpa[1:dim(test_mpa)[1],"CensusMethod"]="Standard belt transect"
test_mpa[1:dim(test_mpa)[1],"CleanHabitat"]="Slope"
test_mpa[1:dim(test_mpa)[1],"DepthCategory"]="4-10m"
test_mpa[1:dim(test_mpa)[1],"sTotal_sampling_area"]=0
test_mpa=droplevels(test_mpa)
B20test_mpa_posterior=posterior_predict(B20_mpa_restricted,newdata=test_mpa,re_formula = NULL)
mpa_restricted$B20_newdata=exp(colMeans(B20test_mpa_posterior))-1
reserves_mpa=test_mpa
reserves_mpa[1:dim(reserves_mpa)[1],"Protection"]="UnfishedHigh"
reserves_mpa=droplevels(reserves_mpa)
B20reserves_mpa_posterior=posterior_predict(B20_mpa_restricted,newdata=reserves_mpa,re_formula = NULL)
mpa_restricted$B20_newdata_reserves=exp(colMeans(B20reserves_mpa_posterior))-1

prob25 <- function(x){length(x[(exp(x)-1)>rpred*0.25])/length(x)} 
prob50 <- function(x){length(x[(exp(x)-1)>rpred*0.50])/length(x)} 
prob75 <- function(x){length(x[(exp(x)-1)>rpred*0.75])/length(x)} 

#for each row, get the probability of the column
newdatampa_25=as.data.frame(apply(B20test_mpa_posterior, 2, prob25) )
colnames(newdatampa_25)="B20newdata25"
newdatampa_50=as.data.frame(apply(B20test_mpa_posterior, 2, prob50) )
colnames(newdatampa_50)="B20newdata50"
newdatampa_75=as.data.frame(apply(B20test_mpa_posterior, 2, prob75) )
colnames(newdatampa_75)="B20newdata75"
newdata_mpa=cbind(mpa_restricted,newdatampa_25,newdatampa_50,newdatampa_75)

reserves25_mpa=as.data.frame(apply(B20reserves_mpa_posterior, 2, prob25) )
colnames(reserves25_mpa)="B20reserves25"
reserves50_mpa=as.data.frame(apply(B20reserves_mpa_posterior, 2, prob50) )
colnames(reserves50_mpa)="B20reserves50"
reserves75_mpa=as.data.frame(apply(B20reserves_mpa_posterior, 2, prob75) )
colnames(reserves75_mpa)="B20reserves75"
newdata_mpa=cbind(newdata_mpa,reserves25_mpa,reserves50_mpa,reserves75_mpa)


#trait diversity
fundivtest_posterior_mpa=posterior_predict(fundiv_mpa_restricted,newdata=test_mpa[!is.na(test_mpa$Trait_diversity),],re_formula = NULL)
newdata_mpa[!is.na(newdata_mpa$Trait_diversity),"fundiv_newdata"]=exp(colMeans(fundivtest_posterior_mpa))
fundivreserves_posterior_mpa=posterior_predict(fundiv_mpa_restricted,newdata=reserves_mpa[!is.na(reserves_mpa$Trait_diversity),],re_formula = NULL)
newdata_mpa[!is.na(newdata_mpa$Trait_diversity),"fundiv_newdata_reserves"]=exp(colMeans(fundivreserves_posterior_mpa))

prob25 <- function(x){length(x[exp(x)>rfunc*0.25])/length(x)} 
prob50 <- function(x){length(x[exp(x)>rfunc*0.50])/length(x)} 
prob75 <- function(x){length(x[exp(x)>rfunc*0.75])/length(x)} 

newdata_mpa[!is.na(newdata_mpa$Trait_diversity),"fundivnewdata25"]=as.data.frame(apply(fundivtest_posterior_mpa, 2, prob25) )
newdata_mpa[!is.na(newdata_mpa$Trait_diversity),"fundivnewdata50"]=as.data.frame(apply(fundivtest_posterior_mpa, 2, prob50) )
newdata_mpa[!is.na(newdata_mpa$Trait_diversity),"fundivnewdata75"]=as.data.frame(apply(fundivtest_posterior_mpa, 2, prob75) )


newdata_mpa[!is.na(newdata_mpa$Trait_diversity),"fundivreserves25"]=as.data.frame(apply(fundivreserves_posterior_mpa, 2, prob25) )
newdata_mpa[!is.na(newdata_mpa$Trait_diversity),"fundivreserves50"]=as.data.frame(apply(fundivreserves_posterior_mpa, 2, prob50) )
newdata_mpa[!is.na(newdata_mpa$Trait_diversity),"fundivreserves75"]=as.data.frame(apply(fundivreserves_posterior_mpa, 2, prob75) )


#parrotfish sraping
hertest_posterior_mpa=posterior_predict(her_mpa_restricted,newdata=test_mpa[!is.na(test_mpa$Scraping_potential),],re_formula = NULL)
newdata_mpa[!is.na(newdata_mpa$Scraping_potential),"her_newdata"]=colMeans(hertest_posterior_mpa)
herreserves_posterior_mpa=posterior_predict(her_mpa_restricted,newdata=reserves_mpa[!is.na(reserves_mpa$Scraping_potential),],re_formula = NULL)
newdata_mpa[!is.na(newdata_mpa$Scraping_potential),"her_newdata_reserves"]=colMeans(herreserves_posterior_mpa)

prob25 <- function(x){length(x[x>rherb*0.25])/length(x)} 
prob50 <- function(x){length(x[x>rherb*0.50])/length(x)} 
prob75 <- function(x){length(x[x>rherb*0.75])/length(x)} 

newdata_mpa[!is.na(newdata_mpa$Scraping_potential),"hernewdata25"]=as.data.frame(apply(hertest_posterior_mpa, 2, prob25) )
newdata_mpa[!is.na(newdata_mpa$Scraping_potential),"hernewdata50"]=as.data.frame(apply(hertest_posterior_mpa, 2, prob50) )
newdata_mpa[!is.na(newdata_mpa$Scraping_potential),"hernewdata75"]=as.data.frame(apply(hertest_posterior_mpa, 2, prob75) )


newdata_mpa[!is.na(newdata_mpa$Scraping_potential),"herreserves25"]=as.data.frame(apply(herreserves_posterior_mpa, 2, prob25) )
newdata_mpa[!is.na(newdata_mpa$Scraping_potential),"herreserves50"]=as.data.frame(apply(herreserves_posterior_mpa, 2, prob50) )
newdata_mpa[!is.na(newdata_mpa$Scraping_potential),"herreserves75"]=as.data.frame(apply(herreserves_posterior_mpa, 2, prob75) )


#all three together
newdata_mpa$totnewdata25=newdata_mpa$B20newdata25*newdata_mpa$hernewdata25*newdata_mpa$fundivnewdata25
newdata_mpa$totnewdata50=newdata_mpa$B20newdata50*newdata_mpa$hernewdata50*newdata_mpa$fundivnewdata50
newdata_mpa$totnewdata75=newdata_mpa$B20newdata75*newdata_mpa$hernewdata75*newdata_mpa$fundivnewdata75
newdata_mpa$totreserves25=newdata_mpa$B20reserves25*newdata_mpa$herreserves25*newdata_mpa$fundivreserves25
newdata_mpa$totreserves50=newdata_mpa$B20reserves50*newdata_mpa$herreserves50*newdata_mpa$fundivreserves50
newdata_mpa$totreserves75=newdata_mpa$B20reserves75*newdata_mpa$herreserves75*newdata_mpa$fundivreserves75

#fished sites
fished_newdata_mpa=newdata_mpa[newdata_mpa$Protection=="Fished",]

#B20
fished_newdata_mpa$predicted_pred25=ifelse(fished_newdata_mpa$B20_newdata<rpred*0.25,0,1)
fished_newdata_mpa$predicted_pred50=ifelse(fished_newdata_mpa$B20_newdata<rpred*0.50,0,1)
fished_newdata_mpa$predicted_pred75=ifelse(fished_newdata_mpa$B20_newdata<rpred*0.75,0,1)
fished_newdata_mpa$unfished_pred25=ifelse(fished_newdata_mpa$B20_newdata_reserves<rpred*0.25,0,1)
fished_newdata_mpa$unfished_pred50=ifelse(fished_newdata_mpa$B20_newdata_reserves<rpred*0.50,0,1)
fished_newdata_mpa$unfished_pred75=ifelse(fished_newdata_mpa$B20_newdata_reserves<rpred*0.75,0,1)

fished_newdata_mpa$unfished_changeB20=ifelse(fished_newdata_mpa$predicted_pred25==0 &fished_newdata_mpa$unfished_pred75==1, "0to75", 
                                             ifelse(fished_newdata_mpa$predicted_pred25==0 & fished_newdata_mpa$unfished_pred50==1 & fished_newdata_mpa$unfished_pred75==0, "0to50",
                                                    ifelse(fished_newdata_mpa$predicted_pred25==0 & fished_newdata_mpa$unfished_pred25==1 & fished_newdata_mpa$unfished_pred50==0, "0to25",
                                                           ifelse(fished_newdata_mpa$predicted_pred25==0 & fished_newdata_mpa$unfished_pred25==0, "0to0",
                                                                  ifelse(fished_newdata_mpa$predicted_pred25==1 &fished_newdata_mpa$predicted_pred50==0 &fished_newdata_mpa$unfished_pred75==1, "25to75",
                                                                         ifelse(fished_newdata_mpa$predicted_pred25==1 &fished_newdata_mpa$predicted_pred50==1&fished_newdata_mpa$predicted_pred75==0 &fished_newdata_mpa$unfished_pred75==1, "50to75",
                                                                                ifelse(fished_newdata_mpa$predicted_pred25==1 &fished_newdata_mpa$predicted_pred50==1&fished_newdata_mpa$predicted_pred75==1 &fished_newdata_mpa$unfished_pred75==1, "75to75",
                                                                                       ifelse(fished_newdata_mpa$predicted_pred25==1 &fished_newdata_mpa$predicted_pred50==0 &fished_newdata_mpa$unfished_pred75==0&fished_newdata_mpa$unfished_pred50==1, "25to50",
                                                                                              ifelse(fished_newdata_mpa$predicted_pred25==1 &fished_newdata_mpa$predicted_pred50==1 &fished_newdata_mpa$unfished_pred75==0&fished_newdata_mpa$unfished_pred50==1, "50to50",
                                                                                                     ifelse(fished_newdata_mpa$predicted_pred25==1 &fished_newdata_mpa$predicted_pred50==0 &fished_newdata_mpa$unfished_pred25==1 &fished_newdata_mpa$unfished_pred75==0&fished_newdata_mpa$unfished_pred50==0, "25to25",-999))))))))))

matrixpred=matrix(NA,nrow=14,ncol=5)
matrixpred[1:4,1]=rep(0,4)
matrixpred[1:4,2]=c(1:4)
matrixpred[1:4,3]=c(length(fished_newdata_mpa$B20_newdata[fished_newdata_mpa$B20_newdata>rpred*0.75]),length(fished_newdata_mpa$B20_newdata[fished_newdata_mpa$B20_newdata>rpred*0.50&fished_newdata_mpa$B20_newdata<rpred*0.75]),length(fished_newdata_mpa$B20_newdata[fished_newdata_mpa$B20_newdata>rpred*0.25&fished_newdata_mpa$B20_newdata<rpred*0.50]),length(fished_newdata_mpa$B20_newdata[fished_newdata_mpa$B20_newdata<rpred*0.25]))
matrixpred[1:4,4]=rep("Initial", 4)
matrixpred[1:4,5]=c("75-100%","50-75%","25-50%","0-25%")
matrixpred[5:8,1]=rep(4,4)
matrixpred[5:8,2]=c(8,7,6,5)
matrixpred[5:8,4]=rep("reserve",4)
matrixpred[5:8,5]=c("0to0","0to25","0to50","0to75")
matrixpred[5:8,3]=c(length(fished_newdata_mpa$unfished_changeB20[fished_newdata_mpa$unfished_changeB20=="0to0"]),length(fished_newdata_mpa$unfished_changeB20[fished_newdata_mpa$unfished_changeB20=="0to25"]),length(fished_newdata_mpa$unfished_changeB20[fished_newdata_mpa$unfished_changeB20=="0to50"]),length(fished_newdata_mpa$unfished_changeB20[fished_newdata_mpa$unfished_changeB20=="0to75"]))
matrixpred[9:11,1]=rep(3,3)
matrixpred[9:11,2]=c(7,6,5)
matrixpred[9:11,4]=rep("reserve",3)
matrixpred[9:11,5]=c("25to25","25to50","25to75")
matrixpred[9:11,3]=c(length(fished_newdata_mpa$unfished_changeB20[fished_newdata_mpa$unfished_changeB20=="25to25"]),length(fished_newdata_mpa$unfished_changeB20[fished_newdata_mpa$unfished_changeB20=="25to50"]),length(fished_newdata_mpa$unfished_changeB20[fished_newdata_mpa$unfished_changeB20=="25to75"]))
matrixpred[12:13,1]=rep(2,2)
matrixpred[12:13,2]=c(6,5)
matrixpred[12:13,4]=rep("reserve",2)
matrixpred[12:13,5]=c("50to50","50to75")
matrixpred[12:13,3]=c(length(fished_newdata_mpa$unfished_changeB20[fished_newdata_mpa$unfished_changeB20=="50to50"]),length(fished_newdata_mpa$unfished_changeB20[fished_newdata_mpa$unfished_changeB20=="50to75"]))
matrixpred[14,1]=1
matrixpred[14,2]=5
matrixpred[14,4]="reserve"
matrixpred[14,5]="75to75"
matrixpred[14,3]=length(fished_newdata_mpa$unfished_changeB20[fished_newdata_mpa$unfished_changeB20=="75to75"])

#alluvial charts
pred=as.data.frame(matrixpred)
colnames(pred)=c("from","to","value","scenario","clarified")
pred$from=as.numeric(as.character(pred$from))
pred$to=as.numeric(as.character(pred$to))
pred$value=as.numeric(as.character(pred$value))
pred=pred[order(pred$from),]
nodes <- data.frame(c("75-100%","50-75%","25-50%","0-25%","75-100%","50-75%","25-50%","0-25%"))
names(nodes) <- "name"
nodes$group = as.factor(c("group1", "group2","group3", "group4", "group5", "group6", "group7", "group8"))

MR <- pred
link <- MR[,c(1:3)]
link <- link[-c(1:4),]
link$from <- link$from-1
link$to <- link$to-1
link$zeroout=ifelse(link$value==0,1,0)
link=link[link$zeroout==0,]
link$zeroout=NULL
link$group=ifelse(link$from==0,"group1", ifelse(link$from==1,"group2", ifelse(link$from==2,"group3", "group4")))
nodes$sites=c(round(MR[1:4,3],digits=3),round(sum(fished_newdata_mpa$unfished_pred75),digits=3),round(sum(fished_newdata_mpa$unfished_pred50[fished_newdata_mpa$unfished_pred75==0]),digits=3),round(sum(fished_newdata_mpa$unfished_pred25[fished_newdata_mpa$unfished_pred50==0]),digits=3),round(length(fished_newdata_mpa$unfished_pred25[fished_newdata_mpa$unfished_pred25==0]),digits=3))
nodes$sites=c(round(MR[1:4,3]/length(fished_newdata_mpa$unfished_pred75),digits=3)*100,round(sum(fished_newdata_mpa$unfished_pred75)/length(fished_newdata_mpa$unfished_pred75),digits=3)*100,round(sum(fished_newdata_mpa$unfished_pred50[fished_newdata_mpa$unfished_pred75==0])/length(fished_newdata_mpa$unfished_pred75),digits=3)*100,round(sum(fished_newdata_mpa$unfished_pred25[fished_newdata_mpa$unfished_pred50==0])/length(fished_newdata_mpa$unfished_pred75),digits=3)*100,round(length(fished_newdata_mpa$unfished_pred25[fished_newdata_mpa$unfished_pred25==0])/length(fished_newdata_mpa$unfished_pred75),digits=3)*100)
nodes$sitesnolabel=rep("",length(nodes$sites))
networkMR_pred_mpa <- sankeyNetwork(Links = link, Nodes = nodes, Source = "from",
                                    Target = "to", Value = "value", NodeID = "sites",
                                    units = , fontSize = 14, nodeWidth = 50,colourScale=my_color, NodeGroup = "group", LinkGroup = "group", width=300)

networkMR_pred_mpa


#trait diversity
fished_newdata_mpa2=fished_newdata_mpa[!is.na(fished_newdata_mpa$fundiv_newdata),]
fished_newdata_mpa2$predicted_func25=ifelse(fished_newdata_mpa2$fundiv_newdata<rfunc*0.25,0,1)
fished_newdata_mpa2$predicted_func50=ifelse(fished_newdata_mpa2$fundiv_newdata<rfunc*0.50,0,1)
fished_newdata_mpa2$predicted_func75=ifelse(fished_newdata_mpa2$fundiv_newdata<rfunc*0.75,0,1)
fished_newdata_mpa2$unfished_func25=ifelse(fished_newdata_mpa2$fundiv_newdata_reserves<rfunc*0.25,0,1)
fished_newdata_mpa2$unfished_func50=ifelse(fished_newdata_mpa2$fundiv_newdata_reserves<rfunc*0.50,0,1)
fished_newdata_mpa2$unfished_func75=ifelse(fished_newdata_mpa2$fundiv_newdata_reserves<rfunc*0.75,0,1)

fished_newdata_mpa2$unfished_changefundiv=ifelse(fished_newdata_mpa2$predicted_func25==0 &fished_newdata_mpa2$unfished_func75==1, "0to75", 
                                                 ifelse(fished_newdata_mpa2$predicted_func25==0 & fished_newdata_mpa2$unfished_func50==1 & fished_newdata_mpa2$unfished_func75==0, "0to50",
                                                        ifelse(fished_newdata_mpa2$predicted_func25==0 & fished_newdata_mpa2$unfished_func25==1 & fished_newdata_mpa2$unfished_func50==0, "0to25",
                                                               ifelse(fished_newdata_mpa2$predicted_func25==0 & fished_newdata_mpa2$unfished_func25==0, "0to0",
                                                                      ifelse(fished_newdata_mpa2$predicted_func25==1 &fished_newdata_mpa2$predicted_func50==0 &fished_newdata_mpa2$unfished_func75==1, "25to75",
                                                                             ifelse(fished_newdata_mpa2$predicted_func25==1 &fished_newdata_mpa2$predicted_func50==1&fished_newdata_mpa2$predicted_func75==0 &fished_newdata_mpa2$unfished_func75==1, "50to75",
                                                                                    ifelse(fished_newdata_mpa2$predicted_func25==1 &fished_newdata_mpa2$predicted_func50==1&fished_newdata_mpa2$predicted_func75==1 &fished_newdata_mpa2$unfished_func75==1, "75to75",
                                                                                           ifelse(fished_newdata_mpa2$predicted_func25==1 &fished_newdata_mpa2$predicted_func50==0 &fished_newdata_mpa2$unfished_func75==0&fished_newdata_mpa2$unfished_func50==1, "25to50",
                                                                                                  ifelse(fished_newdata_mpa2$predicted_func25==1 &fished_newdata_mpa2$predicted_func50==1 &fished_newdata_mpa2$unfished_func75==0&fished_newdata_mpa2$unfished_func50==1, "50to50",
                                                                                                         ifelse(fished_newdata_mpa2$predicted_func25==1 &fished_newdata_mpa2$predicted_func50==0 &fished_newdata_mpa2$unfished_func25==1 &fished_newdata_mpa2$unfished_func75==0&fished_newdata_mpa2$unfished_func50==0, "25to25",-999))))))))))

matrixfunc=matrix(NA,nrow=14,ncol=5)
matrixfunc[1:4,1]=rep(0,4)
matrixfunc[1:4,2]=c(1:4)
matrixfunc[1:4,3]=c(length(fished_newdata_mpa2$fundiv_newdata[fished_newdata_mpa2$fundiv_newdata>rfunc*0.75]),length(fished_newdata_mpa2$fundiv_newdata[fished_newdata_mpa2$fundiv_newdata>rfunc*0.50&fished_newdata_mpa2$fundiv_newdata<rfunc*0.75]),length(fished_newdata_mpa2$fundiv_newdata[fished_newdata_mpa2$fundiv_newdata>rfunc*0.25&fished_newdata_mpa2$fundiv_newdata<rfunc*0.50]),length(fished_newdata_mpa2$fundiv_newdata[fished_newdata_mpa2$fundiv_newdata<rfunc*0.25]))
matrixfunc[1:4,4]=rep("Initial", 4)
matrixfunc[1:4,5]=c("75-100%","50-75%","25-50%","0-25%")
matrixfunc[5:8,1]=rep(4,4)
matrixfunc[5:8,2]=c(8,7,6,5)
matrixfunc[5:8,4]=rep("reserve",4)
matrixfunc[5:8,5]=c("0to0","0to25","0to50","0to75")
matrixfunc[5:8,3]=c(length(fished_newdata_mpa2$unfished_changefundiv[fished_newdata_mpa2$unfished_changefundiv=="0to0"]),length(fished_newdata_mpa2$unfished_changefundiv[fished_newdata_mpa2$unfished_changefundiv=="0to25"]),length(fished_newdata_mpa2$unfished_changefundiv[fished_newdata_mpa2$unfished_changefundiv=="0to50"]),length(fished_newdata_mpa2$unfished_changefundiv[fished_newdata_mpa2$unfished_changefundiv=="0to75"]))
matrixfunc[9:11,1]=rep(3,3)
matrixfunc[9:11,2]=c(7,6,5)
matrixfunc[9:11,4]=rep("reserve",3)
matrixfunc[9:11,5]=c("25to25","25to50","25to75")
matrixfunc[9:11,3]=c(length(fished_newdata_mpa2$unfished_changefundiv[fished_newdata_mpa2$unfished_changefundiv=="25to25"]),length(fished_newdata_mpa2$unfished_changefundiv[fished_newdata_mpa2$unfished_changefundiv=="25to50"]),length(fished_newdata_mpa2$unfished_changefundiv[fished_newdata_mpa2$unfished_changefundiv=="25to75"]))
matrixfunc[12:13,1]=rep(2,2)
matrixfunc[12:13,2]=c(6,5)
matrixfunc[12:13,4]=rep("reserve",2)
matrixfunc[12:13,5]=c("50to50","50to75")
matrixfunc[12:13,3]=c(length(fished_newdata_mpa2$unfished_changefundiv[fished_newdata_mpa2$unfished_changefundiv=="50to50"]),length(fished_newdata_mpa2$unfished_changefundiv[fished_newdata_mpa2$unfished_changefundiv=="50to75"]))
matrixfunc[14,1]=1
matrixfunc[14,2]=5
matrixfunc[14,4]="reserve"
matrixfunc[14,5]="75to75"
matrixfunc[14,3]=length(fished_newdata_mpa2$unfished_changefundiv[fished_newdata_mpa2$unfished_changefundiv=="75to75"])

#alluvial charts
fundiv=as.data.frame(matrixfunc)
colnames(fundiv)=c("from","to","value","scenario","clarified")
fundiv$from=as.numeric(as.character(fundiv$from))
fundiv$to=as.numeric(as.character(fundiv$to))
fundiv$value=as.numeric(as.character(fundiv$value))
fundiv=fundiv[order(fundiv$from),]
nodes <- data.frame(c("75-100%","50-75%","25-50%","0-25%","75-100%","50-75%","25-50%","0-25%"))
names(nodes) <- "name"
nodes$group = as.factor(c("group1", "group2","group3", "group4", "group5", "group6", "group7", "group8"))
MR <- fundiv
link <- MR[,c(1:3)]
link <- link[-c(1:4),]
link$from <- link$from-1
link$to <- link$to-1
link$zeroout=ifelse(link$value==0,1,0)
link=link[link$zeroout==0,]
link$zeroout=NULL
link$group=ifelse(link$from==0,"group1", ifelse(link$from==1,"group2", ifelse(link$from==2,"group3", "group4")))
nodes$sites=c(round(MR[1:4,3],digits=3),round(sum(fished_newdata_mpa2$unfished_func75),digits=3),round(sum(fished_newdata_mpa2$unfished_func50[fished_newdata_mpa2$unfished_func75==0]),digits=3),round(sum(fished_newdata_mpa2$unfished_func25[fished_newdata_mpa2$unfished_func50==0]),digits=3),round(length(fished_newdata_mpa2$unfished_func25[fished_newdata_mpa2$unfished_func25==0]),digits=3))
nodes$sites=c(round(MR[1:4,3]/length(fished_newdata_mpa2$unfished_func75),digits=3)*100,round(sum(fished_newdata_mpa2$unfished_func75)/length(fished_newdata_mpa2$unfished_func75),digits=3)*100,round(sum(fished_newdata_mpa2$unfished_func50[fished_newdata_mpa2$unfished_func75==0])/length(fished_newdata_mpa2$unfished_func75),digits=3)*100,round(sum(fished_newdata_mpa2$unfished_func25[fished_newdata_mpa2$unfished_func50==0])/length(fished_newdata_mpa2$unfished_func75),digits=3)*100,round(length(fished_newdata_mpa2$unfished_func25[fished_newdata_mpa2$unfished_func25==0])/length(fished_newdata_mpa2$unfished_func75),digits=3)*100)
networkMR_func_mpa <- sankeyNetwork(Links = link, Nodes = nodes, Source = "from",
                                    Target = "to", Value = "value", NodeID = "sites",
                                    units = "", fontSize = 13, nodeWidth = 50,colourScale=my_color, NodeGroup = "group", LinkGroup = "group", width=300)


networkMR_func_mpa 


#parrotfish scraping
fished_newdata_mpa3=fished_newdata_mpa2[!is.na(fished_newdata_mpa2$her_newdata),]
fished_newdata_mpa3$predicted_herb25=ifelse(fished_newdata_mpa3$her_newdata<rherb*0.25,0,1)
fished_newdata_mpa3$predicted_herb50=ifelse(fished_newdata_mpa3$her_newdata<rherb*0.50,0,1)
fished_newdata_mpa3$predicted_herb75=ifelse(fished_newdata_mpa3$her_newdata<rherb*0.75,0,1)
fished_newdata_mpa3$unfished_herb25=ifelse(fished_newdata_mpa3$her_newdata_reserves<rherb*0.25,0,1)
fished_newdata_mpa3$unfished_herb50=ifelse(fished_newdata_mpa3$her_newdata_reserves<rherb*0.50,0,1)
fished_newdata_mpa3$unfished_herb75=ifelse(fished_newdata_mpa3$her_newdata_reserves<rherb*0.75,0,1)

fished_newdata_mpa3$unfished_changeherb=ifelse(fished_newdata_mpa3$predicted_herb25==0 &fished_newdata_mpa3$unfished_herb75==1, "0to75", 
                                               ifelse(fished_newdata_mpa3$predicted_herb25==0 & fished_newdata_mpa3$unfished_herb50==1 & fished_newdata_mpa3$unfished_herb75==0, "0to50",
                                                      ifelse(fished_newdata_mpa3$predicted_herb25==0 & fished_newdata_mpa3$unfished_herb25==1 & fished_newdata_mpa3$unfished_herb50==0, "0to25",
                                                             ifelse(fished_newdata_mpa3$predicted_herb25==0 & fished_newdata_mpa3$unfished_herb25==0, "0to0",
                                                                    ifelse(fished_newdata_mpa3$predicted_herb25==1 &fished_newdata_mpa3$predicted_herb50==0 &fished_newdata_mpa3$unfished_herb75==1, "25to75",
                                                                           ifelse(fished_newdata_mpa3$predicted_herb25==1 &fished_newdata_mpa3$predicted_herb50==1&fished_newdata_mpa3$predicted_herb75==0 &fished_newdata_mpa3$unfished_herb75==1, "50to75",
                                                                                  ifelse(fished_newdata_mpa3$predicted_herb25==1 &fished_newdata_mpa3$predicted_herb50==1&fished_newdata_mpa3$predicted_herb75==1 &fished_newdata_mpa3$unfished_herb75==1, "75to75",
                                                                                         ifelse(fished_newdata_mpa3$predicted_herb25==1 &fished_newdata_mpa3$predicted_herb50==0 &fished_newdata_mpa3$unfished_herb75==0&fished_newdata_mpa3$unfished_herb50==1, "25to50",
                                                                                                ifelse(fished_newdata_mpa3$predicted_herb25==1 &fished_newdata_mpa3$predicted_herb50==1 &fished_newdata_mpa3$unfished_herb75==0&fished_newdata_mpa3$unfished_herb50==1, "50to50",
                                                                                                       ifelse(fished_newdata_mpa3$predicted_herb25==1 &fished_newdata_mpa3$predicted_herb50==0 &fished_newdata_mpa3$unfished_herb25==1 &fished_newdata_mpa3$unfished_herb75==0&fished_newdata_mpa3$unfished_herb50==0, "25to25",-999))))))))))
matrixherb=matrix(NA,nrow=14,ncol=5)
matrixherb[1:4,1]=rep(0,4)
matrixherb[1:4,2]=c(1:4)
matrixherb[1:4,3]=c(length(fished_newdata_mpa3$her_newdata[fished_newdata_mpa3$her_newdata>rherb*0.75]),length(fished_newdata_mpa3$her_newdata[fished_newdata_mpa3$her_newdata>rherb*0.50&fished_newdata_mpa3$her_newdata<rherb*0.75]),length(fished_newdata_mpa3$her_newdata[fished_newdata_mpa3$her_newdata>rherb*0.25&fished_newdata_mpa3$her_newdata<rherb*0.50]),length(fished_newdata_mpa3$her_newdata[fished_newdata_mpa3$her_newdata<rherb*0.25]))
matrixherb[1:4,4]=rep("Initial", 4)
matrixherb[1:4,5]=c("75-100%","50-75%","25-50%","0-25%")
matrixherb[5:8,1]=rep(4,4)
matrixherb[5:8,2]=c(8,7,6,5)
matrixherb[5:8,4]=rep("reserve",4)
matrixherb[5:8,5]=c("0to0","0to25","0to50","0to75")
matrixherb[5:8,3]=c(length(fished_newdata_mpa3$unfished_changeherb[fished_newdata_mpa3$unfished_changeherb=="0to0"]),length(fished_newdata_mpa3$unfished_changeherb[fished_newdata_mpa3$unfished_changeherb=="0to25"]),length(fished_newdata_mpa3$unfished_changeherb[fished_newdata_mpa3$unfished_changeherb=="0to50"]),length(fished_newdata_mpa3$unfished_changeherb[fished_newdata_mpa3$unfished_changeherb=="0to75"]))
matrixherb[9:11,1]=rep(3,3)
matrixherb[9:11,2]=c(7,6,5)
matrixherb[9:11,4]=rep("reserve",3)
matrixherb[9:11,5]=c("25to25","25to50","25to75")
matrixherb[9:11,3]=c(length(fished_newdata_mpa3$unfished_changeherb[fished_newdata_mpa3$unfished_changeherb=="25to25"]),length(fished_newdata_mpa3$unfished_changeherb[fished_newdata_mpa3$unfished_changeherb=="25to50"]),length(fished_newdata_mpa3$unfished_changeherb[fished_newdata_mpa3$unfished_changeherb=="25to75"]))
matrixherb[12:13,1]=rep(2,2)
matrixherb[12:13,2]=c(6,5)
matrixherb[12:13,4]=rep("reserve",2)
matrixherb[12:13,5]=c("50to50","50to75")
matrixherb[12:13,3]=c(length(fished_newdata_mpa3$unfished_changeherb[fished_newdata_mpa3$unfished_changeherb=="50to50"]),length(fished_newdata_mpa3$unfished_changeherb[fished_newdata_mpa3$unfished_changeherb=="50to75"]))
matrixherb[14,1]=1
matrixherb[14,2]=5
matrixherb[14,4]="reserve"
matrixherb[14,5]="75to75"
matrixherb[14,3]=length(fished_newdata_mpa3$unfished_changeherb[fished_newdata_mpa3$unfished_changeherb=="75to75"])

#alluvial charts
herb=as.data.frame(matrixherb)
colnames(herb)=c("from","to","value","scenario","clarified")
herb$from=as.numeric(as.character(herb$from))
herb$to=as.numeric(as.character(herb$to))
herb$value=as.numeric(as.character(herb$value))
herb=herb[order(herb$from),]
nodes <- data.frame(c("75-100%","50-75%","25-50%","0-25%","75-100%","50-75%","25-50%","0-25%"))
names(nodes) <- "name"
nodes$group = as.factor(c("group1", "group2","group3", "group4", "group5", "group6", "group7", "group8"))
MR <- herb
link <- MR[,c(1:3)]
link <- link[-c(1:4),]
link$from <- link$from-1
link$to <- link$to-1
link$zeroout=ifelse(link$value==0,1,0)
link=link[link$zeroout==0,]
link$zeroout=NULL
link$group=ifelse(link$from==0,"group1", ifelse(link$from==1,"group2", ifelse(link$from==2,"group3", "group4")))
nodes$sites=c(round(MR[1:4,3],digits=3),round(sum(fished_newdata_mpa3$unfished_herb75),digits=3),round(sum(fished_newdata_mpa3$unfished_herb50[fished_newdata_mpa3$unfished_herb75==0]),digits=3),round(sum(fished_newdata_mpa3$unfished_herb25[fished_newdata_mpa3$unfished_herb50==0]),digits=3),round(length(fished_newdata_mpa3$unfished_herb25[fished_newdata_mpa3$unfished_herb25==0]),digits=3))
nodes$sites=c(round(MR[1:4,3]/length(fished_newdata_mpa3$unfished_herb75),digits=3)*100,round(sum(fished_newdata_mpa3$unfished_herb75)/length(fished_newdata_mpa3$unfished_herb75),digits=3)*100,round(sum(fished_newdata_mpa3$unfished_herb50[fished_newdata_mpa3$unfished_herb75==0])/length(fished_newdata_mpa3$unfished_herb75),digits=3)*100,round(sum(fished_newdata_mpa3$unfished_herb25[fished_newdata_mpa3$unfished_herb50==0])/length(fished_newdata_mpa3$unfished_herb75),digits=3)*100,round(length(fished_newdata_mpa3$unfished_herb25[fished_newdata_mpa3$unfished_herb25==0])/length(fished_newdata_mpa3$unfished_herb75),digits=3)*100)
networkMR_herb_mpa <- sankeyNetwork(Links = link, Nodes = nodes, Source = "from",
                                    Target = "to", Value = "value", NodeID = "sites",
                                    units = "", fontSize = 13, nodeWidth = 50,colourScale=my_color, NodeGroup = "group", LinkGroup = "group", width=300)


networkMR_herb_mpa 

#all three goals
fished_newdata_mpa3$predicted_tot25=ifelse(fished_newdata_mpa3$predicted_pred25==1&fished_newdata_mpa3$predicted_func25==1&fished_newdata_mpa3$predicted_herb25==1,1,0)
fished_newdata_mpa3$predicted_tot50=ifelse(fished_newdata_mpa3$predicted_pred50==1&fished_newdata_mpa3$predicted_func50==1&fished_newdata_mpa3$predicted_herb50==1,1,0)
fished_newdata_mpa3$predicted_tot75=ifelse(fished_newdata_mpa3$predicted_pred75==1&fished_newdata_mpa3$predicted_func75==1&fished_newdata_mpa3$predicted_herb75==1,1,0)
fished_newdata_mpa3$unfished_tot25=ifelse(fished_newdata_mpa3$unfished_pred25==1&fished_newdata_mpa3$unfished_func25==1&fished_newdata_mpa3$unfished_herb25==1,1,0)
fished_newdata_mpa3$unfished_tot50=ifelse(fished_newdata_mpa3$unfished_pred50==1&fished_newdata_mpa3$unfished_func50==1&fished_newdata_mpa3$unfished_herb50==1,1,0)
fished_newdata_mpa3$unfished_tot75=ifelse(fished_newdata_mpa3$unfished_pred75==1&fished_newdata_mpa3$unfished_func75==1&fished_newdata_mpa3$unfished_herb75==1,1,0)

fished_newdata_mpa3$unfished_changetot=ifelse(fished_newdata_mpa3$predicted_tot25==0 &fished_newdata_mpa3$unfished_tot75==1, "0to75", 
                                              ifelse(fished_newdata_mpa3$predicted_tot25==0 & fished_newdata_mpa3$unfished_tot50==1 & fished_newdata_mpa3$unfished_tot75==0, "0to50",
                                                     ifelse(fished_newdata_mpa3$predicted_tot25==0 & fished_newdata_mpa3$unfished_tot25==1 & fished_newdata_mpa3$unfished_tot50==0, "0to25",
                                                            ifelse(fished_newdata_mpa3$predicted_tot25==0 & fished_newdata_mpa3$unfished_tot25==0, "0to0",
                                                                   ifelse(fished_newdata_mpa3$predicted_tot25==1 &fished_newdata_mpa3$predicted_tot50==0 &fished_newdata_mpa3$unfished_tot75==1, "25to75",
                                                                          ifelse(fished_newdata_mpa3$predicted_tot25==1 &fished_newdata_mpa3$predicted_tot50==1&fished_newdata_mpa3$predicted_tot75==0 &fished_newdata_mpa3$unfished_tot75==1, "50to75",
                                                                                 ifelse(fished_newdata_mpa3$predicted_tot25==1 &fished_newdata_mpa3$predicted_tot50==1&fished_newdata_mpa3$predicted_tot75==1 &fished_newdata_mpa3$unfished_tot75==1, "75to75",
                                                                                        ifelse(fished_newdata_mpa3$predicted_tot25==1 &fished_newdata_mpa3$predicted_tot50==0 &fished_newdata_mpa3$unfished_tot75==0&fished_newdata_mpa3$unfished_tot50==1, "25to50",
                                                                                               ifelse(fished_newdata_mpa3$predicted_tot25==1 &fished_newdata_mpa3$predicted_tot50==1 &fished_newdata_mpa3$unfished_tot75==0&fished_newdata_mpa3$unfished_tot50==1, "50to50",
                                                                                                      ifelse(fished_newdata_mpa3$predicted_tot25==1 &fished_newdata_mpa3$predicted_tot50==0 &fished_newdata_mpa3$unfished_tot25==1 &fished_newdata_mpa3$unfished_tot75==0&fished_newdata_mpa3$unfished_tot50==0, "25to25",-999))))))))))
matrixtot=matrix(NA,nrow=14,ncol=5)
matrixtot[1:4,1]=rep(0,4)
matrixtot[1:4,2]=c(1:4)
matrixtot[1:4,3]=c(length(fished_newdata_mpa3$predicted_tot25[fished_newdata_mpa3$predicted_tot25==0]),sum(fished_newdata_mpa3$predicted_tot25)-sum(fished_newdata_mpa3$predicted_tot50),sum(fished_newdata_mpa3$predicted_tot50)-sum(fished_newdata_mpa3$predicted_tot75),sum(fished_newdata_mpa3$predicted_tot75))
matrixtot[1:4,4]=rep("Initial", 4)
matrixtot[1:4,5]=c("75-100%","50-75%","25-50%","0-25%")
matrixtot[5:8,1]=rep(4,4)
matrixtot[5:8,2]=c(8,7,6,5)
matrixtot[5:8,4]=rep("reserve",4)
matrixtot[5:8,5]=c("0to0","0to25","0to50","0to75")
matrixtot[5:8,3]=c(length(fished_newdata_mpa3$unfished_changetot[fished_newdata_mpa3$unfished_changetot=="0to0"]),length(fished_newdata_mpa3$unfished_changetot[fished_newdata_mpa3$unfished_changetot=="0to25"]),length(fished_newdata_mpa3$unfished_changetot[fished_newdata_mpa3$unfished_changetot=="0to50"]),length(fished_newdata_mpa3$unfished_changetot[fished_newdata_mpa3$unfished_changetot=="0to75"]))
matrixtot[9:11,1]=rep(3,3)
matrixtot[9:11,2]=c(7,6,5)
matrixtot[9:11,4]=rep("reserve",3)
matrixtot[9:11,5]=c("25to25","25to50","25to75")
matrixtot[9:11,3]=c(length(fished_newdata_mpa3$unfished_changetot[fished_newdata_mpa3$unfished_changetot=="25to25"]),length(fished_newdata_mpa3$unfished_changetot[fished_newdata_mpa3$unfished_changetot=="25to50"]),length(fished_newdata_mpa3$unfished_changetot[fished_newdata_mpa3$unfished_changetot=="25to75"]))
matrixtot[12:13,1]=rep(2,2)
matrixtot[12:13,2]=c(6,5)
matrixtot[12:13,4]=rep("reserve",2)
matrixtot[12:13,5]=c("50to50","50to75")
matrixtot[12:13,3]=c(length(fished_newdata_mpa3$unfished_changetot[fished_newdata_mpa3$unfished_changetot=="50to50"]),length(fished_newdata_mpa3$unfished_changetot[fished_newdata_mpa3$unfished_changetot=="50to75"]))
matrixtot[14,1]=1
matrixtot[14,2]=5
matrixtot[14,4]="reserve"
matrixtot[14,5]="75to75"
matrixtot[14,3]=length(fished_newdata_mpa3$unfished_changetot[fished_newdata_mpa3$unfished_changetot=="75to75"])

#alluvial charts
tot=as.data.frame(matrixtot)
colnames(tot)=c("from","to","value","scenario","clarified")
tot$from=as.numeric(as.character(tot$from))
tot$to=as.numeric(as.character(tot$to))
tot$value=as.numeric(as.character(tot$value))
tot=tot[order(tot$from),]
nodes <- data.frame(c("75-100%","50-75%","25-50%","0-25%","75-100%","50-75%","25-50%","0-25%"))
names(nodes) <- "name"
nodes$group = as.factor(c("group1", "group2","group3", "group4", "group5", "group6", "group7", "group8"))
MR <- tot
link <- MR[,c(1:3)]
link <- link[-c(1:4),]
link$from <- link$from-1
link$to <- link$to-1
link$zeroout=ifelse(link$value==0,1,0)
link=link[link$zeroout==0,]
link$zeroout=NULL
link$group=ifelse(link$from==0,"group1", ifelse(link$from==1,"group2", ifelse(link$from==2,"group3", "group4")))
nodes$sites=c(round(MR[1:4,3],digits=3),round(sum(fished_newdata_mpa3$unfished_tot75),digits=3),round(sum(fished_newdata_mpa3$unfished_tot50[fished_newdata_mpa3$unfished_tot75==0]),digits=3),round(sum(fished_newdata_mpa3$unfished_tot25[fished_newdata_mpa3$unfished_tot50==0]),digits=3),round(length(fished_newdata_mpa3$unfished_tot25[fished_newdata_mpa3$unfished_tot25==0]),digits=3))
nodes$sites=c(round(MR[1:4,3]/length(fished_newdata_mpa3$unfished_tot75),digits=3)*100,round(sum(fished_newdata_mpa3$unfished_tot75)/length(fished_newdata_mpa3$unfished_tot75),digits=3)*100,round(sum(fished_newdata_mpa3$unfished_tot50[fished_newdata_mpa3$unfished_tot75==0])/length(fished_newdata_mpa3$unfished_tot75),digits=3)*100,round(sum(fished_newdata_mpa3$unfished_tot25[fished_newdata_mpa3$unfished_tot50==0])/length(fished_newdata_mpa3$unfished_tot75),digits=3)*100,round(length(fished_newdata_mpa3$unfished_tot25[fished_newdata_mpa3$unfished_tot25==0])/length(fished_newdata_mpa3$unfished_tot75),digits=3)*100)


networkMR_tot_mpa <- sankeyNetwork(Links = link, Nodes = nodes, Source = "from",
                                   Target = "to", Value = "value", NodeID = "sites",
                                   units = "", fontSize = 13, nodeWidth = 50,colourScale=my_color, NodeGroup = "group", LinkGroup = "group", width=300)


networkMR_tot_mpa

#save everything
#save.image(file='Cinneretal_submitted_afterrevisions.RData')
#load('Cinneretal_submitted_afterrevisions.RData')

