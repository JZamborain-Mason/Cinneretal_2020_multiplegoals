#Code for Cinner et al.2019: Meeting multiple goals for the world's coral reefs
#Created by Jessica Zamborain-Mason and David Moulliot on 15/08/2018
#ARC Centre of Excellence for Coral Reef Studies
#R version 3.4.2 (2017-09-28)

##Remove everything from the environment
rm(list = ls())


##set working directory
#NOTE:Please add your working directory (i.e., location where data is stored)
setwd("c:/users/Jessi/Documents/AUSTRALIA/PhD/Research Worker Josh 2018/Multifunctionality paper/R/FINAL CODE_DATA/FINAL") 
#setwd("") 

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


#load data
alldata<-read.csv("Cinneretal_2019_data_multipleconsgoals.csv",header=T) 

#mpa size and age restrictions
alldata$mpacondition=ifelse(alldata$Protection=="UnfishedHigh" &(alldata$MPAage3>4|alldata$MPAage3==4) & (alldata$NTZarea>2 |alldata$NTZarea==0),1,0)
alldata$keep=ifelse(!alldata$Protection=="UnfishedHigh",1,ifelse(alldata$Protection=="UnfishedHigh" &alldata$mpacondition==1,1,0))
#alldata=alldata[alldata$keep==1,] #Uncomment if you want to use more strick guidelines for the marine reserves used

##Explanatory variables: relevel categorical variables and standardize continuous variables

#Standardize function for continuous variables
standardize <- function(x){(x-mean(x, na.rm=T))/(2*sd(x, na.rm=T))} 

#reef-scale covariates
alldata$DepthCategory<-relevel(alldata$DepthCategory,ref="4-10m")
alldata$CleanHabitat<-relevel(alldata$CleanHabitat,ref="Slope")
alldata$Protection<-relevel(alldata$Protection,ref="Fished")
alldata$CensusMethod<-relevel(alldata$CensusMethod,ref="Standard belt transect")
alldata$sgrav_NC<-standardize(log(alldata$grav_NC+1))
alldata$sgrav_NP<-standardize(log(alldata$grav_NP+1))
alldata$sTotal_sampling_area<-standardize(log(alldata$Total_sampling_area))

#reef cluster-scale covariates
alldata$sOcean_prod<-standardize(log(alldata$Ocean_prod))
alldata$sClimate_stress<-standardize(alldata$Climate_stress)
alldata$sRegional_population_growth<-standardize(alldata$Regional_population_growth)


#nation/state-scale covariates
alldata$sReef_fish_landings_per_km2<-standardize(log(alldata$Reef_fish_landings_per_km2+1))
alldata$sLarger_pop_size<-standardize(log(alldata$Larger_pop_size+1))
alldata$sHDI=standardize(alldata$HDI)



##check that multicolinearity is not  a concern
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

#pairs correlations among covariates:all correlations below 0.6

pairs(~sgrav_NP+sgrav_NC+
        Protection+
        sRegional_population_growth+
        sHDI+
        sLarger_pop_size+
        sReef_fish_landings_per_km2+ 
        sOcean_prod+
        sClimate_stress+
        DepthCategory+
        CleanHabitat+
        sTotal_sampling_area+
        CensusMethod+Geographic_Basin, data= alldata,lower.panel=panel.cor )

#variance inflation factors  function for mixed models
vif.mer <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}
#variance inflation factor
VIF.table=as.data.frame(vif.mer(lmer(log(Biomass_above20cm+1)~
                                       DepthCategory+ CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NP+sgrav_NC+
                                       sRegional_population_growth+sOcean_prod+sClimate_stress+
                                       sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+(1|Geographic_Basin/Larger/ReefCluster), data=alldata)))
colnames(VIF.table)="VIF"

print(VIF.table)

##exploring and transforming response variables 

#Parrotfish scraping: continuous variable with almost 31% of the data containing 0's
#hurdle model: two part model. Part 1 estimates the probability of observing the function (binomial), and part 2 estimates the standardized effect size of the covaraites and the mean herbivore function given that we observed the function (normally distributed when log transformed: gaussian)
summary(alldata$Herbivory_function)# 136 of the total sites do not Parrotfish scraping data
length(alldata$Herbivory_function[alldata$Herbivory_function==0])/length(alldata$Herbivory_function[!is.na(alldata$Herbivory_function)])
hist(log(alldata$Herbivory_function+1))
#presence/absence ofParrotfish scraping (for part 1 of hurdle model)
alldata$PA_Herbivory_function= ifelse(alldata$Herbivory_function>0,1,0)
alldata$tHerbivory_function= log(alldata$Herbivory_function+1)
#filter non zero component of the data and log-transform Parrotfish scraping of non-zero data
nonzeroHFdata=alldata[alldata$Herbivory_function>0,]
nonzeroHFdata$tHerbivory_function=log(nonzeroHFdata$Herbivory_function)
hist(nonzeroHFdata$tHerbivory_function)
nonzeroHFdata=nonzeroHFdata[!is.na(nonzeroHFdata$tHerbivory_function),]
#Trait diversity: continuous variable with no zeros. 
alldata$Functional_entropy=alldata$FDis
summary(alldata$Functional_entropy)# 136 of the total sites do not have Trait diversity data
hist(alldata$Functional_entropy)


#REEF FISH BIOMASS ABOVE 20CM #only 3 % of the data are 0's. Consequently, we log+1 trasform it so it has a normal distribution
summary(alldata$Biomass_above20cm)
length(alldata$Biomass_above20cm[alldata$Biomass_above20cm==0])/length(alldata$Biomass_above20cm)
hist(log(alldata$Biomass_above20cm+1))
hist(sqrt(sqrt(alldata$Biomass_above20cm)))
alldata$tBiomass_above20cm=log(alldata$Biomass_above20cm+1)

#correlations among response variables
nonzero=alldata[alldata$Herbivory_function>0 & alldata$Biomass_above20cm>0,]

pairs(~tBiomass_above20cm+Functional_entropy+tHerbivory_function,data=nonzero,lower.panel=panel.cor, 
      pch = 21, bg = "darkgrey",labels=c("log(Biomass >20cm+1)","Trait diversity", "log(Parrotfish scraping+1)"),cex.labels=1.5,font.labels=2,diag.panel =panel_hist, hist.col="grey")

pairs(~tBiomass_above20cm+Functional_entropy+tHerbivory_function+MidSize,data=nonzero,lower.panel=panel.cor, 
      pch = 21, bg = "darkgrey",labels=c("log(Biomass >20cm+1)","Trait diversity", "log(Parrotfish scraping+1)","Length(cm)"),cex.labels=1.5,font.labels=2,diag.panel =panel_hist, hist.col="grey",main = "Non-zero data")



pairs(~tBiomass_above20cm+Functional_entropy+tHerbivory_function,data=alldata,lower.panel=panel.cor, 
      pch = 21, bg = "darkgrey",labels=c("log(Biomass >20cm+1)","Trait diversity", "log(Parrotfish scraping+1)"),cex.labels=1.5,font.labels=2,diag.panel =panel_hist, hist.col="grey")


pairs(~tBiomass_above20cm+Functional_entropy+tHerbivory_function+MidSize,data=alldata,lower.panel=panel.cor, 
      pch = 21, bg = "darkgrey",labels=c("log(Biomass >20cm+1)","Trait diversity", "log(Parrotfish scraping+1)","Length(cm)"),cex.labels=1.5,font.labels=2,diag.panel =panel_hist, hist.col="grey",main = "All data")

         
#correlations among explanatory variables
pairs(~DepthCategory+Geographic_Basin+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
        sRegional_population_growth+sOcean_prod+sClimate_stress+
        sHDI+sLarger_pop_size+sReef_fish_landings_per_km2,data=alldata,lower.panel=panel.cor, 
      pch = 21, bg = "darkgrey",labels=c("Depth Category","Basin","Habitat type", "Management
protection", "Sampling
method","Sampling area","Market gravity","Nearest
settlement
gravity","Local
population
growth","Ocean
productivity", "Climate stress
index","HDI","Population size",  "Reef fish
landings"),cex.labels=0.9,font.labels=2)

#countries included
countries=as.data.frame(unique(alldata$Larger))
#write.csv(countries, "countries.included.csv",row.names=F)


##Models

#Parrotfish scraping
#hurdle component 1: presence/absence
PA_her_model_NULL<-stan_glmer(PA_Herbivory_function~
                               (1|Larger/ReefCluster),
                              data=alldata[!is.na(alldata$PA_Herbivory_function),],family=binomial(link = "logit"),
                              iter=100000,  warmup=90000,
                              chains=4, 
                              thin=10) 

PA_her_model<-stan_glmer(PA_Herbivory_function~
                           DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                           sRegional_population_growth+sOcean_prod+sClimate_stress+
                           sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                           (1|Larger/ReefCluster),
                           data=alldata[!is.na(alldata$PA_Herbivory_function),],family=binomial(link = "logit"),
                           iter=100000,  warmup=90000,
                           chains=4, 
                           thin=10)
compare_models(loo(PA_her_model,PA_her_model_NULL))

#check model diagnostics

pp_check(PA_her_model)
stan_trace(PA_her_model)
stan_ac(PA_her_model)
plot(PA_her_model, 'mcmc_rhat_hist') 
plot(PA_her_model, 'mcmc_neff_hist') 

#check model fit

hist(resid(PA_her_model)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(PA_her_model), x=fitted(PA_her_model)))
#no heterasdicity
hist(posterior_predict(PA_her_model))
stan_hist(PA_her_model)
posteriorestimates_PA=as.data.frame(PA_her_model$stan_summary)

#hurdle component 2: non-zero data
her_model_NULL<-stan_glmer(tHerbivory_function~
                             (1|Larger/ReefCluster),
                           data=nonzeroHFdata,family=gaussian,
                           iter=100000,  warmup=90000,
                           chains=4, 
                           thin=10) 

her_model<-stan_glmer(tHerbivory_function~
                        DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                        sRegional_population_growth+sOcean_prod+sClimate_stress+
                        sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                        (1|Larger/ReefCluster),
                      data=nonzeroHFdata,family=gaussian,
                      iter=100000,  warmup=90000,
                      chains=4, 
                      thin=10) 

compare_models(loo(her_model_NULL),loo(her_model))

#check model diagnostics
pp_check(her_model)
stan_trace(her_model)
stan_ac(her_model)
plot(her_model, 'mcmc_rhat_hist') 
plot(her_model, 'mcmc_neff_hist') 

#check model fit

hist(resid(her_model)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(her_model), x=fitted(her_model)))
#no heterasdicity
hist(posterior_predict(her_model))
stan_hist(her_model)
posteriorestimates_hf=as.data.frame(her_model$stan_summary)
#write.csv(posteriorestimates_hf,"hf_posterior_allreserves.csv", row.names=F)
#check residuals among geographic basin
nonzeroHFdata$residhf=resid(her_model)

#check residuals do not show patterns with geographic basin
r1=ggplot(nonzeroHFdata, aes(x=Geographic_Basin,y=residhf))+geom_boxplot(fill="grey",col="black")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0),panel.background = element_rect(fill="white",colour="darkgrey"))+xlab("")+ylab("residuals Parrotfish Scraping")

#Trait diversity

nonzerotd=alldata[!is.na(alldata$Functional_entropy),]
fundiv_model_NULL<-stan_glmer(Functional_entropy~
                                (1|Larger/ReefCluster),
                              data=alldata,family=gaussian,
                              iter=100000,  warmup=90000,
                              chains=4, 
                              thin=10) 

fundiv_model<-stan_glmer(Functional_entropy~
                           DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                           sRegional_population_growth+sOcean_prod+sClimate_stress+
                           sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                           (1|Larger/ReefCluster),
                         data=alldata,family=gaussian,
                         iter=100000,  warmup=90000,
                         chains=4, 
                         thin=10) 
compare_models(loo(fundiv_model_NULL),loo(fundiv_model))

#check model diagnostics
pp_check(fundiv_model)
stan_trace(fundiv_model)
stan_ac(fundiv_model)
plot(fundiv_model, 'mcmc_rhat_hist') 
plot(fundiv_model, 'mcmc_neff_hist') 

#check model fit

hist(resid(fundiv_model)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(fundiv_model), x=fitted(fundiv_model)))
#no heterasdicity
hist(posterior_predict(fundiv_model))
stan_hist(fundiv_model)
posteriorestimates_fd=as.data.frame(fundiv_model$stan_summary)
#write.csv(posteriorestimates_fd,"fd_posterior_allreserves.csv", row.names=F)

#check residuals do not show patterns with geographic basin
nonzerotd$residtd=resid(fundiv_model)
r2=ggplot(nonzerotd, aes(x=Geographic_Basin,y=residtd))+geom_boxplot(fill="grey",col="black")+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust=0),panel.background = element_rect(fill="white",colour="darkgrey"))+xlab("")+ylab("residuals Trait Diversity")

#Biomass of reef fish >20cm
B20_model_NULL<-stan_glmer(lBiomass_above20cm~
                             (1|Larger/ReefCluster),
                           data=alldata,family=gaussian,
                           iter=100000,  warmup=90000,
                           chains=4, 
                           thin=10) 

B20_model<-stan_glmer(tBiomass_above20cm~
                        DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                        sRegional_population_growth+sOcean_prod+sClimate_stress+
                        sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                        (1|Larger/ReefCluster),
                      data=alldata,family=gaussian,
                      iter=100000,  warmup=90000,
                      chains=4, 
                      thin=10) 

compare_models(loo(B20_model_NULL),loo(B20_model))

#check model diagnostics
pp_check(B20_model)
stan_trace(B20_model)
stan_ac(B20_model)
plot(B20_model, 'mcmc_rhat_hist') 
plot(B20_model, 'mcmc_neff_hist') 

#check model fit

hist(resid(B20_model)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(B20_model), x=fitted(B20_model)))
#no heterasdicity
hist(posterior_predict(B20_model))
stan_hist(B20_model)
posteriorestimates_b20=as.data.frame(B20_model$stan_summary)
#write.csv(posteriorestimates_b20,"b20_posterior_allreserves.csv", row.names=F)

#check residuals do not show patterns with geographic basin
alldata$residb20=resid(B20_model)
r3=ggplot(alldata, aes(x=Geographic_Basin,y=residb20))+geom_boxplot(fill="grey",col="black")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0),panel.background = element_rect(fill="white",colour="darkgrey"))+xlab("")+ylab("residuals Biomass>20cm")

#Geographic Basin vs residuals for all three models
ggarrange(r1,r2,r3,nrow=1,ncol=3)

##figures##
#coefficient plots

#Parrotfish scraping presence/absence
coefplot_pa=posteriorestimates_PA[2:19,]
coefplot_pa$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Market gravity","Nearest settlement gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_pa$sign=ifelse(coefplot_pa$`2.5%`<0 & coefplot_pa$`97.5%`<0, "negative",ifelse(coefplot_pa$`2.5%`>0 & coefplot_pa$`97.5%`>0, "positive", "no effect"))
coefplot_pa$strength=ifelse(coefplot_pa$sign=="no effect", "open", "closed")
coefplot_pa$environmental=c(1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0)
coefplot_pa$order=c(9,10,11,12,13,2,1,16,17,18,8,7,6,14,15,3,4,5)
coefplot_pa[order(coefplot_pa$order),]
coefplot_pa$variable <- factor(coefplot_pa$variable, levels = coefplot_pa$variable[order(coefplot_pa$order)])
coefplot_pa_hf=
  ggplot(coefplot_pa,aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("log odds ratio")+xlab("")+ggtitle("")
#non-zero Parrotfish scraping
coefplot_hf=posteriorestimates_hf[2:19,]
coefplot_hf$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Market gravity","Nearest settlement gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_hf$sign=ifelse(coefplot_hf$`2.5%`<0 & coefplot_hf$`97.5%`<0, "negative",ifelse(coefplot_hf$`2.5%`>0 & coefplot_hf$`97.5%`>0, "positive", "no effect"))
coefplot_hf$strength=ifelse(coefplot_hf$sign=="no effect", "open", "closed")
coefplot_hf$environmental=c(1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0)
coefplot_hf$order=c(9,10,11,12,13,2,1,16,17,18,8,7,6,14,15,3,4,5)
coefplot_hf[order(coefplot_hf$order),]
coefplot_hf$variable <- factor(coefplot_hf$variable, levels = coefplot_hf$variable[order(coefplot_hf$order)])
coefplot1=
  ggplot(coefplot_hf[coefplot_hf$environmental==1,],aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("Standardized effect size")+xlab("")+ggtitle("Parrotfish scraping")+
  theme(axis.text.y=element_blank())
coefplot1b=
  ggplot(coefplot_hf[coefplot_hf$environmental==0,],aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("Parrotfish scraping")+
  theme(axis.text.y=element_blank())

#BIOMASS REEF FISH >20 CM
coefplot_B20=posteriorestimates_b20[2:19,]
coefplot_B20$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Market gravity","Nearest settlement gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_B20$sign=ifelse(coefplot_B20$`2.5%`<0 & coefplot_B20$`97.5%`<0, "negative",ifelse(coefplot_B20$`2.5%`>0 & coefplot_B20$`97.5%`>0, "positive", "no effect"))
coefplot_B20$strength=ifelse(coefplot_B20$sign=="no effect", "open", "closed")
coefplot_B20$environmental=c(1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0)
coefplot_B20$order=c(9,10,11,12,13,2,1,16,17,18,8,7,6,14,15,3,4,5)
coefplot_B20[order(coefplot_B20$order),]
coefplot_B20$variable <- factor(coefplot_B20$variable, levels = coefplot_B20$variable[order(coefplot_B20$order)])
coefplot2=
  ggplot(coefplot_B20[coefplot_B20$environmental==1,],aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("Biomass > 20cm")+theme_classic()
coefplot2b=
  ggplot(coefplot_B20[coefplot_B20$environmental==0,],aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("Biomass > 20cm")+theme_classic()


#Trait diversity
coefplot_fd=posteriorestimates_fd[2:19,]
coefplot_fd$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Market gravity","Nearest settlement gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_fd$sign=ifelse(coefplot_fd$`2.5%`<0 & coefplot_fd$`97.5%`<0, "negative",ifelse(coefplot_fd$`2.5%`>0 & coefplot_fd$`97.5%`>0, "positive", "no effect"))
coefplot_fd$strength=ifelse(coefplot_fd$sign=="no effect", "open", "closed")
coefplot_fd$environmental=c(1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0)
coefplot_fd$order=c(9,10,11,12,13,2,1,16,17,18,8,7,6,14,15,3,4,5)
coefplot_fd[order(coefplot_fd$order),]
coefplot_fd$variable <- factor(coefplot_fd$variable, levels = coefplot_fd$variable[order(coefplot_fd$order)])
coefplot3=
  ggplot(coefplot_fd[coefplot_fd$environmental==1,],aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("Trait diversity")+
  theme(axis.text.y=element_blank())
coefplot3b=
  ggplot(coefplot_fd[coefplot_fd$environmental==0,],aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("Trait diversity")+
  theme(axis.text.y=element_blank())
Fig.si2=ggarrange(coefplot2,coefplot1,coefplot3,nrow=1,ncol=3, widths=c(1.8,1,1))
Fig.si1=coefplot_pa_hf
Fig.2=ggarrange(coefplot2b,coefplot1b,coefplot3b,nrow=1,ncol=3, widths=c(1.9,1,1))


##maps of response variables
#jitter lat and long points
alldata$Site_Lat2=alldata$Site_Lat+runif(length(alldata$Site_Lat), min=0, max=4)
alldata$Site_Long2=alldata$Site_Long+runif(length(alldata$Site_Long), min=0, max=3)
alldata$Site_Lat2=ifelse(alldata$Site_Lat2>23.5, alldata$Site_Lat,alldata$Site_Lat2)
#extract map and match latitude
mapWorld <- map_data('world', wrap=c(-25,335), ylim=c(-55,75))
alldata$lon2 <- ifelse(alldata$Site_Long2 < -25, alldata$Site_Long2 + 360, alldata$Site_Long2) 
#specify gradient colours
mycolours = c("turquoise4", "darkturquoise", "white", "orangered3", "darkred") 

Fig.1a=ggplot() + 
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group),fill = "grey", color = "grey") +
  coord_fixed(xlim = c(30, 320),  ylim = c(30, -30), ratio = 1.3)+
  geom_point(data=alldata[!is.na(alldata$tHerbivory_function),], aes(x=lon2, y=Site_Lat2, fill=alldata$tHerbivory_function[!is.na(alldata$tHerbivory_function)]),colour="black", pch=21, size=3)+
  scale_fill_gradientn(colors = mycolours, name="log(Herbivory
                       function+1)", guide=F)+  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+scale_x_continuous("",breaks=c(80,180,280),labels=c(-100,0,100))+
  scale_y_continuous("Latitude",breaks=c(-20,-10,0,10,20))+ theme_classic()

Fig.1b=ggplot() + 
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group),fill = "grey", color = "grey") +
  coord_fixed(xlim = c(30, 320),  ylim = c(30, -30), ratio = 1.3)+
  geom_point(data=alldata, aes(x=lon2, y=Site_Lat2, fill=alldata$tBiomass_above20cm),colour="black", pch=21, size=3)+
  scale_fill_gradientn(colors = mycolours, name="log (Biomass
                       > 20 cm +1)", guide=F)+  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+scale_x_continuous("",breaks=c(80,180,280),labels=c(-100,0,100))+
  scale_y_continuous("",breaks=c(-20,-10,0,10,20))+ theme_classic()

Fig.1c=ggplot() + 
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group),fill = "grey", color = "grey") +
  coord_fixed(xlim = c(30, 320),  ylim = c(30, -30), ratio = 1.3)+
  geom_point(data=alldata[!is.na(alldata$Functional_entropy),], aes(x=lon2, y=Site_Lat2, fill=alldata$Functional_entropy[!is.na(alldata$Functional_entropy)] ),colour="black", pch=21, size=3)+
  scale_fill_gradientn(colors = mycolours, name="Trait
                       diversity", guide=F)+  geom_hline(yintercept =23.43695, lty=2)+
  geom_hline(yintercept =-23.43695, lty=2)+scale_x_continuous("Longitude",breaks=c(80,180,280),labels=c(-100,0,100))+
  scale_y_continuous("",breaks=c(-20,-10,0,10,20))+ theme_classic()

Fig.1=ggarrange(Fig.1b,Fig.1a,Fig.1c,nrow=3, ncol=1)


##################################################################################

# Metrics

# select only reef sites that contain all three metrics
alldata=alldata[!is.na(alldata$Functional_entropy)&!is.na(alldata$Herbivory_function),]

#create a matrix to store the thresholds data
res=matrix(NA,nrow=length(alldata$Functional_entropy),ncol=12)
colnames(res)=c("herb25","pred25","func25","tot25","herb50","pred50","func50","tot50","herb75","pred75","func75","tot75")
alldata=cbind(alldata,res)


# Parrotfish scraping
#benchmark
rherb=quantile(alldata$Herbivory_function,p=0.9, na.rm=T)#change p to 0.8 or 0.95 for different benckmarks used
#25%
alldata$herb25=alldata$Herbivory_function
alldata$herb25[alldata$herb25<=(rherb*0.25)]=0
alldata$herb25[alldata$herb25>(rherb*0.25)]=1
length(alldata$herb25[alldata$herb25==1])/length(alldata$herb25)
#50%
alldata$herb50=alldata$Herbivory_function
alldata$herb50[alldata$herb50<=(rherb*0.5)]=0
alldata$herb50[alldata$herb50>(rherb*0.5)]=1
length(alldata$herb50[alldata$herb50==1])/length(alldata$herb50)
#75%
alldata$herb75=alldata$Herbivory_function
alldata$herb75[alldata$herb75<=(rherb*3/4)]=0
alldata$herb75[alldata$herb75>(rherb*3/4)]=1
length(alldata$herb75[alldata$herb75==1])/length(alldata$herb75)
#set white theme for plots
white_theme<-theme(axis.text=element_text(colour="black"),
                   axis.ticks=element_line(colour="black"),
                   panel.grid.minor=element_blank(),
                   panel.background=element_rect(fill="white",colour="black"),
                   legend.justification=c(1,0),legend.position=c(.95, .05),line= element_blank(),
                   plot.background=element_rect(fill="transparent",colour=NA))



#Biomass above 20cm
#benchmark
rpred=quantile(alldata$Biomass_above20cm,p=0.9, na.rm=T) #change p to 0.8 or 0.95 for different benckmarks
#25%
alldata$pred25=alldata$Biomass_above20cm
alldata$pred25[alldata$pred25<=(rpred*0.25)]=0
alldata$pred25[alldata$pred25>(rpred*0.25)]=1
length(alldata$pred25[alldata$pred25==1])/length(alldata$pred25)
#50%
alldata$pred50=alldata$Biomass_above20cm
alldata$pred50[alldata$pred50<=(rpred*0.5)]=0
alldata$pred50[alldata$pred50>(rpred*0.5)]=1
length(alldata$pred50[alldata$pred50==1])/length(alldata$pred50)

#75%
alldata$pred75=alldata$Biomass_above20cm
alldata$pred75[alldata$pred75<=(rpred*3/4)]=0
alldata$pred75[alldata$pred75>(rpred*3/4)]=1
length(alldata$pred75[alldata$pred75==1])/length(alldata$pred75)

# Trait diversity
#benchmark
rfunc=quantile(alldata$Functional_entropy,p=0.9, na.rm=T)#change p to 0.8 or 0.95 for different benckmarks
#25%
alldata$func25=alldata$Functional_entropy
alldata$func25[alldata$func25<=(rfunc*0.25)]=0
alldata$func25[alldata$func25>(rfunc*0.25)]=1
length(alldata$func25[alldata$func25==1])/length(alldata$func25)

#50%
alldata$func50=alldata$Functional_entropy
alldata$func50[alldata$func50<=(rfunc*0.5)]=0
alldata$func50[alldata$func50>(rfunc*0.5)]=1
length(alldata$func50[alldata$func50==1])/length(alldata$func50)
#75%
alldata$func75=alldata$Functional_entropy
alldata$func75[alldata$func75<=(rfunc*3/4)]=0
alldata$func75[alldata$func75>(rfunc*3/4)]=1
length(alldata$func75[alldata$func75==1])/length(alldata$func75)

# Three metrics combined
#25%
alldata$tot25=alldata$func25+alldata$pred25+alldata$herb25
#50%
alldata$tot50=alldata$func50+alldata$pred50+alldata$herb50
#75%
alldata$tot75=alldata$func75+alldata$pred75+alldata$herb75


#visual representation of thresholds and benckmarks
a=ggplot(NULL)+geom_histogram(aes(x=alldata$Herbivory_function), col="black", fill="grey", alpha=0.5)+geom_vline(xintercept=rherb, lwd=2)+geom_vline(xintercept=rherb*0.25, lty=2)+geom_vline(xintercept=rherb*0.5, lty=2,col="blue")+geom_vline(xintercept=rherb*0.75, lty=2,col="red")+xlab("Parrotfish scraping")+xlim(c(0,2000))
b=ggplot(NULL)+geom_histogram(aes(x=alldata$Biomass_above20cm), col="black", fill="grey", alpha=0.5)+geom_vline(xintercept=rpred, lwd=2)+geom_vline(xintercept=rpred*0.25, lty=2)+geom_vline(xintercept=rpred*0.5, lty=2,col="blue")+geom_vline(xintercept=rpred*0.75, lty=2,col="red")+xlab("Biomass >20 cm")+xlim(c(0,10000))
c=ggplot(NULL)+geom_histogram(aes(x=alldata$Functional_entropy), col="black", fill="grey", alpha=0.5)+geom_vline(xintercept=rfunc, lwd=2)+geom_vline(xintercept=rfunc*0.25, lty=2)+geom_vline(xintercept=rfunc*0.5, lty=2,col="blue")+geom_vline(xintercept=rfunc*0.75, lty=2,col="red")+xlab("Trait diversity")

ggarrange(b,a,c,ncol=3,nrow=1)

#openly fished sites:proportion
fished=alldata[alldata$Protection=="Fished",]
length(fished$tot75[fished$tot75==3])/length(fished$tot75)
length(fished$tot50[fished$tot50==3])/length(fished$tot50)
length(fished$tot25[fished$tot25==3])/length(fished$tot25)

#binary data for reefs that have all functions and proportions
alldata$q25<-as.factor(ifelse(alldata$tot25>=3,1,0))
length(alldata$q25[alldata$q25==1])/length(alldata$q25)
alldata$q50<-as.factor(ifelse(alldata$tot50>=3,1,0))
length(alldata$q50[alldata$q50==1])/length(alldata$q50)
alldata$q75<-as.factor(ifelse(alldata$tot75>=3,1,0))
length(alldata$q75[alldata$q75==1])/length(alldata$q75)

#ALL THREE RESPONSE VARIABLES
#Model 25
mod25_3=stan_glmer(q25~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                   (1|Larger/ReefCluster),family = binomial(link=logit),data=alldata,control = list(adapt_delta = 0.999,max_treedepth = 15))

pp_check(mod25_3)
stan_trace(mod25_3)
plot(mod25_3, 'mcmc_rhat_hist')  

#check model fit

hist(resid(mod25_3)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(mod25_3), x=fitted(mod25_3)))

posteriorestimates_3_25=as.data.frame(mod25_3$stan_summary)
#write.csv(posteriorestimates_3_25,"posteriorestimates_3_25.csv")

#Model 50
mod50_3=stan_glmer(q50~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                   (1|Larger/ReefCluster),family = binomial(link=logit),data=alldata,control = list(adapt_delta = 0.999,max_treedepth = 15))

pp_check(mod50_3)
stan_trace(mod50_3)
plot(mod50_3, 'mcmc_rhat_hist')  

#check model fit
hist(resid(mod50_3)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(mod50_3), x=fitted(mod50_3)))
posteriorestimates_3_50=as.data.frame(mod50_3$stan_summary)
#write.csv(posteriorestimates_3_50,"posteriorestimates_3_50.csv")

#Model 75
mod75_3=stan_glmer(q75~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                   (1|Larger/ReefCluster),family = binomial(link=logit),data=alldata,control = list(adapt_delta = 0.999,max_treedepth = 15))
pp_check(mod75_3)
stan_trace(mod75_3)
plot(mod75_3, 'mcmc_rhat_hist')  

#check model fit
hist(resid(mod75_3)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(mod75_3), x=fitted(mod75_3)))

posteriorestimates_3_75=as.data.frame(mod75_3$stan_summary)
#write.csv(posteriorestimates_3_75,"posteriorestimates_3_75.csv")

#coefplots with same structure as manuscript
coefplot_3_25=posteriorestimates_3_25[2:19,]
coefplot_3_25$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Market gravity","Nearest settlement gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_3_25$sign=ifelse(coefplot_3_25$`2.5%`<0 & coefplot_3_25$`97.5%`<0, "negative",ifelse(coefplot_3_25$`2.5%`>0 & coefplot_3_25$`97.5%`>0, "positive", "no effect"))
coefplot_3_25$strength=ifelse(coefplot_3_25$sign=="no effect", "open", "closed")
coefplot_3_25$environmental=c(1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0)
coefplot_3_25$order=c(9,10,11,12,13,2,1,16,17,18,8,7,6,14,15,3,4,5)
coefplot_3_25[order(coefplot_3_25$order),]
coefplot_3_25$variable <- factor(coefplot_3_25$variable, levels = coefplot_3_25$variable[order(coefplot_3_25$order)])
coefplot_3_25_f=
  ggplot(coefplot_3_25,aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("25")+ theme(axis.text.y=element_blank())

coefplot_3_50=posteriorestimates_3_50[2:19,]
coefplot_3_50$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Market gravity","Nearest settlement gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_3_50$sign=ifelse(coefplot_3_50$`2.5%`<0 & coefplot_3_50$`97.5%`<0, "negative",ifelse(coefplot_3_50$`2.5%`>0 & coefplot_3_50$`97.5%`>0, "positive", "no effect"))
coefplot_3_50$strength=ifelse(coefplot_3_50$sign=="no effect", "open", "closed")
coefplot_3_50$environmental=c(1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0)
coefplot_3_50$order=c(9,10,11,12,13,2,1,16,17,18,8,7,6,14,15,3,4,5)
coefplot_3_50[order(coefplot_3_50$order),]
coefplot_3_50$variable <- factor(coefplot_3_50$variable, levels = coefplot_3_50$variable[order(coefplot_3_50$order)])
coefplot_3_50_f=
  ggplot(coefplot_3_50,aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("50")+ theme(axis.text.y=element_blank())

coefplot_3_75=posteriorestimates_3_75[2:19,]
coefplot_3_75$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Market gravity","Nearest settlement gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_3_75$sign=ifelse(coefplot_3_75$`2.5%`<0 & coefplot_3_75$`97.5%`<0, "negative",ifelse(coefplot_3_75$`2.5%`>0 & coefplot_3_75$`97.5%`>0, "positive", "no effect"))
coefplot_3_75$strength=ifelse(coefplot_3_75$sign=="no effect", "open", "closed")
coefplot_3_75$environmental=c(1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0)
coefplot_3_75$order=c(9,10,11,12,13,2,1,16,17,18,8,7,6,14,15,3,4,5)
coefplot_3_75[order(coefplot_3_75$order),]
coefplot_3_75$variable <- factor(coefplot_3_75$variable, levels = coefplot_3_75$variable[order(coefplot_3_75$order)])
coefplot_3_75_f=
  ggplot(coefplot_3_75,aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("75")+ theme(axis.text.y=element_blank())

coefplot3func_q=ggarrange(coefplot_3_25_f,coefplot_3_50_f,coefplot_3_75_f, nrow=1,ncol=3)

#Create modelled data
MyData_3<-expand.grid(sgrav_NC=seq(range(alldata$sgrav_NC)[[1]],range(alldata$sgrav_NC)[[2]],length=101),
                    Protection=c("Fished","Restricted","UnfishedHigh"),
                    DepthCategory=c(">10m"," 0-4m","4-10m"),
                    CleanHabitat=c("Crest","Flat","Lagoon_Back reef","Slope"),
                    sTotal_sampling_area=0,
                    CensusMethod=c("Distance sampling","Point intercept","Standard belt transect"),
                    sRegional_population_growth=0,
                    sOcean_prod=0,
                    sgrav_NP=0,
                    sHDI=0,
                    sReef_fish_landings_per_km2=0,
                    sClimate_stress=0,
                    sLarger_pop_size=0)

X <- model.matrix(~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                    sRegional_population_growth+sOcean_prod+sClimate_stress+
                    sHDI+sLarger_pop_size+sReef_fish_landings_per_km2,data=MyData_3)

coefs25_3=as.matrix(mod25_3)
coefs25_3=coefs25_3[,c("(Intercept)","DepthCategory>10m","DepthCategory0-4m","CleanHabitatCrest","CleanHabitatFlat","CleanHabitatLagoon_Back reef","ProtectionRestricted","ProtectionUnfishedHigh",
                   "CensusMethodDistance sampling","CensusMethodPoint intercept","sTotal_sampling_area","sgrav_NC","sgrav_NP","sRegional_population_growth","sOcean_prod","sClimate_stress",                                                                      
                   "sHDI","sLarger_pop_size","sReef_fish_landings_per_km2" )]
fit25_3=coefs25_3 %*% t(X)
fit25_3= exp(fit25_3) / (1 + exp(fit25_3))

MyData_3=MyData_3 %>% cbind(tidyMCMC(as.mcmc(fit25_3),
                                 conf.int=T, conf.method='HPDinterval'))
setnames(MyData_3, old = c("estimate","std.error","conf.low","conf.high"), new = c("estimate25","std.error25","conf.low25","conf.high25"))

coefs50_3=as.matrix(mod50_3)
coefs50_3=coefs50_3[,c("(Intercept)","DepthCategory>10m","DepthCategory0-4m","CleanHabitatCrest","CleanHabitatFlat","CleanHabitatLagoon_Back reef","ProtectionRestricted","ProtectionUnfishedHigh",
                   "CensusMethodDistance sampling","CensusMethodPoint intercept","sTotal_sampling_area","sgrav_NC","sgrav_NP","sRegional_population_growth","sOcean_prod","sClimate_stress",                                                                      
                   "sHDI","sLarger_pop_size","sReef_fish_landings_per_km2" )]
fit50_3=coefs50_3 %*% t(X)
fit50_3= exp(fit50_3) / (1 + exp(fit50_3))

MyData_3=MyData_3%>% cbind(tidyMCMC(as.mcmc(fit50_3),
                                 conf.int=T, conf.method='HPDinterval'))
setnames(MyData_3, old = c("estimate","std.error","conf.low","conf.high"), new = c("estimate50","std.error50","conf.low50","conf.high50"))

coefs75_3=as.matrix(mod75_3)
coefs75_3=coefs75_3[,c("(Intercept)","DepthCategory>10m","DepthCategory0-4m","CleanHabitatCrest","CleanHabitatFlat","CleanHabitatLagoon_Back reef","ProtectionRestricted","ProtectionUnfishedHigh",
                   "CensusMethodDistance sampling","CensusMethodPoint intercept","sTotal_sampling_area","sgrav_NC","sgrav_NP","sRegional_population_growth","sOcean_prod","sClimate_stress",                                                                      
                   "sHDI","sLarger_pop_size","sReef_fish_landings_per_km2" )]
fit75_3=coefs75_3 %*% t(X)
fit75_3= exp(fit75_3) / (1 + exp(fit75_3))

MyData_3=MyData_3 %>% cbind(tidyMCMC(as.mcmc(fit75_3),
                                 conf.int=T, conf.method='HPDinterval'))
setnames(MyData_3, old = c("estimate","std.error","conf.low","conf.high"), new = c("estimate75","std.error75","conf.low75","conf.high75"))
#write.csv(MyData_3,"MyData_3_allreserves.csv",row.names=F)

# Create Delta
dat_Fished_3 <- MyData_3[which(MyData_3$CleanHabitat=="Slope" & MyData_3$Protection=="Fished" & MyData_3$CensusMethod=="Standard belt transect" &
                             MyData_3$DepthCategory=="4-10m"),]
dat_Reserve_3 <- MyData_3[which(MyData_3$CleanHabitat=="Slope" & MyData_3$Protection=="UnfishedHigh" & MyData_3$CensusMethod=="Standard belt transect" &
                              MyData_3$DepthCategory=="4-10m"),]
dat_Rest_3 <- MyData_3[which(MyData_3$CleanHabitat=="Slope" & MyData_3$Protection=="Restricted" & MyData_3$CensusMethod=="Standard belt transect" &
                           MyData_3$DepthCategory=="4-10m"),]

D1 <- (dat_Reserve_3[,"estimate25"] - dat_Fished_3[,"estimate25"])
D2 <- (dat_Reserve_3[,"estimate50"] - dat_Fished_3[,"estimate50"])
D3 <- (dat_Reserve_3[,"estimate75"] - dat_Fished_3[,"estimate75"])

delta_FR_3 <- cbind(dat_Fished_3,D1,D2,D3)

D1 <- (dat_Rest_3[,"estimate25"] - dat_Fished_3[,"estimate25"])
D2 <- (dat_Rest_3[,"estimate50"] - dat_Fished_3[,"estimate50"])
D3 <- (dat_Rest_3[,"estimate75"] - dat_Fished_3[,"estimate75"])

delta_FRest_3 <- cbind(dat_Fished_3,D1,D2,D3)

D1 <- (dat_Reserve_3[,"estimate25"] - dat_Rest_3[,"estimate25"])
D2 <- (dat_Reserve_3[,"estimate50"] - dat_Rest_3[,"estimate50"])
D3 <- (dat_Reserve_3[,"estimate75"] - dat_Rest_3[,"estimate75"])

delta_RestR_3 <- cbind(dat_Fished_3,D1,D2,D3)



# figures for manuscript
p5_3mrf <- ggplot(delta_FR_3)+
  geom_line(aes(sgrav_NC,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_NC,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_NC,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.15,0.5),labels=NULL) + white_theme 

p6_3rf <- ggplot(delta_FRest_3)+
  geom_line(aes(sgrav_NC,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_NC,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_NC,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.justification=c(1,0), legend.position=c(0,1), line= element_blank(), axis.ticks = element_line()) +
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.15,0.5),labels=NULL) + white_theme

p7_3rr <- ggplot(delta_RestR_3)+
  geom_line(aes(sgrav_NC,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_NC,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_NC,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.justification=c(1,0), legend.position=c(0,1), line= element_blank(), axis.ticks = element_line()) +
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.15,0.5),labels=NULL) + white_theme

p1_3fci <- ggplot(MyData_3[which(MyData_3$CleanHabitat=="Slope" & MyData_3$Protection=="Fished" & MyData_3$CensusMethod=="Standard belt transect" &
                                 MyData_3$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_NC,estimate25),colour="#fbb4b9",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high25,ymin=conf.low25),fill="#fbb4b9",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate50),colour="#f768a1",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high50,ymin=conf.low50),fill="#f768a1",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high75,ymin=conf.low75),fill="#7a0177",alpha=0.3)+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + white_theme

p2_3mrci <- ggplot(MyData_3[which(MyData_3$CleanHabitat=="Slope" & MyData_3$Protection=="UnfishedHigh" & MyData_3$CensusMethod=="Standard belt transect" &
                                  MyData_3$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_NC,estimate25),colour="#fbb4b9",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high25,ymin=conf.low25),fill="#fbb4b9",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate50),colour="#f768a1",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high50,ymin=conf.low50),fill="#f768a1",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high75,ymin=conf.low75),fill="#7a0177",alpha=0.3)+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + white_theme


p3_3rci <- ggplot(MyData_3[which(MyData_3$CleanHabitat=="Slope" & MyData_3$Protection=="Restricted" & MyData_3$CensusMethod=="Standard belt transect" &
                                 MyData_3$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_NC,estimate25),colour="#fbb4b9",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high25,ymin=conf.low25),fill="#fbb4b9",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate50),colour="#f768a1",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high50,ymin=conf.low50),fill="#f768a1",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high75,ymin=conf.low75),fill="#7a0177",alpha=0.3)+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + white_theme




# Parrotfish scraping
alldata$herb25=as.factor(alldata$herb25)
alldata$herb50=as.factor(alldata$herb50)
alldata$herb75=as.factor(alldata$herb75)
#Model 25
mod25_hf=stan_glmer(herb25~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                   (1|Larger/ReefCluster),family = binomial(link=logit),data=alldata,control = list(adapt_delta = 0.999,max_treedepth = 15))

pp_check(mod25_hf)
stan_trace(mod25_hf)
plot(mod25_hf, 'mcmc_rhat_hist')  

#check model fit
hist(resid(mod25_hf)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(mod25_hf), x=fitted(mod25_hf)))
posteriorestimates_hf_25=as.data.frame(mod25_hf$stan_summary)
#write.csv(posteriorestimates_hf_25,"posteriorestimates_hf_25_allreserves.csv")

#Model 50
mod50_hf=stan_glmer(herb50~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                   (1|Larger/ReefCluster),family = binomial(link=logit),data=alldata,control = list(adapt_delta = 0.999,max_treedepth = 15))

pp_check(mod50_hf)
stan_trace(mod50_hf)
plot(mod50_hf, 'mcmc_rhat_hist')  

#check model fit
hist(resid(mod50_hf)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(mod50_hf), x=fitted(mod50_hf)))
posteriorestimates_hf_50=as.data.frame(mod50_hf$stan_summary)
#write.csv(posteriorestimates_hf_50,"posteriorestimates_hf_50_allreserves.csv")

#Model 75
mod75_hf=stan_glmer(herb75~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                   (1|Larger/ReefCluster),family = binomial(link=logit),data=alldata,control = list(adapt_delta = 0.999,max_treedepth = 15))
pp_check(mod75_hf)
stan_trace(mod75_hf)
plot(mod75_hf, 'mcmc_rhat_hist')  

#check model fit
hist(resid(mod75_hf)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(mod75_hf), x=fitted(mod75_hf)))
posteriorestimates_hf_75=as.data.frame(mod75_hf$stan_summary)
#write.csv(posteriorestimates_hf_75,"posteriorestimates_hf_75_allreserves.csv")



#coefplots with same structure as manuscript
coefplot_hf_25=posteriorestimates_hf_25[2:19,]
coefplot_hf_25$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Market gravity","Nearest settlement gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_hf_25$sign=ifelse(coefplot_hf_25$`2.5%`<0 & coefplot_hf_25$`97.5%`<0, "negative",ifelse(coefplot_hf_25$`2.5%`>0 & coefplot_hf_25$`97.5%`>0, "positive", "no effect"))
coefplot_hf_25$strength=ifelse(coefplot_hf_25$sign=="no effect", "open", "closed")
coefplot_hf_25$environmental=c(1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0)
coefplot_hf_25$order=c(9,10,11,12,13,2,1,16,17,18,8,7,6,14,15,3,4,5)
coefplot_hf_25[order(coefplot_hf_25$order),]
coefplot_hf_25$variable <- factor(coefplot_hf_25$variable, levels = coefplot_hf_25$variable[order(coefplot_hf_25$order)])
coefplot_hf_25_f=
  ggplot(coefplot_hf_25,aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("25")+ theme(axis.text.y=element_blank())

coefplot_hf_50=posteriorestimates_hf_50[2:19,]
coefplot_hf_50$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Market gravity","Nearest settlement gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_hf_50$sign=ifelse(coefplot_hf_50$`2.5%`<0 & coefplot_hf_50$`97.5%`<0, "negative",ifelse(coefplot_hf_50$`2.5%`>0 & coefplot_hf_50$`97.5%`>0, "positive", "no effect"))
coefplot_hf_50$strength=ifelse(coefplot_hf_50$sign=="no effect", "open", "closed")
coefplot_hf_50$environmental=c(1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0)
coefplot_hf_50$order=c(9,10,11,12,13,2,1,16,17,18,8,7,6,14,15,3,4,5)
coefplot_hf_50[order(coefplot_hf_50$order),]
coefplot_hf_50$variable <- factor(coefplot_hf_50$variable, levels = coefplot_hf_50$variable[order(coefplot_hf_50$order)])
coefplot_hf_50_f=
  ggplot(coefplot_hf_50,aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("50")+ theme(axis.text.y=element_blank())

coefplot_hf_75=posteriorestimates_hf_75[2:19,]
coefplot_hf_75$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Market gravity","Nearest settlement gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_hf_75$sign=ifelse(coefplot_hf_75$`2.5%`<0 & coefplot_hf_75$`97.5%`<0, "negative",ifelse(coefplot_hf_75$`2.5%`>0 & coefplot_hf_75$`97.5%`>0, "positive", "no effect"))
coefplot_hf_75$strength=ifelse(coefplot_hf_75$sign=="no effect", "open", "closed")
coefplot_hf_75$environmental=c(1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0)
coefplot_hf_75$order=c(9,10,11,12,13,2,1,16,17,18,8,7,6,14,15,3,4,5)
coefplot_hf_75[order(coefplot_hf_75$order),]
coefplot_hf_75$variable <- factor(coefplot_hf_75$variable, levels = coefplot_hf_75$variable[order(coefplot_hf_75$order)])
coefplot_hf_75_f=
  ggplot(coefplot_hf_75,aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("75")+ theme(axis.text.y=element_blank())

coefplothf_q=ggarrange(coefplot_hf_25_f,coefplot_hf_50_f,coefplot_hf_75_f, nrow=1,ncol=3)

#Create modeled data
MyData_hf<-expand.grid(sgrav_NC=seq(range(alldata$sgrav_NC)[[1]],range(alldata$sgrav_NC)[[2]],length=101),
                    Protection=c("Fished","Restricted","UnfishedHigh"),
                    DepthCategory=c(">10m"," 0-4m","4-10m"),
                    CleanHabitat=c("Crest","Flat","Lagoon_Back reef","Slope"),
                    sTotal_sampling_area=0,
                    CensusMethod=c("Distance sampling","Point intercept","Standard belt transect"),
                    sRegional_population_growth=0,
                    sOcean_prod=0,
                    sgrav_NP=0,
                    sHDI=0,
                    sReef_fish_landings_per_km2=0,
                    sClimate_stress=0,
                    sLarger_pop_size=0)

X <- model.matrix(~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                    sRegional_population_growth+sOcean_prod+sClimate_stress+
                    sHDI+sLarger_pop_size+sReef_fish_landings_per_km2,data=MyData_hf)

coefs25_hf=as.matrix(mod25_hf)
coefs25_hf=coefs25_hf[,c("(Intercept)","DepthCategory>10m","DepthCategory0-4m","CleanHabitatCrest","CleanHabitatFlat","CleanHabitatLagoon_Back reef","ProtectionRestricted","ProtectionUnfishedHigh",
                   "CensusMethodDistance sampling","CensusMethodPoint intercept","sTotal_sampling_area","sgrav_NC","sgrav_NP","sRegional_population_growth","sOcean_prod","sClimate_stress",                                                                      
                   "sHDI","sLarger_pop_size","sReef_fish_landings_per_km2" )]
fit25_hf=coefs25_hf %*% t(X)
fit25_hf= exp(fit25_hf) / (1 + exp(fit25_hf))

MyData_hf=MyData_hf %>% cbind(tidyMCMC(as.mcmc(fit25_hf),
                                 conf.int=T, conf.method='HPDinterval'))
setnames(MyData_hf, old = c("estimate","std.error","conf.low","conf.high"), new = c("estimate25","std.error25","conf.low25","conf.high25"))

coefs50_hf=as.matrix(mod50_hf)
coefs50_hf=coefs50_hf[,c("(Intercept)","DepthCategory>10m","DepthCategory0-4m","CleanHabitatCrest","CleanHabitatFlat","CleanHabitatLagoon_Back reef","ProtectionRestricted","ProtectionUnfishedHigh",
                   "CensusMethodDistance sampling","CensusMethodPoint intercept","sTotal_sampling_area","sgrav_NC","sgrav_NP","sRegional_population_growth","sOcean_prod","sClimate_stress",                                                                      
                   "sHDI","sLarger_pop_size","sReef_fish_landings_per_km2" )]
fit50_hf=coefs50_hf %*% t(X)
fit50_hf= exp(fit50_hf) / (1 + exp(fit50_hf))

MyData_hf=MyData_hf %>% cbind(tidyMCMC(as.mcmc(fit50_hf),
                                 conf.int=T, conf.method='HPDinterval'))
setnames(MyData_hf, old = c("estimate","std.error","conf.low","conf.high"), new = c("estimate50","std.error50","conf.low50","conf.high50"))

coefs75_hf=as.matrix(mod75_hf)
coefs75_hf=coefs75_hf[,c("(Intercept)","DepthCategory>10m","DepthCategory0-4m","CleanHabitatCrest","CleanHabitatFlat","CleanHabitatLagoon_Back reef","ProtectionRestricted","ProtectionUnfishedHigh",
                   "CensusMethodDistance sampling","CensusMethodPoint intercept","sTotal_sampling_area","sgrav_NC","sgrav_NP","sRegional_population_growth","sOcean_prod","sClimate_stress",                                                                      
                   "sHDI","sLarger_pop_size","sReef_fish_landings_per_km2" )]
fit75_hf=coefs75_hf %*% t(X)
fit75_hf= exp(fit75_hf) / (1 + exp(fit75_hf))

MyData_hf=MyData_hf %>% cbind(tidyMCMC(as.mcmc(fit75_hf),
                                 conf.int=T, conf.method='HPDinterval'))
setnames(MyData_hf, old = c("estimate","std.error","conf.low","conf.high"), new = c("estimate75","std.error75","conf.low75","conf.high75"))
#write.csv(MyData_hf,"MyData_hf_allreserves.csv",row.names=F)

# Create Delta
dat_Fished_hf <- MyData_hf[which(MyData_hf$CleanHabitat=="Slope" & MyData_hf$Protection=="Fished" & MyData_hf$CensusMethod=="Standard belt transect" &
                             MyData_hf$DepthCategory=="4-10m"),]
dat_Reserve_hf <- MyData_hf[which(MyData_hf$CleanHabitat=="Slope" & MyData_hf$Protection=="UnfishedHigh" & MyData_hf$CensusMethod=="Standard belt transect" &
                              MyData_hf$DepthCategory=="4-10m"),]
dat_Rest_hf <- MyData_hf[which(MyData_hf$CleanHabitat=="Slope" & MyData_hf$Protection=="Restricted" & MyData_hf$CensusMethod=="Standard belt transect" &
                           MyData_hf$DepthCategory=="4-10m"),]

D1 <- (dat_Reserve_hf[,"estimate25"] - dat_Fished_hf[,"estimate25"])
D2 <- (dat_Reserve_hf[,"estimate50"] - dat_Fished_hf[,"estimate50"])
D3 <- (dat_Reserve_hf[,"estimate75"] - dat_Fished_hf[,"estimate75"])

delta_FR_hf <- cbind(dat_Fished_hf,D1,D2,D3)

D1 <- (dat_Rest_hf[,"estimate25"] - dat_Fished_hf[,"estimate25"])
D2 <- (dat_Rest_hf[,"estimate50"] - dat_Fished_hf[,"estimate50"])
D3 <- (dat_Rest_hf[,"estimate75"] - dat_Fished_hf[,"estimate75"])

delta_FRest_hf <- cbind(dat_Fished_hf,D1,D2,D3)

D1 <- (dat_Reserve_hf[,"estimate25"] - dat_Rest_hf[,"estimate25"])
D2 <- (dat_Reserve_hf[,"estimate50"] - dat_Rest_hf[,"estimate50"])
D3 <- (dat_Reserve_hf[,"estimate75"] - dat_Rest_hf[,"estimate75"])

delta_RestR_hf <- cbind(dat_Fished_hf,D1,D2,D3)


#for manuscript
p5_hmrf <- ggplot(delta_FR_hf)+
  geom_line(aes(sgrav_NC,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_NC,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_NC,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.15,0.5),labels=NULL) + white_theme 

p6_hrf <- ggplot(delta_FRest_hf)+
  geom_line(aes(sgrav_NC,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_NC,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_NC,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.justification=c(1,0), legend.position=c(0,1), line= element_blank(), axis.ticks = element_line()) +
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.15,0.5),labels=NULL) + white_theme

p7_hrr <- ggplot(delta_RestR_hf)+
  geom_line(aes(sgrav_NC,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_NC,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_NC,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.justification=c(1,0), legend.position=c(0,1), line= element_blank(), axis.ticks = element_line()) +
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.15,0.5),labels=NULL) + white_theme

p1_hfci <- ggplot(MyData_hf[which(MyData_hf$CleanHabitat=="Slope" & MyData_hf$Protection=="Fished" & MyData_hf$CensusMethod=="Standard belt transect" &
                                 MyData_hf$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_NC,estimate25),colour="#fbb4b9",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high25,ymin=conf.low25),fill="#fbb4b9",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate50),colour="#f768a1",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high50,ymin=conf.low50),fill="#f768a1",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high75,ymin=conf.low75),fill="#7a0177",alpha=0.3)+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + white_theme

p2_hmrci <- ggplot(MyData_hf[which(MyData_hf$CleanHabitat=="Slope" & MyData_hf$Protection=="UnfishedHigh" & MyData_hf$CensusMethod=="Standard belt transect" &
                                  MyData_hf$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_NC,estimate25),colour="#fbb4b9",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high25,ymin=conf.low25),fill="#fbb4b9",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate50),colour="#f768a1",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high50,ymin=conf.low50),fill="#f768a1",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high75,ymin=conf.low75),fill="#7a0177",alpha=0.3)+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + white_theme


p3_hrci <- ggplot(MyData_hf[which(MyData_hf$CleanHabitat=="Slope" & MyData_hf$Protection=="Restricted" & MyData_hf$CensusMethod=="Standard belt transect" &
                                 MyData_hf$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_NC,estimate25),colour="#fbb4b9",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high25,ymin=conf.low25),fill="#fbb4b9",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate50),colour="#f768a1",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high50,ymin=conf.low50),fill="#f768a1",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high75,ymin=conf.low75),fill="#7a0177",alpha=0.3)+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + white_theme


# Reef fish biomass above 20cm
alldata$pred25=as.factor(alldata$pred25)
alldata$pred50=as.factor(alldata$pred50)
alldata$pred75=as.factor(alldata$pred75)

#Model 25
mod25_b20=stan_glmer(pred25~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                   (1|Larger/ReefCluster),family = binomial(link=logit),data=alldata,control = list(adapt_delta = 0.999,max_treedepth = 15))

pp_check(mod25_b20)
stan_trace(mod25_b20)
plot(mod25_b20, 'mcmc_rhat_hist')  

#check model fit
hist(resid(mod25_b20)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(mod25_b20), x=fitted(mod25_b20)))
posteriorestimates_b20_25=as.data.frame(mod25_b20$stan_summary)
#write.csv(posteriorestimates_b20_25,"posteriorestimates_b20_25.csv")


#Model 50
mod50_b20=stan_glmer(pred50~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                   (1|Larger/ReefCluster),family = binomial(link=logit),data=alldata,control = list(adapt_delta = 0.999,max_treedepth = 15))

pp_check(mod50_b20)
stan_trace(mod50_b20)
plot(mod50_b20, 'mcmc_rhat_hist')  

#check model fit
hist(resid(mod50_b20)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(mod50_b20), x=fitted(mod50_b20)))
posteriorestimates_b20_50=as.data.frame(mod50_b20$stan_summary)
#write.csv(posteriorestimates_b20_50,"posteriorestimates_b20_50.csv")


#Model 75
mod75_b20=stan_glmer(pred75~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                   (1|Larger/ReefCluster),family = binomial(link=logit),data=alldata,control = list(adapt_delta = 0.999,max_treedepth = 15))
pp_check(mod75_b20)
stan_trace(mod75_b20)
plot(mod75_b20, 'mcmc_rhat_hist')  

#check model fit
hist(resid(mod75_b20)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(mod75_b20), x=fitted(mod75_b20)))
posteriorestimates_b20_75=as.data.frame(mod75_b20$stan_summary)
#write.csv(posteriorestimates_b20_75,"posteriorestimates_b20_75.csv")


#coefplots with same structure as manuscript
coefplot_b20_25=posteriorestimates_b20_25[2:19,]
coefplot_b20_25$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Market gravity","Nearest settlement gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_b20_25$sign=ifelse(coefplot_b20_25$`2.5%`<0 & coefplot_b20_25$`97.5%`<0, "negative",ifelse(coefplot_b20_25$`2.5%`>0 & coefplot_b20_25$`97.5%`>0, "positive", "no effect"))
coefplot_b20_25$strength=ifelse(coefplot_b20_25$sign=="no effect", "open", "closed")
coefplot_b20_25$environmental=c(1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0)
coefplot_b20_25$order=c(9,10,11,12,13,2,1,16,17,18,8,7,6,14,15,3,4,5)
coefplot_b20_25[order(coefplot_b20_25$order),]
coefplot_b20_25$variable <- factor(coefplot_b20_25$variable, levels = coefplot_b20_25$variable[order(coefplot_b20_25$order)])
coefplot_b20_25_f=
  ggplot(coefplot_b20_25,aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("25")+ theme(axis.text.y=element_blank())

coefplot_b20_50=posteriorestimates_b20_50[2:19,]
coefplot_b20_50$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Market gravity","Nearest settlement gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_b20_50$sign=ifelse(coefplot_b20_50$`2.5%`<0 & coefplot_b20_50$`97.5%`<0, "negative",ifelse(coefplot_b20_50$`2.5%`>0 & coefplot_b20_50$`97.5%`>0, "positive", "no effect"))
coefplot_b20_50$strength=ifelse(coefplot_b20_50$sign=="no effect", "open", "closed")
coefplot_b20_50$environmental=c(1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0)
coefplot_b20_50$order=c(9,10,11,12,13,2,1,16,17,18,8,7,6,14,15,3,4,5)
coefplot_b20_50[order(coefplot_b20_50$order),]
coefplot_b20_50$variable <- factor(coefplot_b20_50$variable, levels = coefplot_b20_50$variable[order(coefplot_b20_50$order)])
coefplot_b20_50_f=
  ggplot(coefplot_b20_50,aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("50")+ theme(axis.text.y=element_blank())

coefplot_b20_75=posteriorestimates_b20_75[2:19,]
coefplot_b20_75$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Market gravity","Nearest settlement gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_b20_75$sign=ifelse(coefplot_b20_75$`2.5%`<0 & coefplot_b20_75$`97.5%`<0, "negative",ifelse(coefplot_b20_75$`2.5%`>0 & coefplot_b20_75$`97.5%`>0, "positive", "no effect"))
coefplot_b20_75$strength=ifelse(coefplot_b20_75$sign=="no effect", "open", "closed")
coefplot_b20_75$environmental=c(1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0)
coefplot_b20_75$order=c(9,10,11,12,13,2,1,16,17,18,8,7,6,14,15,3,4,5)
coefplot_b20_75[order(coefplot_b20_75$order),]
coefplot_b20_75$variable <- factor(coefplot_b20_75$variable, levels = coefplot_b20_75$variable[order(coefplot_b20_75$order)])
coefplot_b20_75_f=
  ggplot(coefplot_b20_75,aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("75")+ theme(axis.text.y=element_blank())

coefplotb20_q=ggarrange(coefplot_b20_25_f,coefplot_b20_50_f,coefplot_b20_75_f, nrow=1,ncol=3)

#Create modelled data
MyData_b20<-expand.grid(sgrav_NC=seq(range(alldata$sgrav_NC)[[1]],range(alldata$sgrav_NC)[[2]],length=101),
                    Protection=c("Fished","Restricted","UnfishedHigh"),
                    DepthCategory=c(">10m"," 0-4m","4-10m"),
                    CleanHabitat=c("Crest","Flat","Lagoon_Back reef","Slope"),
                    sTotal_sampling_area=0,
                    CensusMethod=c("Distance sampling","Point intercept","Standard belt transect"),
                    sRegional_population_growth=0,
                    sOcean_prod=0,
                    sgrav_NP=0,
                    sHDI=0,
                    sReef_fish_landings_per_km2=0,
                    sClimate_stress=0,
                    sLarger_pop_size=0)

X <- model.matrix(~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                    sRegional_population_growth+sOcean_prod+sClimate_stress+
                    sHDI+sLarger_pop_size+sReef_fish_landings_per_km2,data=MyData_b20)


coefs25_b20=as.matrix(mod25_b20)
coefs25_b20=coefs25_b20[,c("(Intercept)","DepthCategory>10m","DepthCategory0-4m","CleanHabitatCrest","CleanHabitatFlat","CleanHabitatLagoon_Back reef","ProtectionRestricted","ProtectionUnfishedHigh",
                   "CensusMethodDistance sampling","CensusMethodPoint intercept","sTotal_sampling_area","sgrav_NC","sgrav_NP","sRegional_population_growth","sOcean_prod","sClimate_stress",                                                                      
                   "sHDI","sLarger_pop_size","sReef_fish_landings_per_km2" )]
fit25_b20=coefs25_b20 %*% t(X)
fit25_b20= exp(fit25_b20) / (1 + exp(fit25_b20))

MyData_b20=MyData_b20 %>% cbind(tidyMCMC(as.mcmc(fit25_b20),
                                 conf.int=T, conf.method='HPDinterval'))
setnames(MyData_b20, old = c("estimate","std.error","conf.low","conf.high"), new = c("estimate25","std.error25","conf.low25","conf.high25"))

coefs50_b20=as.matrix(mod50_b20)
coefs50_b20=coefs50_b20[,c("(Intercept)","DepthCategory>10m","DepthCategory0-4m","CleanHabitatCrest","CleanHabitatFlat","CleanHabitatLagoon_Back reef","ProtectionRestricted","ProtectionUnfishedHigh",
                   "CensusMethodDistance sampling","CensusMethodPoint intercept","sTotal_sampling_area","sgrav_NC","sgrav_NP","sRegional_population_growth","sOcean_prod","sClimate_stress",                                                                      
                   "sHDI","sLarger_pop_size","sReef_fish_landings_per_km2" )]
fit50_b20=coefs50_b20 %*% t(X)
fit50_b20= exp(fit50_b20) / (1 + exp(fit50_b20))

MyData_b20=MyData_b20 %>% cbind(tidyMCMC(as.mcmc(fit50_b20),
                                 conf.int=T, conf.method='HPDinterval'))
setnames(MyData_b20, old = c("estimate","std.error","conf.low","conf.high"), new = c("estimate50","std.error50","conf.low50","conf.high50"))

coefs75_b20=as.matrix(mod75_b20)
coefs75_b20=coefs75_b20[,c("(Intercept)","DepthCategory>10m","DepthCategory0-4m","CleanHabitatCrest","CleanHabitatFlat","CleanHabitatLagoon_Back reef","ProtectionRestricted","ProtectionUnfishedHigh",
                   "CensusMethodDistance sampling","CensusMethodPoint intercept","sTotal_sampling_area","sgrav_NC","sgrav_NP","sRegional_population_growth","sOcean_prod","sClimate_stress",                                                                      
                   "sHDI","sLarger_pop_size","sReef_fish_landings_per_km2" )]
fit75_b20=coefs75_b20 %*% t(X)
fit75_b20= exp(fit75_b20) / (1 + exp(fit75_b20))

MyData_b20=MyData_b20 %>% cbind(tidyMCMC(as.mcmc(fit75_b20),
                                 conf.int=T, conf.method='HPDinterval'))
setnames(MyData_b20, old = c("estimate","std.error","conf.low","conf.high"), new = c("estimate75","std.error75","conf.low75","conf.high75"))
#write.csv(MyData_b20,"MyData_b20_allreserves.csv",row.names=F)

# Create Delta
dat_Fished_b20 <- MyData_b20[which(MyData_b20$CleanHabitat=="Slope" & MyData_b20$Protection=="Fished" & MyData_b20$CensusMethod=="Standard belt transect" &
                             MyData_b20$DepthCategory=="4-10m"),]
dat_Reserve_b20<- MyData_b20[which(MyData_b20$CleanHabitat=="Slope" & MyData_b20$Protection=="UnfishedHigh" & MyData_b20$CensusMethod=="Standard belt transect" &
                              MyData_b20$DepthCategory=="4-10m"),]
dat_Rest_b20 <- MyData_b20[which(MyData_b20$CleanHabitat=="Slope" & MyData_b20$Protection=="Restricted" & MyData_b20$CensusMethod=="Standard belt transect" &
                           MyData_b20$DepthCategory=="4-10m"),]

D1 <- (dat_Reserve_b20[,"estimate25"] - dat_Fished_b20[,"estimate25"])
D2 <- (dat_Reserve_b20[,"estimate50"] - dat_Fished_b20[,"estimate50"])
D3 <- (dat_Reserve_b20[,"estimate75"] - dat_Fished_b20[,"estimate75"])

delta_FR_b20<- cbind(dat_Fished_b20,D1,D2,D3)

D1 <- (dat_Rest_b20[,"estimate25"] - dat_Fished_b20[,"estimate25"])
D2 <- (dat_Rest_b20[,"estimate50"] - dat_Fished_b20[,"estimate50"])
D3 <- (dat_Rest_b20[,"estimate75"] - dat_Fished_b20[,"estimate75"])

delta_FRest_b20 <- cbind(dat_Fished_b20,D1,D2,D3)

D1 <- (dat_Reserve_b20[,"estimate25"] - dat_Rest_b20[,"estimate25"])
D2 <- (dat_Reserve_b20[,"estimate50"] - dat_Rest_b20[,"estimate50"])
D3 <- (dat_Reserve_b20[,"estimate75"] - dat_Rest_b20[,"estimate75"])

delta_RestR_b20 <- cbind(dat_Fished_b20,D1,D2,D3)



# figures for manuscript

p5_bmrf <- ggplot(delta_FR_b20)+
  geom_line(aes(sgrav_NC,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_NC,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_NC,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.15,0.5),labels=NULL) + white_theme 

p6_brf <- ggplot(delta_FRest_b20)+
  geom_line(aes(sgrav_NC,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_NC,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_NC,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.justification=c(1,0), legend.position=c(0,1), line= element_blank(), axis.ticks = element_line()) +
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.15,0.5),labels=NULL) + white_theme

p7_brr <- ggplot(delta_RestR_b20)+
  geom_line(aes(sgrav_NC,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_NC,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_NC,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.justification=c(1,0), legend.position=c(0,1), line= element_blank(), axis.ticks = element_line()) +
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.15,0.5),labels=NULL) + white_theme

#plots with confidence intervals
p1_bfci <- ggplot(MyData_b20[which(MyData_b20$CleanHabitat=="Slope" & MyData_b20$Protection=="Fished" & MyData_b20$CensusMethod=="Standard belt transect" &
                                 MyData_b20$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_NC,estimate25),colour="#fbb4b9",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high25,ymin=conf.low25),fill="#fbb4b9",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate50),colour="#f768a1",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high50,ymin=conf.low50),fill="#f768a1",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high75,ymin=conf.low75),fill="#7a0177",alpha=0.3)+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + white_theme

p2_bmrci <- ggplot(MyData_b20[which(MyData_b20$CleanHabitat=="Slope" & MyData_b20$Protection=="UnfishedHigh" & MyData_b20$CensusMethod=="Standard belt transect" &
                                  MyData_b20$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_NC,estimate25),colour="#fbb4b9",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high25,ymin=conf.low25),fill="#fbb4b9",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate50),colour="#f768a1",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high50,ymin=conf.low50),fill="#f768a1",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high75,ymin=conf.low75),fill="#7a0177",alpha=0.3)+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + white_theme


p3_brci <- ggplot(MyData_b20[which(MyData_b20$CleanHabitat=="Slope" & MyData_b20$Protection=="Restricted" & MyData_b20$CensusMethod=="Standard belt transect" &
                                 MyData_b20$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_NC,estimate25),colour="#fbb4b9",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high25,ymin=conf.low25),fill="#fbb4b9",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate50),colour="#f768a1",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high50,ymin=conf.low50),fill="#f768a1",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high75,ymin=conf.low75),fill="#7a0177",alpha=0.3)+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + white_theme




# Trait diversity
alldata$func25=as.factor(alldata$func25)
alldata$func50=as.factor(alldata$func50)
alldata$func75=as.factor(alldata$func75)

#Model 25
mod25_td=stan_glmer(func25~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                   (1|Larger/ReefCluster),family = binomial(link=logit),data=alldata,control = list(adapt_delta = 0.999,max_treedepth = 15))

pp_check(mod25_td)
stan_trace(mod25_td)
plot(mod25_td, 'mcmc_rhat_hist')  

#check model fit
hist(resid(mod25_td)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(mod25_td), x=fitted(mod25_td)))
posteriorestimates_rao_25=as.data.frame(mod25_td$stan_summary)
#write.csv(posteriorestimates_rao_25,"posteriorestimates_rao_25.csv")

#model 50
mod50_td=stan_glmer(func50~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                   (1|Larger/ReefCluster),family = binomial(link=logit),data=alldata,control = list(adapt_delta = 0.999,max_treedepth = 15))

pp_check(mod50_td)
stan_trace(mod50_td)
plot(mod50_td, 'mcmc_rhat_hist')  

#check model fit
hist(resid(mod50_td)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(mod50_td), x=fitted(mod50_td)))
posteriorestimates_rao_50=as.data.frame(mod50_td$stan_summary)
#write.csv(posteriorestimates_rao_50,"posteriorestimates_rao_50.csv")

#Model 75
mod75_td=stan_glmer(func75~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                   sRegional_population_growth+sOcean_prod+sClimate_stress+
                   sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+
                   (1|Larger/ReefCluster),family = binomial(link=logit),data=alldata,control = list(adapt_delta = 0.999,max_treedepth = 15))
pp_check(mod75_td)
stan_trace(mod75_td)
plot(mod75_td, 'mcmc_rhat_hist')  

#check model fit
hist(resid(mod75_td)) 
ggplot(data=NULL)+
  geom_point(aes(y=resid(mod75_td), x=fitted(mod75_td)))
posteriorestimates_rao_75=as.data.frame(mod75_td$stan_summary)
#write.csv(posteriorestimates_rao_75,"posteriorestimates_rao_75.csv")


#coefplots with same structure as manuscript
coefplot_rao_25=posteriorestimates_rao_25[2:19,]
coefplot_rao_25$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Market gravity","Nearest settlement gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_rao_25$sign=ifelse(coefplot_rao_25$`2.5%`<0 & coefplot_rao_25$`97.5%`<0, "negative",ifelse(coefplot_rao_25$`2.5%`>0 & coefplot_rao_25$`97.5%`>0, "positive", "no effect"))
coefplot_rao_25$strength=ifelse(coefplot_rao_25$sign=="no effect", "open", "closed")
coefplot_rao_25$environmental=c(1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0)
coefplot_rao_25$order=c(9,10,11,12,13,2,1,16,17,18,8,7,6,14,15,3,4,5)
coefplot_rao_25[order(coefplot_rao_25$order),]
coefplot_rao_25$variable <- factor(coefplot_rao_25$variable, levels = coefplot_rao_25$variable[order(coefplot_rao_25$order)])
coefplot_rao_25_f=
  ggplot(coefplot_rao_25,aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("25")+ theme(axis.text.y=element_blank())

coefplot_rao_50=posteriorestimates_rao_50[2:19,]
coefplot_rao_50$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Market gravity","Nearest settlement gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_rao_50$sign=ifelse(coefplot_rao_50$`2.5%`<0 & coefplot_rao_50$`97.5%`<0, "negative",ifelse(coefplot_rao_50$`2.5%`>0 & coefplot_rao_50$`97.5%`>0, "positive", "no effect"))
coefplot_rao_50$strength=ifelse(coefplot_rao_50$sign=="no effect", "open", "closed")
coefplot_rao_50$environmental=c(1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0)
coefplot_rao_50$order=c(9,10,11,12,13,2,1,16,17,18,8,7,6,14,15,3,4,5)
coefplot_rao_50[order(coefplot_rao_50$order),]
coefplot_rao_50$variable <- factor(coefplot_rao_50$variable, levels = coefplot_rao_50$variable[order(coefplot_rao_50$order)])
coefplot_rao_50_f=
  ggplot(coefplot_rao_50,aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("50")+ theme(axis.text.y=element_blank())

coefplot_rao_75=posteriorestimates_rao_75[2:19,]
coefplot_rao_75$variable=c("> 10m Depth","< 4m Depth","Reef crest","Reef flat", "Reef lagoon", "Fishing restricted","High compliance reserve", "Distance sampling method","Point intercept method","Sampling area","Market gravity","Nearest settlement gravity","Local population growth","Ocean productivity", "Climate stress index","Human development index","Population size",  "Reef fish landings")
coefplot_rao_75$sign=ifelse(coefplot_rao_75$`2.5%`<0 & coefplot_rao_75$`97.5%`<0, "negative",ifelse(coefplot_rao_75$`2.5%`>0 & coefplot_rao_75$`97.5%`>0, "positive", "no effect"))
coefplot_rao_75$strength=ifelse(coefplot_rao_75$sign=="no effect", "open", "closed")
coefplot_rao_75$environmental=c(1,1,1,1,1,0,0,1,1,1,0,0,0,1,1,0,0,0)
coefplot_rao_75$order=c(9,10,11,12,13,2,1,16,17,18,8,7,6,14,15,3,4,5)
coefplot_rao_75[order(coefplot_rao_75$order),]
coefplot_rao_75$variable <- factor(coefplot_rao_75$variable, levels = coefplot_rao_75$variable[order(coefplot_rao_75$order)])
coefplot_rao_75_f=
  ggplot(coefplot_rao_75,aes(x=variable,y=mean,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange(aes(colour=sign,shape=strength),size=0.7)+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+coord_flip()+
  geom_hline(aes(yintercept=0),colour="dark grey",linetype="dashed")+
  ylab("")+xlab("")+ggtitle("75")+ theme(axis.text.y=element_blank())

coefplotrao_q=ggarrange(coefplot_rao_25_f,coefplot_rao_50_f,coefplot_rao_75_f, nrow=1,ncol=3)

#Create modelled data
MyData_td<-expand.grid(sgrav_NC=seq(range(alldata$sgrav_NC)[[1]],range(alldata$sgrav_NC)[[2]],length=101),
                    Protection=c("Fished","Restricted","UnfishedHigh"),
                    DepthCategory=c(">10m"," 0-4m","4-10m"),
                    CleanHabitat=c("Crest","Flat","Lagoon_Back reef","Slope"),
                    sTotal_sampling_area=0,
                    CensusMethod=c("Distance sampling","Point intercept","Standard belt transect"),
                    sRegional_population_growth=0,
                    sOcean_prod=0,
                    sgrav_NP=0,
                    sHDI=0,
                    sReef_fish_landings_per_km2=0,
                    sClimate_stress=0,
                    sLarger_pop_size=0)

X <- model.matrix(~DepthCategory+CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_NC+sgrav_NP+
                    sRegional_population_growth+sOcean_prod+sClimate_stress+
                    sHDI+sLarger_pop_size+sReef_fish_landings_per_km2,data=MyData_td)

coefs25_td=as.matrix(mod25_td)
coefs25_td=coefs25_td[,c("(Intercept)","DepthCategory>10m","DepthCategory0-4m","CleanHabitatCrest","CleanHabitatFlat","CleanHabitatLagoon_Back reef","ProtectionRestricted","ProtectionUnfishedHigh",
                   "CensusMethodDistance sampling","CensusMethodPoint intercept","sTotal_sampling_area","sgrav_NC","sgrav_NP","sRegional_population_growth","sOcean_prod","sClimate_stress",                                                                      
                   "sHDI","sLarger_pop_size","sReef_fish_landings_per_km2" )]
fit25_td=coefs25_td %*% t(X)
fit25_td= exp(fit25_td) / (1 + exp(fit25_td))

MyData_td=MyData_td %>% cbind(tidyMCMC(as.mcmc(fit25_td),
                                 conf.int=T, conf.method='HPDinterval'))
setnames(MyData_td, old = c("estimate","std.error","conf.low","conf.high"), new = c("estimate25","std.error25","conf.low25","conf.high25"))

coefs50_td=as.matrix(mod50_td)
coefs50_td=coefs50_td[,c("(Intercept)","DepthCategory>10m","DepthCategory0-4m","CleanHabitatCrest","CleanHabitatFlat","CleanHabitatLagoon_Back reef","ProtectionRestricted","ProtectionUnfishedHigh",
                   "CensusMethodDistance sampling","CensusMethodPoint intercept","sTotal_sampling_area","sgrav_NC","sgrav_NP","sRegional_population_growth","sOcean_prod","sClimate_stress",                                                                      
                   "sHDI","sLarger_pop_size","sReef_fish_landings_per_km2" )]
fit50_td=coefs50_td %*% t(X)
fit50_td= exp(fit50_td) / (1 + exp(fit50_td))

MyData_td=MyData_td %>% cbind(tidyMCMC(as.mcmc(fit50_td),
                                 conf.int=T, conf.method='HPDinterval'))
setnames(MyData_td, old = c("estimate","std.error","conf.low","conf.high"), new = c("estimate50","std.error50","conf.low50","conf.high50"))

coefs75_td=as.matrix(mod75_td)
coefs75_td=coefs75_td[,c("(Intercept)","DepthCategory>10m","DepthCategory0-4m","CleanHabitatCrest","CleanHabitatFlat","CleanHabitatLagoon_Back reef","ProtectionRestricted","ProtectionUnfishedHigh",
                   "CensusMethodDistance sampling","CensusMethodPoint intercept","sTotal_sampling_area","sgrav_NC","sgrav_NP","sRegional_population_growth","sOcean_prod","sClimate_stress",                                                                      
                   "sHDI","sLarger_pop_size","sReef_fish_landings_per_km2" )]
fit75_td=coefs75_td %*% t(X)
fit75_td= exp(fit75_td) / (1 + exp(fit75_td))

MyData_td=MyData_td %>% cbind(tidyMCMC(as.mcmc(fit75_td),
                                 conf.int=T, conf.method='HPDinterval'))
setnames(MyData_td, old = c("estimate","std.error","conf.low","conf.high"), new = c("estimate75","std.error75","conf.low75","conf.high75"))
#write.csv(MyData_td,"MyData_td_allreserves.csv",row.names=F)

# Create Delta
dat_Fished_td <- MyData_td[which(MyData_td$CleanHabitat=="Slope" & MyData_td$Protection=="Fished" & MyData_td$CensusMethod=="Standard belt transect" &
                             MyData_td$DepthCategory=="4-10m"),]
dat_Reserve_td <- MyData_td[which(MyData_td$CleanHabitat=="Slope" & MyData_td$Protection=="UnfishedHigh" & MyData_td$CensusMethod=="Standard belt transect" &
                              MyData_td$DepthCategory=="4-10m"),]
dat_Rest_td <- MyData_td[which(MyData_td$CleanHabitat=="Slope" & MyData_td$Protection=="Restricted" & MyData_td$CensusMethod=="Standard belt transect" &
                           MyData_td$DepthCategory=="4-10m"),]

D1 <- (dat_Reserve_td[,"estimate25"] - dat_Fished_td[,"estimate25"])
D2 <- (dat_Reserve_td[,"estimate50"] - dat_Fished_td[,"estimate50"])
D3 <- (dat_Reserve_td[,"estimate75"] - dat_Fished_td[,"estimate75"])

delta_FR_td<- cbind(dat_Fished_td,D1,D2,D3)

D1 <- (dat_Rest_td[,"estimate25"] - dat_Fished_td[,"estimate25"])
D2 <- (dat_Rest_td[,"estimate50"] - dat_Fished_td[,"estimate50"])
D3 <- (dat_Rest_td[,"estimate75"] - dat_Fished_td[,"estimate75"])

delta_FRest_td <- cbind(dat_Fished_td,D1,D2,D3)

D1 <- (dat_Reserve_td[,"estimate25"] - dat_Rest_td[,"estimate25"])
D2 <- (dat_Reserve_td[,"estimate50"] - dat_Rest_td[,"estimate50"])
D3 <- (dat_Reserve_td[,"estimate75"] - dat_Rest_td[,"estimate75"])

delta_RestR_td <- cbind(dat_Fished_td,D1,D2,D3)




#for manuscript
p5_rmrf <- ggplot(delta_FR_td)+
  geom_line(aes(sgrav_NC,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_NC,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_NC,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.15,0.5),labels=NULL) + white_theme 

p6_rrf <- ggplot(delta_FRest_td)+
  geom_line(aes(sgrav_NC,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_NC,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_NC,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.justification=c(1,0), legend.position=c(0,1), line= element_blank(), axis.ticks = element_line()) +
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.15,0.5),labels=NULL) + white_theme

p7_rrr <- ggplot(delta_RestR_td)+
  geom_line(aes(sgrav_NC,D1),colour="#fbb4b9",size=1)+
  geom_line(aes(sgrav_NC,D2),colour="#f768a1",size=1)+
  geom_line(aes(sgrav_NC,D3),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.justification=c(1,0), legend.position=c(0,1), line= element_blank(), axis.ticks = element_line()) +
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(-0.15,0.5),labels=NULL) + white_theme

#plots with confidence intervals
p1_rfci <- ggplot(MyData_td[which(MyData_td$CleanHabitat=="Slope" & MyData_td$Protection=="Fished" & MyData_td$CensusMethod=="Standard belt transect" &
                                 MyData_td$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_NC,estimate25),colour="#fbb4b9",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high25,ymin=conf.low25),fill="#fbb4b9",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate50),colour="#f768a1",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high50,ymin=conf.low50),fill="#f768a1",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high75,ymin=conf.low75),fill="#7a0177",alpha=0.3)+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + white_theme

p2_rmrci <- ggplot(MyData_td[which(MyData_td$CleanHabitat=="Slope" & MyData_td$Protection=="UnfishedHigh" & MyData_td$CensusMethod=="Standard belt transect" &
                                  MyData_td$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_NC,estimate25),colour="#fbb4b9",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high25,ymin=conf.low25),fill="#fbb4b9",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate50),colour="#f768a1",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high50,ymin=conf.low50),fill="#f768a1",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high75,ymin=conf.low75),fill="#7a0177",alpha=0.3)+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + white_theme


p3_rrci <- ggplot(MyData_td[which(MyData_td$CleanHabitat=="Slope" & MyData_td$Protection=="Restricted" & MyData_td$CensusMethod=="Standard belt transect" &
                                 MyData_td$DepthCategory=="4-10m"),])+
  geom_line(aes(sgrav_NC,estimate25),colour="#fbb4b9",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high25,ymin=conf.low25),fill="#fbb4b9",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate50),colour="#f768a1",size=1)+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high50,ymin=conf.low50),fill="#f768a1",alpha=0.3)+
  geom_line(aes(sgrav_NC,estimate75),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  geom_ribbon(aes(x=sgrav_NC,ymax=conf.high75,ymin=conf.low75),fill="#7a0177",alpha=0.3)+
  scale_x_continuous("",labels=NULL)+
  scale_y_continuous("",limits=c(0,1),labels=NULL) + white_theme

#manuscript
Fig.grouped=ggarrange(p1_bfci, p1_hfci,p1_rfci,p1_3fci,p2_bmrci,p2_hmrci,p2_rmrci,p2_3mrci,p3_brci,p3_hrci,p3_rrci,p3_3rci,p5_bmrf,p5_hmrf,p5_rmrf,p5_3mrf,p6_brf,p6_hrf,p6_rrf,p6_3rf, p7_brr,p7_hrr,p7_rrr,p7_3rr,nrow=6, ncol=4, widths = c(1.1,1,1,1))
annotate_figure(Fig.grouped,left = text_grob("Probability", rot = 90),
                bottom=text_grob("Std Gravity (log +1 transformed)"))


Fig.sup_quartile=ggarrange(coefplotb20_q,coefplothf_q,coefplotrao_q,coefplot3func_q, nrow=1,ncol=4,widths=c(1,1,1,1))
annotate_figure(Fig.sup_quartile,
                bottom=text_grob("Log odds ratio"))


############################################################

#plot of conditioned reserves minus all reserves:extract predicted reserve data (do it restrciting MPA and not doing it) and then substract one from the other
#Reserves_td=dat_Reserve_td[,c("sgrav_NC","estimate25","estimate50","estimate75")]
#colnames(Reserves_td)=c("sgrav_NC","estimate25_cond","estimate50_cond","estimate75_cond")
#Reserves_hf=dat_Reserve_hf[,c("sgrav_NC","estimate25","estimate50","estimate75")]
#colnames(Reserves_hf)=c("sgrav_NC","estimate25_cond","estimate50_cond","estimate75_cond")
#Reserves_b20=dat_Reserve_b20[,c("sgrav_NC","estimate25","estimate50","estimate75")]
#colnames(Reserves_b20)=c("sgrav_NC","estimate25_cond","estimate50_cond","estimate75_cond")
#Reserves_3=dat_Reserve_3[,c("sgrav_NC","estimate25","estimate50","estimate75")]
#colnames(Reserves_3)=c("sgrav_NC","estimate25_cond","estimate50_cond","estimate75_cond")
#write.csv(Reserves_hf,"Reserves_hf.csv",row.names=F)
#write.csv(Reserves_td,"Reserves_td.csv",row.names=F)
#write.csv(Reserves_b20,"Reserves_b20.csv",row.names=F)
#write.csv(Reserves_3,"Reserves_3.csv",row.names=F)

Reserves_td=dat_Reserve_td[,c("sgrav_NC","estimate25","estimate50","estimate75")]
colnames(Reserves_td)=c("sgrav_NC","estimate25_cond","estimate50_cond","estimate75_cond")
Reserves_hf=dat_Reserve_hf[,c("sgrav_NC","estimate25","estimate50","estimate75")]
colnames(Reserves_hf)=c("sgrav_NC","estimate25_cond","estimate50_cond","estimate75_cond")
Reserves_b20=dat_Reserve_b20[,c("sgrav_NC","estimate25","estimate50","estimate75")]
colnames(Reserves_b20)=c("sgrav_NC","estimate25_cond","estimate50_cond","estimate75_cond")
Reserves_3=dat_Reserve_3[,c("sgrav_NC","estimate25","estimate50","estimate75")]
colnames(Reserves_3)=c("sgrav_NC","estimate25_cond","estimate50_cond","estimate75_cond")


Reserves_hf=read.csv("Reserves_conditioned_hf.csv",head=T)
Reserves_td=read.csv("Reserves_conditioned_td.csv",head=T)
Reserves_b20=read.csv("Reserves_conditioned_b20.csv",head=T)
Reserves_3=read.csv("Reserves_conditioned_3.csv",head=T)


mpa_cond_td <- ggplot(NULL)+
  geom_line(aes(x=Reserves_conditioned_td$sgrav_NC,y=(Reserves_conditioned_td$estimate25-Reserves_td$estimate25)),colour="#fbb4b9",size=1)+
  geom_line(aes(x=Reserves_conditioned_td$sgrav_NC,y=(Reserves_conditioned_td$estimate50-Reserves_td$estimate50)),colour="#f768a1",size=1)+
  geom_line(aes(x=Reserves_conditioned_td$sgrav_NC,y=(Reserves_conditioned_td$estimate75-Reserves_td$estimate75)),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",limits=c(-0.5,2))+
  scale_y_continuous("",limits=c(-0.4,0.5),labels=NULL) + white_theme 

mpa_cond_hf <- ggplot(NULL)+
  geom_line(aes(x=Reserves_conditioned_hf$sgrav_NC,y=(Reserves_conditioned_hf$estimate25-Reserves_hf$estimate25)),colour="#fbb4b9",size=1)+
  geom_line(aes(x=Reserves_conditioned_hf$sgrav_NC,y=(Reserves_conditioned_hf$estimate50-Reserves_hf$estimate50)),colour="#f768a1",size=1)+
  geom_line(aes(x=Reserves_conditioned_hf$sgrav_NC,y=(Reserves_conditioned_hf$estimate75-Reserves_hf$estimate75)),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",limits=c(-0.5,2))+
  scale_y_continuous("",limits=c(-0.4,0.5),labels=NULL) + white_theme 

mpa_cond_b20 <- ggplot(NULL)+
  geom_line(aes(x=Reserves_conditioned_b20$sgrav_NC,y=(Reserves_conditioned_b20$estimate25-Reserves_b20$estimate25)),colour="#fbb4b9",size=1)+
  geom_line(aes(x=Reserves_conditioned_b20$sgrav_NC,y=(Reserves_conditioned_b20$estimate50-Reserves_b20$estimate50)),colour="#f768a1",size=1)+
  geom_line(aes(x=Reserves_conditioned_b20$sgrav_NC,y=(Reserves_conditioned_b20$estimate75-Reserves_b20$estimate75)),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",limits=c(-0.5,2))+
  scale_y_continuous("",limits=c(-0.4,0.5)) + white_theme 


mpa_cond_3<- ggplot(NULL)+
  geom_line(aes(x=Reserves_conditioned_3$sgrav_NC,y=(Reserves_conditioned_3$estimate25-Reserves_3$estimate25)),colour="#fbb4b9",size=1)+
  geom_line(aes(x=Reserves_conditioned_3$sgrav_NC,y=(Reserves_conditioned_3$estimate50-Reserves_3$estimate50)),colour="#f768a1",size=1)+
  geom_line(aes(x=Reserves_conditioned_3$sgrav_NC,y=(Reserves_conditioned_3$estimate75-Reserves_3$estimate75)),colour="#7a0177",size=1) + theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous("",limits=c(-0.5,2))+
  scale_y_continuous("",limits=c(-0.4,0.5),labels=NULL) + white_theme 

ggarrange(mpa_cond_b20,mpa_cond_hf,mpa_cond_td ,mpa_cond_3 ,nrow=1,ncol=4, widths=c(1.1,1,1,1))

###################################################################################

#maps of thresholds


FigureX=ggplot() + 
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group),fill = "grey", color = "grey") +
  coord_fixed(xlim = c(30, 320),  ylim = c(30, -30), ratio = 1.3)+
  geom_point(data=alldata[!is.na(alldata$tot25),], aes(x=lon2, y=Site_Lat2),colour="black", pch=16, size=3)+
  geom_point(data=alldata[alldata$tot25==3,], aes(x=lon2, y=Site_Lat2),colour="#f768a1", pch=16, size=3)+
  geom_point(data=alldata[alldata$tot50==3,], aes(x=lon2, y=Site_Lat2),colour="#c51b8a", pch=16, size=3)+
  geom_point(data=alldata[alldata$tot75==3,], aes(x=lon2, y=Site_Lat2),colour="#7a0177", pch=16, size=3)+
  
  geom_hline(yintercept =23.43695, lty=2)+
  
  geom_hline(yintercept =-23.43695, lty=2)+scale_x_continuous("Longitude",breaks=c(80,180,280),labels=c(-100,0,100))+
  scale_y_continuous("Latitude",breaks=c(-20,-10,0,10,20))+ theme_classic()


##############################################################################

#histograms of potential conservation gains: marine reserves minus openly fished
#categorize the probabilities acording to gravity and then do the same categories for all the data


#b20
delta_FR_b20$probgroup25=ifelse(delta_FR_b20$D1>0 & delta_FR_b20$D1<0.1,"<0.1", ifelse(delta_FR_b20$D1>=0.1 &delta_FR_b20$D1<0.2,"0.1-0.2",ifelse(delta_FR_b20$D1>=0.2 &delta_FR_b20$D1<0.3,"0.2-0.3",">0.3")))
mingrav25b20_cat1=delta_FR_b20$sgrav_NC[delta_FR_b20$probgroup25=="<0.1"][1]
mingrav25b20_cat2=delta_FR_b20$sgrav_NC[delta_FR_b20$probgroup25=="0.1-0.2"][1]
mingrav25b20_cat3=delta_FR_b20$sgrav_NC[delta_FR_b20$probgroup25=="0.2-0.3"][1]
mingrav25b20_cat4=delta_FR_b20$sgrav_NC[delta_FR_b20$probgroup25==">0.3"][1]
alldata$probb20_25=ifelse(alldata$sgrav_NC>=mingrav25b20_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav25b20_cat1&alldata$sgrav_NC>=mingrav25b20_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav25b20_cat2&alldata$sgrav_NC>=mingrav25b20_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav25b20_cat3,">0.3",NA))))

delta_FR_b20$probgroup50=ifelse(delta_FR_b20$D2>0 & delta_FR_b20$D2<0.1,"<0.1", ifelse(delta_FR_b20$D2>=0.1 &delta_FR_b20$D2<0.2,"0.1-0.2",ifelse(delta_FR_b20$D2>=0.2 &delta_FR_b20$D2<0.3,"0.2-0.3",">0.3")))
mingrav50b20_cat1=delta_FR_b20$sgrav_NC[delta_FR_b20$probgroup50=="<0.1"][1]
mingrav50b20_cat2=delta_FR_b20$sgrav_NC[delta_FR_b20$probgroup50=="0.1-0.2"][1]
mingrav50b20_cat3=delta_FR_b20$sgrav_NC[delta_FR_b20$probgroup50=="0.2-0.3"][1]
mingrav50b20_cat4=delta_FR_b20$sgrav_NC[delta_FR_b20$probgroup50==">0.3"][1]
alldata$probb20_50=ifelse(alldata$sgrav_NC>=mingrav50b20_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav50b20_cat1&alldata$sgrav_NC>=mingrav50b20_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav50b20_cat2&alldata$sgrav_NC>=mingrav50b20_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav50b20_cat3,">0.3",NA))))

delta_FR_b20$probgroup75=ifelse(delta_FR_b20$D3>0 & delta_FR_b20$D3<0.1,"<0.1", ifelse(delta_FR_b20$D3>=0.1 &delta_FR_b20$D3<0.2,"0.1-0.2",ifelse(delta_FR_b20$D3>=0.2 &delta_FR_b20$D3<0.3,"0.2-0.3",">0.3")))
mingrav75b20_cat1=delta_FR_b20$sgrav_NC[delta_FR_b20$probgroup75=="<0.1"][1]
mingrav75b20_cat2=delta_FR_b20$sgrav_NC[delta_FR_b20$probgroup75=="0.1-0.2"][1]
mingrav75b20_cat3=delta_FR_b20$sgrav_NC[delta_FR_b20$probgroup75=="0.2-0.3"][1]
mingrav75b20_cat4=delta_FR_b20$sgrav_NC[delta_FR_b20$probgroup75==">0.3"][1]
alldata$probb20_75=ifelse(alldata$sgrav_NC>=mingrav75b20_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav75b20_cat1&alldata$sgrav_NC>=mingrav75b20_cat2,"0.1-0.2",NA))


b20a=ggplot(alldata, aes(x=probb20_25,y=..count../sum(..count..))) +
  geom_bar(fill="#fbb4b9")+theme_classic()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1)) +
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=NULL)+xlab("")+ ylab("25 target")

b20b=ggplot(alldata, aes(x=probb20_50,y=..count../sum(..count..))) +
  geom_bar(fill="#f768a1")+theme_classic()+ 
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1)) +
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=NULL)+xlab("")+ ylab("50 target")

b20c=ggplot(alldata, aes(x=probb20_75,y=..count../sum(..count..))) +
  geom_bar(fill="#7a0177")+theme_classic()+ 
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1)) +
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4"))+xlab("")+ ylab("75 target")
ggarrange(b20a,b20b,b20c,nrow=3, ncol=1, heights=c(1,1,1.1))

#3
delta_FR_3$probgroup25=ifelse(delta_FR_3$D1>0 & delta_FR_3$D1<0.1,"<0.1", ifelse(delta_FR_3$D1>=0.1 &delta_FR_3$D1<0.2,"0.1-0.2",ifelse(delta_FR_3$D1>=0.2 &delta_FR_3$D1<0.3,"0.2-0.3",">0.3")))
mingrav253_cat1=delta_FR_3$sgrav_NC[delta_FR_3$probgroup25=="<0.1"][1]
mingrav253_cat2=delta_FR_3$sgrav_NC[delta_FR_3$probgroup25=="0.1-0.2"][1]
mingrav253_cat3=delta_FR_3$sgrav_NC[delta_FR_3$probgroup25=="0.2-0.3"][1]
mingrav253_cat4=delta_FR_3$sgrav_NC[delta_FR_3$probgroup25==">0.3"][1]
alldata$prob3_25=ifelse(alldata$sgrav_NC>=mingrav253_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav253_cat1&alldata$sgrav_NC>=mingrav253_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav253_cat2&alldata$sgrav_NC>=mingrav253_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav253_cat3,">0.3",NA))))

delta_FR_3$probgroup50=ifelse(delta_FR_3$D2>0 &delta_FR_3$D2<0.1,"<0.1", ifelse(delta_FR_3$D2>=0.1 &delta_FR_3$D2<0.2,"0.1-0.2",ifelse(delta_FR_3$D2>=0.2 &delta_FR_3$D2<0.3,"0.2-0.3",">0.3")))
mingrav503_cat1=delta_FR_3$sgrav_NC[delta_FR_3$probgroup50=="<0.1"][1]
mingrav503_cat2=delta_FR_3$sgrav_NC[delta_FR_3$probgroup50=="0.1-0.2"][1]
mingrav503_cat3=delta_FR_3$sgrav_NC[delta_FR_3$probgroup50=="0.2-0.3"][1]
mingrav503_cat4=delta_FR_3$sgrav_NC[delta_FR_3$probgroup50==">0.3"][1]
alldata$prob3_50=ifelse(alldata$sgrav_NC>=mingrav503_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav503_cat1&alldata$sgrav_NC>=mingrav503_cat2,"0.1-0.2",NA))

delta_FR_3$probgroup75=ifelse(delta_FR_3$D3>0 &delta_FR_3$D3<0.1,"<0.1", ifelse(delta_FR_3$D3>=0.1 &delta_FR_3$D3<0.2,"0.1-0.2",ifelse(delta_FR_3$D3>=0.2 &delta_FR_3$D3<0.3,"0.2-0.3",">0.3")))
mingrav753_cat1=delta_FR_3$sgrav_NC[delta_FR_3$probgroup75=="<0.1"][1]
mingrav753_cat2=delta_FR_3$sgrav_NC[delta_FR_3$probgroup75=="0.1-0.2"][1]
mingrav753_cat3=delta_FR_3$sgrav_NC[delta_FR_3$probgroup75=="0.2-0.3"][1]
mingrav753_cat4=delta_FR_3$sgrav_NC[delta_FR_3$probgroup75==">0.3"][1]
alldata$prob3_75=ifelse(alldata$sgrav_NC>=mingrav753_cat1,"<0.1","<0.1")


g3a=ggplot(alldata, aes(x=prob3_25,y=..count../sum(..count..))) +
  geom_bar(fill="#fbb4b9")+theme_classic()+scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=NULL)+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) +
  xlab("")+ ylab("")
g3b=ggplot(alldata, aes(x=prob3_50,y=..count../sum(..count..))) +
  geom_bar(fill="#f768a1")+theme_classic()+ 
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=NULL)+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) +
  xlab("")+ ylab("")
g3c=ggplot(alldata, aes(x=prob3_75,y=..count../sum(..count..))) +
  geom_bar(fill="#7a0177")+theme_classic()+ 
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4"))+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) +
  xlab("")+ ylab("")
ggarrange(g3a,g3b,g3c,nrow=3, ncol=1, heights=c(1,1,1.1))


#td
delta_FR_td$probgroup25=ifelse(delta_FR_td$D1<0.1,"<0.1", ifelse(delta_FR_td$D1>=0.1 &delta_FR_td$D1<0.2,"0.1-0.2",ifelse(delta_FR_td$D1>=0.2 &delta_FR_td$D1<0.3,"0.2-0.3",">0.3")))
mingrav25td_cat1=delta_FR_td$sgrav_NC[delta_FR_td$probgroup25=="<0.1"][1]
mingrav25td_cat2=delta_FR_td$sgrav_NC[delta_FR_td$probgroup25=="0.1-0.2"][1]
mingrav25td_cat3=delta_FR_td$sgrav_NC[delta_FR_td$probgroup25=="0.2-0.3"][1]
mingrav25td_cat4=delta_FR_td$sgrav_NC[delta_FR_td$probgroup25==">0.3"][1]
alldata$probtd_25=ifelse(alldata$sgrav_NC>=mingrav25td_cat1,"<0.1","<0.1")

delta_FR_td$probgroup50=ifelse(delta_FR_td$D2>0 &delta_FR_td$D2<0.1,"<0.1", ifelse(delta_FR_td$D2>=0.1 &delta_FR_td$D2<0.2,"0.1-0.2",ifelse(delta_FR_td$D2>=0.2 &delta_FR_td$D2<0.3,"0.2-0.3",">0.3")))
mingrav50td_cat1=delta_FR_td$sgrav_NC[delta_FR_td$probgroup50=="<0.1"][1]
mingrav50td_cat2=delta_FR_td$sgrav_NC[delta_FR_td$probgroup50=="0.1-0.2"][1]
mingrav50td_cat3=delta_FR_td$sgrav_NC[delta_FR_td$probgroup50=="0.2-0.3"][1]
mingrav50td_cat4=delta_FR_td$sgrav_NC[delta_FR_td$probgroup50==">0.3"][1]
alldata$probtd_50=ifelse(alldata$sgrav_NC>=mingrav50td_cat2,"0.1-0.2","<0.1")

delta_FR_td$probgroup75=ifelse(delta_FR_td$D3>0 &delta_FR_td$D3<0.1,"<0.1", ifelse(delta_FR_td$D3>=0.1 &delta_FR_td$D3<0.2,"0.1-0.2",ifelse(delta_FR_td$D3>=0.2 &delta_FR_td$D3<0.3,"0.2-0.3",">0.3")))
mingrav75td_cat1=delta_FR_td$sgrav_NC[delta_FR_td$probgroup75=="<0.1"][1]
mingrav75td_cat2=delta_FR_td$sgrav_NC[delta_FR_td$probgroup75=="0.1-0.2"][1]
mingrav75td_cat3=delta_FR_td$sgrav_NC[delta_FR_td$probgroup75=="0.2-0.3"][1]
mingrav75td_cat4=delta_FR_td$sgrav_NC[delta_FR_td$probgroup75==">0.3"][1]
alldata$probtd_75=ifelse(alldata$sgrav_NC>=mingrav75td_cat1,"<0.1","<0.1")


tda=ggplot(alldata, aes(x=probtd_25,y=..count../sum(..count..))) +
  geom_bar(fill="#fbb4b9")+theme_classic()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) +
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=NULL)+xlab("")+ ylab("")
tdb=ggplot(alldata, aes(x=probtd_50,y=..count../sum(..count..))) +
  geom_bar(fill="#f768a1")+theme_classic()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) + 
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=NULL)+xlab("")+ ylab("")
tdc=ggplot(alldata, aes(x=probtd_75,y=..count../sum(..count..))) +
  geom_bar(fill="#7a0177")+theme_classic()+ 
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) +
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4"))+xlab("")+ ylab("")
ggarrange(tda,tdb,tdc,nrow=3, ncol=1, heights=c(1,1,1.1))

#hf
delta_FR_hf$probgroup25=ifelse(delta_FR_hf$D1>0 &delta_FR_hf$D1<0.1,"<0.1", ifelse(delta_FR_hf$D1>=0.1 &delta_FR_hf$D1<0.2,"0.1-0.2",ifelse(delta_FR_hf$D1>=0.2 &delta_FR_hf$D1<0.3,"0.2-0.3",">0.3")))
mingrav25hf_cat1=delta_FR_hf$sgrav_NC[delta_FR_hf$probgroup25=="<0.1"][1]
mingrav25hf_cat2=delta_FR_hf$sgrav_NC[delta_FR_hf$probgroup25=="0.1-0.2"][1]
mingrav25hf_cat3=delta_FR_hf$sgrav_NC[delta_FR_hf$probgroup25=="0.2-0.3"][1]
mingrav25hf_cat4=delta_FR_hf$sgrav_NC[delta_FR_hf$probgroup25==">0.3"][1]
alldata$probhf_25=ifelse(alldata$sgrav_NC>=mingrav25hf_cat3,"0.2-0.3","0.1-0.2")


delta_FR_hf$probgroup50=ifelse(delta_FR_hf$D2>0 &delta_FR_hf$D2<0.1,"<0.1", ifelse(delta_FR_hf$D2>=0.1 &delta_FR_hf$D2<0.2,"0.1-0.2",ifelse(delta_FR_hf$D2>=0.2 &delta_FR_hf$D2<0.3,"0.2-0.3",">0.3")))
mingrav50hf_cat1=delta_FR_hf$sgrav_NC[delta_FR_hf$probgroup50=="<0.1"][1]
mingrav50hf_cat2=delta_FR_hf$sgrav_NC[delta_FR_hf$probgroup50=="0.1-0.2"][1]
mingrav50hf_cat3=delta_FR_hf$sgrav_NC[delta_FR_hf$probgroup50=="0.2-0.3"][1]
mingrav50hf_cat4=delta_FR_hf$sgrav_NC[delta_FR_hf$probgroup50==">0.3"][1]
alldata$probhf_50=ifelse(alldata$sgrav_NC>=mingrav50hf_cat4,">0.3",ifelse(alldata$sgrav_NC<mingrav50hf_cat4 &alldata$sgrav_NC>=mingrav50hf_cat3,"0.2-0.3","0.1-0.2"))


delta_FR_hf$probgroup75=ifelse(delta_FR_hf$D3>0 &delta_FR_hf$D3<0.1,"<0.1", ifelse(delta_FR_hf$D3>=0.1 &delta_FR_hf$D3<0.2,"0.1-0.2",ifelse(delta_FR_hf$D3>=0.2 &delta_FR_hf$D3<0.3,"0.2-0.3",">0.3")))
mingrav75hf_cat1=delta_FR_hf$sgrav_NC[delta_FR_hf$probgroup75=="<0.1"][1]
mingrav75hf_cat2=delta_FR_hf$sgrav_NC[delta_FR_hf$probgroup75=="0.1-0.2"][1]
mingrav75hf_cat3=delta_FR_hf$sgrav_NC[delta_FR_hf$probgroup75=="0.2-0.3"][1]
mingrav75hf_cat4=delta_FR_hf$sgrav_NC[delta_FR_hf$probgroup75==">0.3"][1]
alldata$probhf_75=ifelse(alldata$sgrav_NC>=mingrav75hf_cat3,"0.2-0.3",">0.3")


hfa=ggplot(alldata, aes(x=probhf_25,y=..count../sum(..count..))) +
  geom_bar(fill="#fbb4b9")+theme_classic()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) +
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=NULL)+xlab("")+ ylab("")
hfb=ggplot(alldata, aes(x=probhf_50,y=..count../sum(..count..))) +
  geom_bar(fill="#f768a1")+theme_classic()+ 
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) +
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=NULL)+xlab("")+ ylab("")
hfc=ggplot(alldata, aes(x=probhf_75,y=..count../sum(..count..))) +
  geom_bar(fill="#7a0177")+theme_classic()+ 
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) +
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4"))+xlab("")+ ylab("")
ggarrange(hfa,hfb,hfc,nrow=3, ncol=1, heights=c(1,1,1.1))


gains=ggarrange(b20a,hfa,tda,g3a,b20b,hfb,tdb,g3b,b20c,hfc,tdc,g3c, labels=c("A","B","C","D","E","F","G","H","I","J","K","L"),nrow=3,ncol=4, widths=c(1.2,1,1,1),font.label = list(size = 10, color = "black", face = "bold"))
annotate_figure(gains,
                bottom=text_grob("Difference in probability between openly fished reefs and marine reserves"))

#histograms of potential conservation gains: restricted reefs minus openly fished

#b20
delta_FRest_b20$probgroup25=ifelse( delta_FRest_b20$D1<0.1,"<0.1", ifelse(delta_FRest_b20$D1>=0.1 &delta_FRest_b20$D1<0.2,"0.1-0.2",ifelse(delta_FRest_b20$D1>=0.2 &delta_FRest_b20$D1<0.3,"0.2-0.3",">0.3")))
mingrav25b20_cat1=delta_FRest_b20$sgrav_NC[delta_FRest_b20$probgroup25=="<0.1"][1]
mingrav25b20_cat2=delta_FRest_b20$sgrav_NC[delta_FRest_b20$probgroup25=="0.1-0.2"][1]
mingrav25b20_cat3=delta_FRest_b20$sgrav_NC[delta_FRest_b20$probgroup25=="0.2-0.3"][1]
mingrav25b20_cat4=delta_FRest_b20$sgrav_NC[delta_FRest_b20$probgroup25==">0.3"][1]
alldata$probb20_25=ifelse(alldata$sgrav_NC>=mingrav25b20_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav25b20_cat1&alldata$sgrav_NC>=mingrav25b20_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav25b20_cat2,"0.2-0.3",">0.3")))

delta_FRest_b20$probgroup50=ifelse( delta_FRest_b20$D2<0.1,"<0.1", ifelse(delta_FRest_b20$D2>=0.1 &delta_FRest_b20$D2<0.2,"0.1-0.2",ifelse(delta_FRest_b20$D2>=0.2 &delta_FRest_b20$D2<0.3,"0.2-0.3",">0.3")))
mingrav50b20_cat1=delta_FRest_b20$sgrav_NC[delta_FRest_b20$probgroup50=="<0.1"][1]
mingrav50b20_cat2=delta_FRest_b20$sgrav_NC[delta_FRest_b20$probgroup50=="0.1-0.2"][1]
mingrav50b20_cat3=delta_FRest_b20$sgrav_NC[delta_FRest_b20$probgroup50=="0.2-0.3"][1]
mingrav50b20_cat4=delta_FRest_b20$sgrav_NC[delta_FRest_b20$probgroup50==">0.3"][1]
alldata$probb20_50=ifelse(alldata$sgrav_NC>=mingrav50b20_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav50b20_cat1&alldata$sgrav_NC>=mingrav50b20_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav50b20_cat2&alldata$sgrav_NC>=mingrav50b20_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav50b20_cat3,">0.3",NA))))


delta_FRest_b20$probgroup75=ifelse( delta_FRest_b20$D3<0.1,"<0.1", ifelse(delta_FRest_b20$D3>=0.1 &delta_FRest_b20$D3<0.2,"0.1-0.2",ifelse(delta_FRest_b20$D3>=0.2 &delta_FRest_b20$D3<0.3,"0.2-0.3",">0.3")))
mingrav75b20_cat1=delta_FRest_b20$sgrav_NC[delta_FRest_b20$probgroup75=="<0.1"][1]
mingrav75b20_cat2=delta_FRest_b20$sgrav_NC[delta_FRest_b20$probgroup75=="0.1-0.2"][1]
mingrav75b20_cat3=delta_FRest_b20$sgrav_NC[delta_FRest_b20$probgroup75=="0.2-0.3"][1]
mingrav75b20_cat4=delta_FRest_b20$sgrav_NC[delta_FRest_b20$probgroup75==">0.3"][1]
alldata$probb20_75=ifelse(alldata$sgrav_NC>=mingrav75b20_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav75b20_cat1&alldata$sgrav_NC>=mingrav75b20_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav75b20_cat2&alldata$sgrav_NC>=mingrav75b20_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav75b20_cat3,">0.3",NA))))



b20a=ggplot(alldata, aes(x=probb20_25,y=..count../sum(..count..))) +
  geom_bar(fill="#fbb4b9")+theme_classic()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1)) +
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=NULL)+xlab("")+ ylab("25 threshold")

b20b=ggplot(alldata, aes(x=probb20_50,y=..count../sum(..count..))) +
  geom_bar(fill="#f768a1")+theme_classic()+ 
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1)) +
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=NULL)+xlab("")+ ylab("50 threshold")

b20c=ggplot(alldata, aes(x=probb20_75,y=..count../sum(..count..))) +
  geom_bar(fill="#7a0177")+theme_classic()+ 
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1)) +
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4"))+xlab("")+ ylab("75 threshold")
ggarrange(b20a,b20b,b20c,nrow=3, ncol=1, heights=c(1,1,1.1))

#3
delta_FRest_3$probgroup25=ifelse( delta_FRest_3$D1<0.1,"<0.1", ifelse(delta_FRest_3$D1>=0.1 &delta_FRest_3$D1<0.2,"0.1-0.2",ifelse(delta_FRest_3$D1>=0.2 &delta_FRest_3$D1<0.3,"0.2-0.3",">0.3")))
mingrav253_cat1=delta_FRest_3$sgrav_NC[delta_FRest_3$probgroup25=="<0.1"][1]
mingrav253_cat2=delta_FRest_3$sgrav_NC[delta_FRest_3$probgroup25=="0.1-0.2"][1]
mingrav253_cat3=delta_FRest_3$sgrav_NC[delta_FRest_3$probgroup25=="0.2-0.3"][1]
mingrav253_cat4=delta_FRest_3$sgrav_NC[delta_FRest_3$probgroup25==">0.3"][1]
alldata$prob3_25=ifelse(alldata$sgrav_NC>=mingrav253_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav253_cat1&alldata$sgrav_NC>=mingrav253_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav253_cat2,"0.2-0.3",">0.3")))

delta_FRest_3$probgroup50=ifelse(delta_FRest_3$D2<0.1,"<0.1", ifelse(delta_FRest_3$D2>=0.1 &delta_FRest_3$D2<0.2,"0.1-0.2",ifelse(delta_FRest_3$D2>=0.2 &delta_FRest_3$D2<0.3,"0.2-0.3",">0.3")))
mingrav503_cat1=delta_FRest_3$sgrav_NC[delta_FRest_3$probgroup50=="<0.1"][1]
mingrav503_cat2=delta_FRest_3$sgrav_NC[delta_FRest_3$probgroup50=="0.1-0.2"][1]
mingrav503_cat3=delta_FRest_3$sgrav_NC[delta_FRest_3$probgroup50=="0.2-0.3"][1]
mingrav503_cat4=delta_FRest_3$sgrav_NC[delta_FRest_3$probgroup50==">0.3"][1]
alldata$prob3_50=ifelse(alldata$sgrav_NC>=mingrav503_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav503_cat1&alldata$sgrav_NC>=mingrav503_cat2,"0.1-0.2","0.2-0.3"))

delta_FRest_3$probgroup75=ifelse(delta_FRest_3$D3<0.1,"<0.1", ifelse(delta_FRest_3$D3>=0.1 &delta_FRest_3$D3<0.2,"0.1-0.2",ifelse(delta_FRest_3$D3>=0.2 &delta_FRest_3$D3<0.3,"0.2-0.3",">0.3")))
mingrav753_cat1=delta_FRest_3$sgrav_NC[delta_FRest_3$probgroup75=="<0.1"][1]
mingrav753_cat2=delta_FRest_3$sgrav_NC[delta_FRest_3$probgroup75=="0.1-0.2"][1]
mingrav753_cat3=delta_FRest_3$sgrav_NC[delta_FRest_3$probgroup75=="0.2-0.3"][1]
mingrav753_cat4=delta_FRest_3$sgrav_NC[delta_FRest_3$probgroup75==">0.3"][1]
alldata$prob3_75=ifelse(alldata$sgrav_NC>=mingrav753_cat1,"<0.1","<0.1")

g3a=ggplot(alldata, aes(x=prob3_25,y=..count../sum(..count..))) +
  geom_bar(fill="#fbb4b9")+theme_classic()+scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=NULL)+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) +
  xlab("")+ ylab("")
g3b=ggplot(alldata, aes(x=prob3_50,y=..count../sum(..count..))) +
  geom_bar(fill="#f768a1")+theme_classic()+ 
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=NULL)+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) +
  xlab("")+ ylab("")
g3c=ggplot(alldata, aes(x=prob3_75,y=..count../sum(..count..))) +
  geom_bar(fill="#7a0177")+theme_classic()+ 
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4"))+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) +
  xlab("")+ ylab("")
ggarrange(g3a,g3b,g3c,nrow=3, ncol=1, heights=c(1,1,1.1))


#td
delta_FRest_td$probgroup25=ifelse(delta_FRest_td$D1<0.1,"<0.1", ifelse(delta_FRest_td$D1>=0.1 &delta_FRest_td$D1<0.2,"0.1-0.2",ifelse(delta_FRest_td$D1>=0.2 &delta_FRest_td$D1<0.3,"0.2-0.3",">0.3")))
mingrav25td_cat1=delta_FRest_td$sgrav_NC[delta_FRest_td$probgroup25=="<0.1"][1]
mingrav25td_cat2=delta_FRest_td$sgrav_NC[delta_FRest_td$probgroup25=="0.1-0.2"][1]
mingrav25td_cat3=delta_FRest_td$sgrav_NC[delta_FRest_td$probgroup25=="0.2-0.3"][1]
mingrav25td_cat4=delta_FRest_td$sgrav_NC[delta_FRest_td$probgroup25==">0.3"][1]
alldata$probtd_25=ifelse(alldata$sgrav_NC>=mingrav25td_cat0,"<=0",ifelse(alldata$sgrav_NC<mingrav25td_cat0&alldata$sgrav_NC>=mingrav25td_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav25td_cat1&alldata$sgrav_NC>=mingrav25td_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav25td_cat2&alldata$sgrav_NC>=mingrav25td_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav25td_cat3,">0.3",NA)))))
alldata$probtd_25=ifelse(alldata$sgrav_NC>=mingrav25td_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav25td_cat1&alldata$sgrav_NC>=mingrav25td_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav25td_cat2&alldata$sgrav_NC>=mingrav25td_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav25td_cat3,">0.3",NA))))
alldata$probtd_25=ifelse(alldata$sgrav_NC>=mingrav25td_cat1,"<0.1","<0.1")


delta_FRest_td$probgroup50=ifelse(delta_FRest_td$D2<0.1,"<0.1", ifelse(delta_FRest_td$D2>=0.1 &delta_FRest_td$D2<0.2,"0.1-0.2",ifelse(delta_FRest_td$D2>=0.2 &delta_FRest_td$D2<0.3,"0.2-0.3",">0.3")))
mingrav50td_cat1=delta_FRest_td$sgrav_NC[delta_FRest_td$probgroup50=="<0.1"][1]
mingrav50td_cat2=delta_FRest_td$sgrav_NC[delta_FRest_td$probgroup50=="0.1-0.2"][1]
mingrav50td_cat3=delta_FRest_td$sgrav_NC[delta_FRest_td$probgroup50=="0.2-0.3"][1]
mingrav50td_cat4=delta_FRest_td$sgrav_NC[delta_FRest_td$probgroup50==">0.3"][1]
alldata$probtd_50=ifelse(alldata$sgrav_NC>=mingrav50td_cat0,"<=0",ifelse(alldata$sgrav_NC<mingrav50td_cat0&alldata$sgrav_NC>=mingrav50td_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav50td_cat1&alldata$sgrav_NC>=mingrav50td_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav50td_cat2&alldata$sgrav_NC>=mingrav50td_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav50td_cat3,">0.3",NA)))))
alldata$probtd_50=ifelse(alldata$sgrav_NC>=mingrav50td_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav50td_cat1&alldata$sgrav_NC>=mingrav50td_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav50td_cat2&alldata$sgrav_NC>=mingrav50td_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav50td_cat3,">0.3",NA))))
alldata$probtd_50=ifelse(alldata$sgrav_NC>=mingrav50td_cat1,"<0.1","<0.1")


delta_FRest_td$probgroup75=ifelse(delta_FRest_td$D3<0.1,"<0.1", ifelse(delta_FRest_td$D3>=0.1 &delta_FRest_td$D3<0.2,"0.1-0.2",ifelse(delta_FRest_td$D3>=0.2 &delta_FRest_td$D3<0.3,"0.2-0.3",">0.3")))
mingrav75td_cat1=delta_FRest_td$sgrav_NC[delta_FRest_td$probgroup75=="<0.1"][1]
mingrav75td_cat2=delta_FRest_td$sgrav_NC[delta_FRest_td$probgroup75=="0.1-0.2"][1]
mingrav75td_cat3=delta_FRest_td$sgrav_NC[delta_FRest_td$probgroup75=="0.2-0.3"][1]
mingrav75td_cat4=delta_FRest_td$sgrav_NC[delta_FRest_td$probgroup75==">0.3"][1]
alldata$probtd_75=ifelse(alldata$sgrav_NC>=mingrav75td_cat0,"<=0",ifelse(alldata$sgrav_NC<mingrav75td_cat0&alldata$sgrav_NC>=mingrav75td_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav75td_cat1&alldata$sgrav_NC>=mingrav75td_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav75td_cat2&alldata$sgrav_NC>=mingrav75td_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav75td_cat3,">0.3",NA)))))
alldata$probtd_75=ifelse(alldata$sgrav_NC>=mingrav75td_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav75td_cat1&alldata$sgrav_NC>=mingrav75td_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav75td_cat2&alldata$sgrav_NC>=mingrav75td_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav75td_cat3,">0.3",NA))))
alldata$probtd_75=ifelse(alldata$sgrav_NC>=mingrav75td_cat1,"<0.1","<0.1")


tda=ggplot(alldata, aes(x=probtd_25,y=..count../sum(..count..))) +
  geom_bar(fill="#fbb4b9")+theme_classic()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) +
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=NULL)+xlab("")+ ylab("")
tdb=ggplot(alldata, aes(x=probtd_50,y=..count../sum(..count..))) +
  geom_bar(fill="#f768a1")+theme_classic()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) + 
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=NULL)+xlab("")+ ylab("")
tdc=ggplot(alldata, aes(x=probtd_75,y=..count../sum(..count..))) +
  geom_bar(fill="#7a0177")+theme_classic()+ 
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) +
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4"))+xlab("")+ ylab("")
ggarrange(tda,tdb,tdc,nrow=3, ncol=1, heights=c(1,1,1.1))

#hf
delta_FRest_hf$probgroup25=ifelse(delta_FRest_hf$D1<0.1,"<0.1", ifelse(delta_FRest_hf$D1>=0.1 &delta_FRest_hf$D1<0.2,"0.1-0.2",ifelse(delta_FRest_hf$D1>=0.2 &delta_FRest_hf$D1<0.3,"0.2-0.3",">0.3")))
mingrav25hf_cat1=delta_FRest_hf$sgrav_NC[delta_FRest_hf$probgroup25=="<0.1"][1]
mingrav25hf_cat2=delta_FRest_hf$sgrav_NC[delta_FRest_hf$probgroup25=="0.1-0.2"][1]
mingrav25hf_cat3=delta_FRest_hf$sgrav_NC[delta_FRest_hf$probgroup25=="0.2-0.3"][1]
mingrav25hf_cat4=delta_FRest_hf$sgrav_NC[delta_FRest_hf$probgroup25==">0.3"][1]
alldata$probhf_25=ifelse(alldata$sgrav_NC>=mingrav25hf_cat4,">0.3",ifelse(alldata$sgrav_NC<mingrav25hf_cat4&alldata$sgrav_NC>=mingrav25hf_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav25hf_cat3&alldata$sgrav_NC>=mingrav25hf_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav25hf_cat2&alldata$sgrav_NC>=mingrav25hf_cat1,"<0.1",,ifelse(alldata$sgrav_NC<mingrav25hf_cat1&alldata$sgrav_NC>=mingrav25hf_cat0,"<=0",NA)))))
alldata$probhf_25=ifelse(alldata$sgrav_NC>=mingrav25hf_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav25hf_cat3&alldata$sgrav_NC>=mingrav25hf_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav25hf_cat2,"<0.1",NA)))
alldata$probhf_25=ifelse(alldata$sgrav_NC>=mingrav25hf_cat3,"0.2-0.3","0.1-0.2")
alldata$probhf_25=ifelse(alldata$sgrav_NC>=mingrav25hf_cat4,">0.3",ifelse(alldata$sgrav_NC<mingrav25hf_cat4 &alldata$sgrav_NC>=mingrav25hf_cat3,"0.2-0.3","0.1-0.2"))

#problem with R and logical statements(https://www.reddit.com/r/rstats/comments/34p0bm/logical_operators_not_working_as_expected_r_says/)
alldata$sgrav_NC[alldata$sgrav_NC<mingrav25hf_cat2]

delta_FRest_hf$probgroup50=ifelse(delta_FRest_hf$D2<0.1,"<0.1", ifelse(delta_FRest_hf$D2>=0.1 &delta_FRest_hf$D2<0.2,"0.1-0.2",ifelse(delta_FRest_hf$D2>=0.2 &delta_FRest_hf$D2<0.3,"0.2-0.3",">0.3")))
mingrav50hf_cat1=delta_FRest_hf$sgrav_NC[delta_FRest_hf$probgroup50=="<0.1"][1]
mingrav50hf_cat2=delta_FRest_hf$sgrav_NC[delta_FRest_hf$probgroup50=="0.1-0.2"][1]
mingrav50hf_cat3=delta_FRest_hf$sgrav_NC[delta_FRest_hf$probgroup50=="0.2-0.3"][1]
mingrav50hf_cat4=delta_FRest_hf$sgrav_NC[delta_FRest_hf$probgroup50==">0.3"][1]
alldata$probhf_50=ifelse(alldata$sgrav_NC>=mingrav50hf_cat4,">0.3",ifelse(alldata$sgrav_NC<mingrav50hf_cat4&alldata$sgrav_NC>=mingrav50hf_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav50hf_cat3&alldata$sgrav_NC>=mingrav50hf_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav50hf_cat2&alldata$sgrav_NC>=mingrav50hf_cat1,"<0.1",,ifelse(alldata$sgrav_NC<mingrav50hf_cat1&alldata$sgrav_NC>=mingrav50hf_cat0,"<=0",NA)))))
alldata$probhf_50=ifelse(alldata$sgrav_NC>=mingrav50hf_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav50hf_cat3&alldata$sgrav_NC>=mingrav50hf_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav50hf_cat2,"<0.1",NA)))
alldata$probhf_50=ifelse(alldata$sgrav_NC>=mingrav50hf_cat3,"0.2-0.3","0.1-0.2")
alldata$probhf_50=ifelse(alldata$sgrav_NC>=mingrav50hf_cat4,">0.3",ifelse(alldata$sgrav_NC<mingrav50hf_cat4&alldata$sgrav_NC>=mingrav50hf_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav50hf_cat3&alldata$sgrav_NC>=mingrav50hf_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav50hf_cat2&alldata$sgrav_NC>=mingrav50hf_cat1,"<0.1",NA))))

delta_FRest_hf$probgroup75=ifelse(delta_FRest_hf$D3<0.1,"<0.1", ifelse(delta_FRest_hf$D3>=0.1 &delta_FRest_hf$D3<0.2,"0.1-0.2",ifelse(delta_FRest_hf$D3>=0.2 &delta_FRest_hf$D3<0.3,"0.2-0.3",">0.3")))
mingrav75hf_cat1=delta_FRest_hf$sgrav_NC[delta_FRest_hf$probgroup75=="<0.1"][1]
mingrav75hf_cat2=delta_FRest_hf$sgrav_NC[delta_FRest_hf$probgroup75=="0.1-0.2"][1]
mingrav75hf_cat3=delta_FRest_hf$sgrav_NC[delta_FRest_hf$probgroup75=="0.2-0.3"][1]
mingrav75hf_cat4=delta_FRest_hf$sgrav_NC[delta_FRest_hf$probgroup75==">0.3"][1]
alldata$probhf_75=ifelse(alldata$sgrav_NC>=mingrav75hf_cat4,">0.3",ifelse(alldata$sgrav_NC<mingrav75hf_cat4&alldata$sgrav_NC>=mingrav75hf_cat3,"0.2-0.3",ifelse(alldata$sgrav_NC<mingrav75hf_cat3&alldata$sgrav_NC>=mingrav75hf_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav75hf_cat2&alldata$sgrav_NC>=mingrav75hf_cat1,"<0.1",,ifelse(alldata$sgrav_NC<mingrav75hf_cat1&alldata$sgrav_NC>=mingrav75hf_cat0,"<=0",NA)))))
alldata$probhf_75=ifelse(alldata$sgrav_NC>=mingrav75hf_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav75hf_cat1 & alldata$sgrav_NC>=mingrav75hf_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav75hf_cat2,"0.2-0.3",NA)))
alldata$probhf_75=ifelse(alldata$sgrav_NC>=mingrav75hf_cat1,"<0.1",ifelse(alldata$sgrav_NC<mingrav75hf_cat1 & alldata$sgrav_NC>=mingrav75hf_cat2,"0.1-0.2",ifelse(alldata$sgrav_NC<mingrav75hf_cat2,"0.2-0.3",NA)))
alldata$probhf_75=ifelse(alldata$sgrav_NC>=mingrav75hf_cat2,"0.1-0.2","0.2-0.3")


hfa=ggplot(alldata, aes(x=probhf_25,y=..count../sum(..count..))) +
  geom_bar(fill="#fbb4b9")+theme_classic()+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) +
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=NULL)+xlab("")+ ylab("")
hfb=ggplot(alldata, aes(x=probhf_50,y=..count../sum(..count..))) +
  geom_bar(fill="#f768a1")+theme_classic()+ 
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) +
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=NULL)+xlab("")+ ylab("")
hfc=ggplot(alldata, aes(x=probhf_75,y=..count../sum(..count..))) +
  geom_bar(fill="#7a0177")+theme_classic()+ 
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0, 1),labels=NULL) +
  scale_x_discrete(limits=c("<0.1","0.1-0.2","0.2-0.3",">0.3"),labels=c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4"))+xlab("")+ ylab("")
ggarrange(hfa,hfb,hfc,nrow=3, ncol=1, heights=c(1,1,1.1))

windows()
gains=ggarrange(b20a,hfa,tda,g3a,b20b,hfb,tdb,g3b,b20c,hfc,tdc,g3c, labels=c("A","B","C","D","E","F","G","H","I","J","K","L"),nrow=3,ncol=4, widths=c(1.2,1,1,1),font.label = list(size = 10, color = "black", face = "bold"))
annotate_figure(gains,
                bottom=text_grob("Difference in probability"),
                left=text_grob("Proportion of reefs passing the:", rot=90))
