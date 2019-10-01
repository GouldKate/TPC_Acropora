##Photosynthesis and Respiration code

rm(list=ls())


library(tidyverse)
library(magrittr)
library(ggplot2)

library(nls.multstart)
library(broom)
library(purrr)
library(plyr)
library(dplyr)
library(boot)

##Install packages
# load packages
library(nls.multstart)
library(broom)
library(purrr)
library(tidyverse)
library(nlstools)
library(proto)
library(nls2)
library(here)
library(lme4)
library(lmerTest)
library(data.table)

#set wd
# get the file path

setwd( "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output")
getwd()


#resp.datafin2<-read.csv("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/final_NP_R_GP_edited.csv")

# added Fragment.ID.full to final_NP_R_GP_Edited.csv and saved as  FINAL_resp_data)
#read in your rawdata
raw_data <- read.csv("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/Data_forTPC.csv")#read in Gross Photosynthesis, respiration rates, and Pnet
View(raw_data)



 #filtering out NP, removes it from the list 
mydata<-raw_data
mydata<-filter(mydata, rate.type !="NP") 
View(mydata)



#convert your C temperature to K
mydata%<>%
  mutate(K=mydata$Temp.C + 273.15)


mydataP<-subset(mydata,rate.type=="GP")
mydataR<-subset(mydata,rate.type=="R")

#### Take absolute value of negative Respiration rates umol.cm2.hr then log both rates for schoolfield equation and NLS fit
#### 

mydataR$umol.cm2.hr<-abs(mydataR$umol.cm2.hr)

Mydata<-rbind(mydataP,mydataR)



Mydata$log.rate <- log(Mydata$umol.cm2.hr + 1)  #logging and adding 0.3(-2 was smallest in data set) because a log of zero does not exist



filtdata<-Mydata
##############################################-----



#Define the Schoolfield equation:
schoolfield_high <- function(lnc, E, Eh, Th, temp, Tc) {
  Tc <- 273.15 + Tc
  k <- 8.62e-5
  boltzmann.term <- lnc + log(exp(E/k*(1/Tc - 1/temp)))
  inactivation.term <- log(1/(1 + exp(Eh/k*(1/Th - 1/temp))))
  return(boltzmann.term + inactivation.term)
}

#step one subset data: first make subsets for all of your treatment groups 
#-----subset data--------

# by genotype, K2, M5, M6# 

K2.df<-filtdata%>%
    filter(genotype=="K2")
M5.df<-filtdata%>%
    filter(genotype=="M5")
M6.df<-filtdata%>%
    filter(genotype=="M6")

#by Rate.Type GP / R#

#---- K2------
K2.df.GP<-K2.df%>%
  filter(rate.type=="GP")
K2.df.R<-K2.df%>%
  filter(rate.type=="R")

#----M5------
M5.df.GP<-M5.df%>%
  filter(rate.type=="GP")
M5.df.R<-M5.df%>%
  filter(rate.type=="R")

#------M6------
M6.df.GP<-M6.df%>%
  filter(rate.type=="GP")
M6.df.R<-M6.df%>%
  filter(rate.type=="R")


#-finally by treatment T1, T2, C & T3

K2.df.GP.T1<-K2.df.GP%>%filter(treatment=="T1")
M5.df.GP.T1<-M5.df.GP%>%filter(treatment=="T1")
M6.df.GP.T1<-M6.df.GP%>%filter(treatment=="T1")
K2.df.GP.T2<-K2.df.GP%>%filter(treatment=="T2")
M5.df.GP.T2<-M5.df.GP%>%filter(treatment=="T2")
M6.df.GP.T2<-M6.df.GP%>%filter(treatment=="T2")
 K2.df.GP.C<-K2.df.GP%>%filter(treatment=="C")
 M5.df.GP.C<-M5.df.GP%>%filter(treatment=="C")
 M6.df.GP.C<-M6.df.GP%>%filter(treatment=="C")
K2.df.GP.T3<-K2.df.GP%>%filter(treatment=="T3")
M5.df.GP.T3<-M5.df.GP%>%filter(treatment=="T3")
M6.df.GP.T3<-M6.df.GP%>%filter(treatment=="T3")


K2.df.R.T1<-K2.df.R%>%filter(treatment=="T1")
M5.df.R.T1<-M5.df.R%>%filter(treatment=="T1")
M6.df.R.T1<-M6.df.R%>%filter(treatment=="T1")
K2.df.R.T2<-K2.df.R%>%filter(treatment=="T2")
M5.df.R.T2<-M5.df.R%>%filter(treatment=="T2")
M6.df.R.T2<-M6.df.R%>%filter(treatment=="T2")
 K2.df.R.C<-K2.df.R%>%filter(treatment=="C")
 M5.df.R.C<-M5.df.R%>%filter(treatment=="C")
 M6.df.R.C<-M6.df.R%>%filter(treatment=="C")
K2.df.R.T3<-K2.df.R%>%filter(treatment=="T3")
M5.df.R.T3<-M5.df.R%>%filter(treatment=="T3")
M6.df.R.T3<-M6.df.R%>%filter(treatment=="T3")

  
#-------calculate parameter estimates per genotype-------
#Step one: make empty DF to fill with fitted values and run function
##"""it is very important to empty the dataframes before running the function multiple times"""

All.fittedT<-data.frame()#create an empty df to fill with fitted values
All.CIT<-data.frame()#create an empty df. to fill with confidence intervals

mult.fit.curves<-function(Data){
  fit2 <- nls_multstart(log.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 27),
                        data = Data,
                        iter = 500,
                        start_lower = c(lnc = -10, E = 0.1, Eh = 0.2, Th = 285),
                        start_upper = c(lnc = 10, E = 2, Eh = 5, Th = 330),
                        supp_errors = 'Y',
                        na.action = na.omit,
                        lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))
  #print(fit2)
  preds <- augment(fit2)
  Data%<>%
    mutate(fitted=preds$.fitted,
           residuals=preds$.resid)
  
  All.fittedT<<-rbind(All.fittedT,Data)
  
  fit_boots <- Data %>% 
    modelr::bootstrap(n = 10, id = 'boot_num') %>% #change the number for more straps
    group_by(boot_num) %>%
    mutate(fit = map(strap, ~nls_multstart(log.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 27),
                                           data = data.frame(.),
                                           iter = 100,
                                           start_lower = c(lnc = -10, E = 0.1, Eh = 0.2, Th = 285),
                                           start_upper = c(lnc = 10, E = 2, Eh = 5, Th = 330),
                                           lower = c(lnc=-10, E=0, Eh=0, Th=0),
                                           supp_errors = 'Y')
    ))
  fit_boots
  
  # get predictions
  preds_boot <- fit_boots %>%
    tidyr::  unnest(fit %>% map(augment)) %>%
    ungroup()
  
  new_preds <- Data %>%
    do(., data.frame(K = seq(min(.$K), max(.$K), length.out = 250), stringsAsFactors = FALSE))
  
  preds <- augment(fit2, newdata = new_preds)
  
  df.genotype<-(Data$genotype[1])
  df.rate.type<-(Data$rate.type[1])
  df.treatment<-(Data$treatment[1])
 
  
  preds <- fit_boots %>%
   tidyr:: unnest(fit %>% map(augment, newdata = new_preds)) %>%
    # group by each value of K and get quantiles
   dplyr:: group_by(., K) %>%
    summarise(lwr_CI = quantile(.fitted, 0.025),
              upr_CI = quantile(.fitted, 0.975)) %>%
    ungroup() %>%
    merge(., preds, by = 'K')%>%
    mutate(genotype=factor(df.genotype),
           treatment=factor(df.treatment),
           rate.type=factor(df.rate.type))
  All.CIT<<-rbind(All.CIT, preds)
  
}
#"""it is very important to empty the dataframes before running the function multiple times"""


#-----put all the subsetted data through the function--------
#First run all of these and get All.fitted. and All. CI

  
 mult.fit.curves(K2.df.GP.T1)
 mult.fit.curves(M5.df.GP.T1)
 mult.fit.curves(M6.df.GP.T1)
 mult.fit.curves(K2.df.GP.T2)
 mult.fit.curves(M5.df.GP.T2)
 mult.fit.curves(M6.df.GP.T2)
 mult.fit.curves(K2.df.GP.C)
 mult.fit.curves(M5.df.GP.C)
 mult.fit.curves(M6.df.GP.C)
 mult.fit.curves(K2.df.GP.T3)
 mult.fit.curves(M5.df.GP.T3)
 mult.fit.curves(M6.df.GP.T3)
 
 mult.fit.curves(K2.df.R.T1)
 mult.fit.curves(M5.df.R.T1)
 mult.fit.curves(M6.df.R.T1)
 mult.fit.curves(K2.df.R.T2)
 mult.fit.curves(M5.df.R.T2)
 mult.fit.curves(M6.df.R.T2)
 mult.fit.curves(K2.df.R.C)
 mult.fit.curves(M5.df.R.C)
 mult.fit.curves(M6.df.R.C)
 mult.fit.curves(K2.df.R.T3)
 mult.fit.curves(M5.df.R.T3)
 mult.fit.curves(M6.df.R.T3)
 
  


#----Graph the outputs--------
#graph with your wobbly predictions
All.fittedr<-All.fittedT%<>%
  mutate(group=factor(paste(rate.type,genotype)))

All.CIr<-All.CIT%<>%
  mutate(group=factor(paste(rate.type,genotype)))


#-----Save plot in PP-------##
library(esquisse)
library(rvg)

TPCall<-ggplot() +
  geom_ribbon(data=subset(All.CIT, lwr_CI>0),
              aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedT, aes(x=(K - 273.15), y=log.rate, shape=genotype)) +
  geom_line(data=All.fittedT, aes(x=(K - 273.15), y=fitted, colour=group)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~treatment)
esquisse:::ggplot_to_ppt("TPCall")



#----Graph the outputs--------
#graph with your wobbly predictions
All.fittedg<-All.fittedT%<>%
  mutate(group=factor(paste(treatment,rate.type)))

All.CIg<-All.CIT%<>%
  mutate(group=factor(paste(treatment,rate.type)))
####
TPCallg<-ggplot() +
  geom_ribbon(data=subset(All.CIg, lwr_CI>0),
              aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedg, aes(x=(K - 273.15), y=log.rate, shape=treatment)) +
  geom_line(data=All.fittedg, aes(x=(K - 273.15), y=fitted, colour=group)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~genotype)
esquisse:::ggplot_to_ppt("TPCallg")

#----Graph the outputs--------
#graph with your wobbly predictions
All.fittedi<-All.fittedT%<>%
  mutate(group=factor(paste(rate.type,individual.ID)))

All.CIi<-All.CIT%<>%
  mutate(group=factor(paste(rate.type,genotype)))
####
TPCallg<-ggplot() +
  geom_ribbon(data=subset(All.CIg, lwr_CI>0),
              aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedg, aes(x=(K - 273.15), y=log.rate, shape=treatment)) +
  geom_line(data=All.fittedg, aes(x=(K - 273.15), y=fitted, colour=group)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~genotype)
esquisse:::ggplot_to_ppt("TPCallg")


#graph using a smoothing function
STPCall<-ggplot() +
  geom_ribbon(data=subset(All.CIT, lwr_CI>0),aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedT, aes(x=(K - 273.15), y=log.rate, shape=genotype)) +
  geom_smooth(data=All.fittedT, aes(x=(K - 273.15), y=fitted, colour=group, se=F)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~treatment)

esquisse:::ggplot_to_ppt("STPCall")

#graph using a smoothing function
STPCallg<-ggplot() +
  geom_ribbon(data=subset(All.CIg, lwr_CI>0),aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedg, aes(x=(K - 273.15), y=log.rate, shape=treatment)) +
  geom_smooth(data=All.fittedg, aes(x=(K - 273.15), y=fitted, colour=group,se=F)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~genotype)

esquisse:::ggplot_to_ppt("STPCallg")


#graph using a smoothing function and fancy up this graphical representation of biological processes
TPC<-ggplot() +
  geom_ribbon(data=subset(All.CIT, lwr_CI>0),aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedT, aes(x=(K - 273.15), y=log.rate, shape=genotype))+ 
  geom_smooth(data=All.fittedT, aes(x=(K - 273.15), y=fitted, colour=group,se=F)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  facet_wrap(~treatment)+
  theme_classic()+
  labs(title="Thermal Performance Curves Across treatments")
esquisse:::ggplot_to_ppt("TPC")

TPCg<-ggplot() +
  geom_ribbon(data=subset(All.CIg, lwr_CI>0),aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedg, aes(x=(K - 273.15), y=log.rate, shape=treatment))+ 
  geom_smooth(data=All.fittedg, aes(x=(K - 273.15), y=fitted, colour=group,se=F)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  facet_wrap(~genotype)+
  theme_classic()+
  labs(title="Thermal Performance Curves Across treatments")
esquisse:::ggplot_to_ppt("TPCg")




#### How to find Cmin and Cmax

###------Blank the All.parameters to get subsetted data to genotype for Topt calculation and ANOVA comparison-----

#----run this second to get parameters for each genotype
All.parameters<-data.frame()#create an empty df to fill with calculated parameters


mult.fit.curves2<-function(Data){
  fit2 <- nls_multstart(log.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 27),
                        data = Data,
                        iter = 500,
                        start_lower = c(lnc = -10, E = 0.1, Eh = 0.2, Th = 285),
                        start_upper = c(lnc = 10, E = 2, Eh = 5, Th = 330),
                        supp_errors = 'Y',
                        na.action = na.omit,
                        lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))
  #print(fit2)
  params <- tidy(fit2)
  
  params%<>%
    mutate(genotype=Data$genotype[1],
           treatment=Data$treatment[1],
           rate.type=Data$rate.type[1],
           temp.Cat= Data$temp.Cat[1],
           individual.ID=Data$individual.ID[1])
  print(params)
  All.parameters<<-rbind(params, All.parameters)
}

#Parameters for Topt and Pmax
#(subset
mult.fit.curves2(subset(K2.df.GP.T1, individual.ID== "K201a"))
mult.fit.curves2(subset(M5.df.GP.T1, individual.ID== "M501a"))
mult.fit.curves2(subset(M6.df.GP.T1, individual.ID== "M601a"))
mult.fit.curves2(subset(K2.df.GP.T2, individual.ID== "K202a"))
mult.fit.curves2(subset(M5.df.GP.T2, individual.ID== "M502a"))
mult.fit.curves2(subset(M6.df.GP.T2, individual.ID== "M602a"))
mult.fit.curves2(subset(K2.df.GP.C, individual.ID== "K2Ca"))
mult.fit.curves2(subset(M5.df.GP.C, individual.ID== "M5Ca"))
mult.fit.curves2(subset(M6.df.GP.C, individual.ID== "M6Ca"))
mult.fit.curves2(subset(K2.df.GP.T3, individual.ID== "K203a"))
mult.fit.curves2(subset(M5.df.GP.T3, individual.ID== "M503a"))
mult.fit.curves2(subset(M6.df.GP.T3, individual.ID== "M603a"))
mult.fit.curves2(subset(K2.df.R.T1, individual.ID==  "K201a"))
mult.fit.curves2(subset(M5.df.R.T1, individual.ID==  "M501a"))
mult.fit.curves2(subset(M6.df.R.T1, individual.ID==  "M601a"))
mult.fit.curves2(subset(K2.df.R.T2, individual.ID==  "K202a"))
mult.fit.curves2(subset(M5.df.R.T2, individual.ID==  "M502a"))
mult.fit.curves2(subset(M6.df.R.T2, individual.ID==  "M602a"))
mult.fit.curves2(subset(K2.df.R.C, individual.ID==  "K20Ca"))
mult.fit.curves2(subset(M5.df.R.C, individual.ID==  "M50Ca"))
mult.fit.curves2(subset(M6.df.R.C, individual.ID==  "M60Ca"))
mult.fit.curves2(subset(K2.df.R.T3, individual.ID==  "K203a"))
mult.fit.curves2(subset(M5.df.R.T3, individual.ID==  "M503a"))
mult.fit.curves2(subset(M6.df.R.T3, individual.ID==  "M603a"))
mult.fit.curves2(subset(K2.df.GP.T1, individual.ID== "K201b"))
mult.fit.curves2(subset(M5.df.GP.T1, individual.ID== "M501b"))
mult.fit.curves2(subset(M6.df.GP.T1, individual.ID== "M601b"))
mult.fit.curves2(subset(K2.df.GP.T2, individual.ID== "K202b"))
mult.fit.curves2(subset(M5.df.GP.T2, individual.ID== "M502b"))
mult.fit.curves2(subset(M6.df.GP.T2, individual.ID== "M602b"))
mult.fit.curves2(subset(K2.df.GP.C, individual.ID== "K20Cb"))
mult.fit.curves2(subset(M5.df.GP.C, individual.ID== "M50Cb"))
mult.fit.curves2(subset(M6.df.GP.C, individual.ID== "M60Cb"))
mult.fit.curves2(subset(K2.df.GP.T3, individual.ID== "K203b"))
mult.fit.curves2(subset(M5.df.GP.T3, individual.ID== "M503b"))
mult.fit.curves2(subset(M6.df.GP.T3, individual.ID== "M603b"))
mult.fit.curves2(subset(K2.df.R.T1, individual.ID== "K201b"))
mult.fit.curves2(subset(M5.df.R.T1, individual.ID== "M501b"))
mult.fit.curves2(subset(M6.df.R.T1, individual.ID== "M601b"))
mult.fit.curves2(subset(K2.df.R.T2, individual.ID== "K202b"))
mult.fit.curves2(subset(M5.df.R.T2, individual.ID== "M502b"))
mult.fit.curves2(subset(M6.df.R.T2, individual.ID== "M602b"))
mult.fit.curves2(subset(K2.df.R.C, individual.ID== "K20Cb"))
mult.fit.curves2(subset(M5.df.R.C, individual.ID== "M50Cb"))
mult.fit.curves2(subset(M6.df.R.C, individual.ID== "M60Cb"))
mult.fit.curves2(subset(K2.df.R.T3, individual.ID==  "K203b"))
mult.fit.curves2(subset(M5.df.R.T3, individual.ID==  "M503b"))
mult.fit.curves2(subset(M6.df.R.T3, individual.ID==  "M603b"))
mult.fit.curves2(subset(K2.df.GP.T1, individual.ID== "K201c"))
mult.fit.curves2(subset(M5.df.GP.T1, individual.ID== "M501c"))
mult.fit.curves2(subset(M6.df.GP.T1, individual.ID== "M601c"))
mult.fit.curves2(subset(K2.df.GP.T2, individual.ID== "K202c"))
mult.fit.curves2(subset(M5.df.GP.T2, individual.ID== "M502c"))
mult.fit.curves2(subset(M6.df.GP.T2, individual.ID== "M602c"))
mult.fit.curves2(subset(K2.df.GP.C, individual.ID== "K20Cc"))
mult.fit.curves2(subset(M5.df.GP.C, individual.ID== "M50Cc"))
mult.fit.curves2(subset(M6.df.GP.C, individual.ID== "M60Cc"))
mult.fit.curves2(subset(K2.df.GP.T3, individual.ID== "K203c"))
mult.fit.curves2(subset(M5.df.GP.T3, individual.ID== "M503c"))
mult.fit.curves2(subset(M6.df.GP.T3, individual.ID== "M603c"))
mult.fit.curves2(subset(K2.df.R.T1, individual.ID==  "K201c"))
mult.fit.curves2(subset(M5.df.R.T1, individual.ID==  "M501c"))
mult.fit.curves2(subset(M6.df.R.T1, individual.ID==  "M601c"))
mult.fit.curves2(subset(K2.df.R.T2, individual.ID==  "K202c"))
mult.fit.curves2(subset(M5.df.R.T2, individual.ID==  "M502c"))
mult.fit.curves2(subset(M6.df.R.T2, individual.ID==  "M602c"))
mult.fit.curves2(subset(K2.df.R.C, individual.ID== "K20Cc"))
mult.fit.curves2(subset(M5.df.R.C, individual.ID== "M50Cc"))
mult.fit.curves2(subset(M6.df.R.C, individual.ID== "M60Cc"))
mult.fit.curves2(subset(K2.df.R.T3, individual.ID==  "K203c"))
mult.fit.curves2(subset(M5.df.R.T3, individual.ID==  "M503c"))
mult.fit.curves2(subset(M6.df.R.T3, individual.ID==  "M603c"))
mult.fit.curves2(subset(K2.df.GP.T1, individual.ID== "K201d"))
mult.fit.curves2(subset(M5.df.GP.T1, individual.ID== "M501d"))
mult.fit.curves2(subset(M6.df.GP.T1, individual.ID== "M601d"))
mult.fit.curves2(subset(K2.df.GP.T2, individual.ID== "K202d"))
mult.fit.curves2(subset(M5.df.GP.T2, individual.ID== "M502d"))
mult.fit.curves2(subset(M6.df.GP.T2, individual.ID== "M602d"))
mult.fit.curves2(subset(K2.df.GP.C, individual.ID== "K20Cd"))
mult.fit.curves2(subset(M5.df.GP.C, individual.ID== "M50Cd"))
mult.fit.curves2(subset(M6.df.GP.C, individual.ID== "M60Cd"))
mult.fit.curves2(subset(K2.df.GP.T3, individual.ID== "K203d"))
mult.fit.curves2(subset(M5.df.GP.T3, individual.ID== "M503d"))
mult.fit.curves2(subset(M6.df.GP.T3, individual.ID== "M603d"))
mult.fit.curves2(subset(K2.df.R.T1, individual.ID==  "K201d"))
mult.fit.curves2(subset(M5.df.R.T1, individual.ID==  "M501d"))
mult.fit.curves2(subset(M6.df.R.T1, individual.ID==  "M601d"))
mult.fit.curves2(subset(K2.df.R.T2, individual.ID==  "K202d"))
mult.fit.curves2(subset(M5.df.R.T2, individual.ID==  "M502d"))
mult.fit.curves2(subset(M6.df.R.T2, individual.ID==  "M602d"))
mult.fit.curves2(subset(K2.df.R.C, individual.ID== "K20Cd"))
mult.fit.curves2(subset(M5.df.R.C, individual.ID== "M50Cd"))
mult.fit.curves2(subset(M6.df.R.C, individual.ID== "M60Cd"))
mult.fit.curves2(subset(K2.df.R.T3, individual.ID==  "K203d"))
mult.fit.curves2(subset(M5.df.R.T3, individual.ID==  "M503d"))
mult.fit.curves2(subset(M6.df.R.T3, individual.ID==  "M603d"))


#mult.fit.curves2(K2.df.GP.T1)
#mult.fit.curves2(M5.df.GP.T1)
#mult.fit.curves2(M6.df.GP.T1)
#mult.fit.curves2(K2.df.GP.T2)
#mult.fit.curves2(M5.df.GP.T2)
#mult.fit.curves2(M6.df.GP.T2)
#mult.fit.curves2(K2.df.GP.C)
#mult.fit.curves2(M5.df.GP.C)
#mult.fit.curves2(M6.df.GP.C)
#mult.fit.curves2(K2.df.GP.T3)
#mult.fit.curves2(M5.df.GP.T3)
#mult.fit.curves2(M6.df.GP.T3)
#
#mult.fit.curves2(K2.df.R.T1)
#mult.fit.curves2(M5.df.R.T1)
#mult.fit.curves2(M6.df.R.T1)
#mult.fit.curves2(K2.df.R.T2)
#mult.fit.curves2(M5.df.R.T2)
#mult.fit.curves2(M6.df.R.T2)
#mult.fit.curves2(K2.df.R.C)
#mult.fit.curves2(M5.df.R.C)
#mult.fit.curves2(M6.df.R.C)
#mult.fit.curves2(K2.df.R.T3)
#mult.fit.curves2(M5.df.R.T3)
#mult.fit.curves2(M6.df.R.T3)
  



# plot mean of parameters  across all genotypes (not very attractive)
MeanParams<-left_join(All.CIT, All.parameters)
View(MeanParams)


para<-ggplot(MeanParams, aes(col = rate.type)) +
  geom_point(aes(genotype, estimate)) +
  facet_wrap(~ term, scale = 'free_x', ncol = 4) +
  geom_linerange(aes(genotype, ymin = lwr_CI, ymax = upr_CI)) +
  coord_flip() +
  scale_color_manual(values = c('green4', 'black')) +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  theme(legend.position = 'top') +
  xlab('curve') +
  ylab('parameter estimate')#+
facet_wrap(rate.type)


esquisse:::ggplot_to_ppt("para")


###################### Obtaining Topt


#48 GROUPS
#192 ESTIAMTES (4 x 48)
#96 estimates/4 treatments= 24/2 (light and Dark)= 12/temperature /3 genotypes= 4 (treatments)
#------T opt from parameters-------
get_topt<-function(E, Th, Eh){
  return((Eh*Th)/(Eh+(8.62e-05*Th*(log((Eh/E))-1))))
}

topt_data<-All.parameters

topt_data<-All.parameters%>%
  select(individual.ID, genotype, rate.type, treatment, estimate, term)%>%
  group_by(individual.ID, genotype, rate.type, treatment)%>%
  spread(term,estimate)

write.csv(topt_data, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/TPC_param_est_ind.csv")


topt_data$topt<-get_topt(E=topt_data$E, Th=topt_data$Th, Eh=topt_data$Eh)

topt_data$topt_C<-topt_data$topt-273.15

write.csv(topt_data, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/Topt_estimates_ind.csv")

#-------linear Models does treatment or genotype affect Topt----
library(car)
library(lme4)
library(stats)
library(Rmisc)
library(Hmisc)

topt_nonegGP<-subset(topt_data, rate.type=="GP")
topt_nonegGP<-topt_nonegGP%>%
  filter(topt_C >0)
#M602a, M602c, M60Cb

topt_nonegR<-subset(topt_data, rate.type=="R")
topt_nonegR<-topt_nonegR%>%
  filter(topt_C<0)
#K201b, K20Cb




topt_noneg<-rbind(topt_nonegGP,topt_nonegR)
View(topt_noneg)

write.csv(topt_noneg, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/Topt_estimates_noneg.csv")


# plot distribution of Topt ##### ------ Parameter graphs--- 
Toptmeans<-ggplot(subset(topt_noneg, rate.type=="GP"), aes(topt_C)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~rate.type)+
  xlab('Optimum Temperature (ºC)') +
  ggtitle('Distribution of optimum temperatures')

esquisse:::ggplot_to_ppt("Toptmeans")

Toptmeans<-ggplot(subset(topt_noneg, rate.type=="R"), aes(topt_C)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~rate.type)+
  xlab('Optimum Temperature (ºC)') +
  ggtitle('Distribution of optimum temperatures')

esquisse:::ggplot_to_ppt("Toptmeans")


# plot distribution of estimated parameters E, Eh, lnc, Th
p1 <- ggplot(topt_noneg, aes(E)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ rate.type, scales = 'free_x')
p1
esquisse:::ggplot_to_ppt("p1")
# plot distribution of estimated parameters E, Eh, lnc, Th
p2 <- ggplot(topt_noneg, aes(Eh)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ rate.type, scales = 'free_x')
p2
esquisse:::ggplot_to_ppt("p2")

# plot distribution of estimated parameters E, Eh, lnc, Th
p3 <- ggplot(topt_noneg, aes(lnc)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ rate.type)
p3
esquisse:::ggplot_to_ppt("p3")

# plot distribution of estimated parameters E, Eh, lnc, Th
p4 <- ggplot(subset(topt_noneg, rate.type=="GP"), aes(Th)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ rate.type)
p4
esquisse:::ggplot_to_ppt("p4")

p4r <- ggplot(subset(topt_noneg, rate.type=="R"), aes(Th)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ rate.type)
esquisse:::ggplot_to_ppt("p4r")
# plot distribution of estimated parameters E, Eh, lnc, Th by genotype
# 
p6T <- ggplot(subset(topt_noneg, rate.type=="GP"), aes(topt_C)) +
geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ genotype, scales = 'free_x')
p6
esquisse:::ggplot_to_ppt("p6T")

p6r <- ggplot(subset(topt_noneg, rate.type=="R"), aes(topt_C)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ genotype, scales = 'free_x')
p6
esquisse:::ggplot_to_ppt("p6r")

p6Tt <- ggplot(subset(topt_noneg, rate.type=="GP"), aes(topt_C)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ treatment, scales = 'free_x')
p6
esquisse:::ggplot_to_ppt("p6Tt")

p6rt <- ggplot(subset(topt_noneg, rate.type=="R"), aes(topt_C)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ treatment, scales = 'free_x')
p6
esquisse:::ggplot_to_ppt("p6rt")

p6 <- ggplot(topt_noneg, aes(E)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ genotype, scales = 'free_x')
p6
esquisse:::ggplot_to_ppt("p6")
# plot distribution of estimated parameters E, Eh, lnc, Th
p7 <- ggplot(topt_noneg, aes(Eh)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ genotype, scales = 'free_x')
p7
esquisse:::ggplot_to_ppt("p7")

# plot distribution of estimated parameters E, Eh, lnc, Th
p8<- ggplot(topt_noneg, aes(lnc)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ genotype)
p8
esquisse:::ggplot_to_ppt("p8")

# plot distribution of estimated parameters E, Eh, lnc, Th 
p5G <- ggplot(subset(topt_noneg, rate.type=="GP"), aes(Th)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ genotype)
p5
esquisse:::ggplot_to_ppt("p5G")

p5r <- ggplot(subset(topt_noneg, rate.type=="R"), aes(Th)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ genotype)
p5
esquisse:::ggplot_to_ppt("p5r")



meanDF2<-summary(topt_noneg)
View(meanDF2)
summary(topt_noneg)

#individual.ID genotype rate.type treatment       E                  Eh                 lnc       
#K201a  : 2    K2:29    GP:41     C :17     Min.   :0.000000   Min.   :  0.00000   Min.   :1.008  
#K201c  : 2    M5:29    R :45     T1:23     1st Qu.:0.006218   1st Qu.:  0.00000   1st Qu.:1.362  
#K201d  : 2    M6:28              T2:22     Median :0.040731   Median :  0.07545   Median :1.594  
#K202a  : 2                       T3:24     Mean   :0.228450   Mean   :  8.18734   Mean   :1.641  
#K202b  : 2                                 3rd Qu.:0.110749   3rd Qu.: 12.77611   3rd Qu.:1.641  
#K202c  : 2                                 Max.   :4.815770   Max.   :174.04987   Max.   :5.204  
#(Other):74                                                                                       
#Th                topt               topt_C       
#Min.   :     0.1   Min.   : -0.02846   Min.   :-273.18  
#1st Qu.:   299.1   1st Qu.: -0.00003   1st Qu.:-273.15  
#Median :   306.8   Median : 37.68001   Median :-235.47  
#Mean   :  4363.6   Mean   :147.59133   Mean   :-125.56  
#3rd Qu.:   307.4   3rd Qu.:304.61827   3rd Qu.:  31.47  
#Max.   :222026.8   Max.   :306.15136   Max.   :  33.00




### Table of Averages per genotype

topt_nonegAVG<-topt_noneg

topt_nonegAVG$avgE<-topt_nonegAVG$E

data.summarytoptE<-topt_noneg %>%
  group_by(genotype, treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(meanE=mean(E), seE=sd(E)/sqrt(n()))#%>%
data.summarytoptEh<-topt_noneg %>%
  group_by(genotype, treatment, rate.type) %>% 
  dplyr::summarise(meanEh=mean(Eh), seEh=sd(Eh)/sqrt(n()))
data.summarytoptlnc<-topt_noneg %>%
  group_by(genotype, treatment, rate.type) %>% 
  dplyr::summarise(meanlnc=mean(lnc), selnc=sd(lnc)/sqrt(n()))
data.summarytoptTh<-topt_noneg %>%
  group_by(genotype, treatment, rate.type) %>% 
  dplyr::summarise(meanTh=mean(Th), seTh=sd(Th)/sqrt(n()))
data.summarytoptTopt<-topt_noneg %>%
  group_by(genotype, treatment, rate.type) %>% 
  dplyr::summarise(meanTopt=mean(topt_C), seTopt=sd(topt_C)/sqrt(n()))

Allsum<-left_join(data.summarytoptTopt, data.summarytoptE)
Allsum2<-left_join(Allsum, data.summarytoptTh)
Allsum3<-left_join(Allsum2,data.summarytoptEh)
Allsumfinal<-left_join(Allsum3, data.summarytoptlnc)

write.csv(Allsumfinal, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/TPC_param_mean_se.csv")

all<-Allsumfinal


#visualize groups GP
#
GPraw<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=treatment, topt_C, colour=treatment))+
  geom_boxplot()+
  geom_point()+
  ylab("Topt (ºC)")

esquisse:::ggplot_to_ppt("GPraw")

GPoptS<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=genotype, topt_C, colour=genotype))+
  geom_boxplot()+
  geom_point()+
  ylab("Topt (ºC)")
esquisse::ggplot_to_ppt("GPoptS")

GPtbS<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=genotype, topt_C, colour=treatment, shape=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Topt ºC)")+
  labs(title="Mean Topt across Genotypes (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("GPtbS")

GPtbS<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=treatment, topt_C, colour=genotype, shape=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Topt ºC)")+
  labs(title="Mean Topt across Treatments (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("GPtbS")

library(car)

topt_nonegRep<-read.csv("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/Topt_estimates_noneg_repl.csv")


#linear model
gp.mod<-Anova(lm(topt_C~treatment+genotype, data=subset(topt_noneg,rate.type=="GP")))
summary(gp.mod)
gp.mod


# No significant differnece in Topt across treatment or genotype
m1N<-aov(topt_C~treatment+genotype, data=subset(topt_noneg,rate.type=="GP"))
TukeyHSD(m1N)
summary(m1N)

### No difference in Topt
### 
##visualize groups R

RTopt<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=treatment, topt_C, colour=treatment))+
  geom_boxplot()+
  geom_point()+
  ylab("Topt (ºC)")

esquisse:::ggplot_to_ppt("RTopt")

RtoptS<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=genotype, topt_C, colour=genotype))+
  geom_boxplot()+
  geom_point()+
  ylab("Topt (ºC)")
esquisse::ggplot_to_ppt("RtoptS")

RtbT<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=treatment, topt_C, colour=genotype, shape=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Topt (ºC)")+
  labs(title="Mean Topt across Treatment (R)")+
  coord_flip()

esquisse:::ggplot_to_ppt("RtbT")

RtbS<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=genotype, topt_C, colour=treatment, shape=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Topt (ºC)")+
  labs(title="Mean Topt across genotype (R)")+
  coord_flip()

esquisse:::ggplot_to_ppt("RtbS")
TukeyHSD(m2N)

#linear model
R.mod<-Anova(lm(topt_C~treatment+genotype, data=subset(topt_noneg,rate.type=="R")))
summary(R.mod)
R.mod

m2N<-aov(topt_C~treatment+genotype, data=subset(topt_noneg,rate.type=="R")) #treatment significant (PAST) 0.0108
TukeyHSD(m2N, "treatment")
TukeyHSD(m2N, "genotype")

TukeyHSD(m2N)

K2modT<-aov(topt_C~treatment, data=subset(topt_noneg,genotype=="K2")) # treatment significant across depth for PAST 0.0234
summary(K2modT)
TukeyHSD(K2modT)

K2mod<-aov(topt_C~replicas, data=subset(topt_nonegRep,genotype=="K2")) # treatment significant across depth for PAST 0.0234
summary(K2mod)
TukeyHSD(K2mod)

M5modT<-aov(topt_C~treatment, data=subset(topt_noneg,genotype=="M5")) # treatment significant across depth for PAST 0.0234
summary(M5modT)
TukeyHSD(M5modT)

M5mod<-aov(topt_C~replicas, data=subset(topt_nonegRep,genotype=="M5")) # treatment significant across depth for PAST 0.0234
summary(M5mod)
TukeyHSD(M5mod)

M6modT<-aov(topt_C~treatment, data=subset(topt_noneg,genotype=="M6")) # treatment significant across depth for PAST 0.0234
summary(M6modT)
TukeyHSD(M6modT)

M6mod<-aov(topt_C~replicas, data=subset(topt_nonegRep,genotype=="M6")) # treatment significant across depth for PAST 0.0234
summary(M6mod)
TukeyHSD(M6mod)

#-------linear models does treatment or genotype lnc----
#visualize groups GPs
GPtLT<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=treatment, lnc, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("lnc")+
  labs(title="Mean lnc across Treatment (GP)")+
  coord_flip()

esquisse:::ggplot_to_ppt("GPtLT")

GPtbLT<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=genotype, lnc, colour=treatment, shape=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("lnc")+
  labs(title="Mean lnc across Genotype (GP)")+
  coord_flip()


esquisse:::ggplot_to_ppt("GPtbLT")


#linear model
IncGP.mod<-Anova(lm(lnc~treatment+genotype, data=subset(topt_noneg,rate.type=="GP")))
summary(IncGP.mod)
IncGP.mod

m3N<-aov(lnc~treatment+genotype, data=subset(topt_noneg,rate.type=="GP"))
TukeyHSD(m3N)
m3N

#visualize groups R
RtbLT<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=treatment, lnc, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("lnc")+
  labs(title="Mean lnc across Treatment (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("RtbLT")


RtbLS<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=genotype, lnc, shape=genotype,colour=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("lnc")+
  labs(title="Mean lnc across Genotype (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("RtbLS")


#linear model
Inc.R.mod<-Anova(lm(lnc~treatment+genotype, data=subset(topt_noneg,rate.type=="R")))
Inc.R.mod
summary(Inc.R.mod)

m4N<-aov(lnc~treatment+genotype, data=subset(topt_noneg,rate.type=="R"))
TukeyHSD(m4N)
summary(m4N)




##### Statistical test for E GP

EGP<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=treatment, E, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("E")+
  labs(title="Mean E across Treatment (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("EGP")


EsGP<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=genotype, E, shape=genotype,colour=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("E")+
  labs(title="Mean E across Genotype (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("EsGP")


#linear model
E.GP.mod<-Anova(lm(E~treatment+genotype, data=subset(topt_noneg,rate.type=="GP")))
E.GP.mod
summary(E.GP.mod)

mEN<-aov(E~treatment+genotype, data=subset(topt_noneg,rate.type=="GP"))
TukeyHSD(mEN)
summary(mEN)


##### Statistical test for E R
#visualize groups R
ER<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=treatment, E, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("E")+
  labs(title="Mean E across Treatment (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("ER")


EsR<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=genotype, E, shape=genotype,colour=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("E")+
  labs(title="Mean E across Genotype (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("EsR")


#linear model
E.R.mod<-Anova(lm(E~treatment+genotype, data=subset(topt_noneg,rate.type=="R")))
E.R.mod
summary(E.R.mod)

mER<-aov(E~treatment+genotype, data=subset(topt_noneg,rate.type=="R"))
TukeyHSD(mER)
summary(mER)






##### Statistical test for Eh GP

EhGP<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=treatment, Eh, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Eh")+
  labs(title="Mean Eh across Treatment (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("EhGP")


EhsGP<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=genotype, Eh, shape=genotype,colour=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Eh")+
  labs(title="Mean Eh across Genotype (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("EhsGP")


#linear model
Eh.GP.mod<-Anova(lm(Eh~treatment+genotype, data=subset(topt_noneg,rate.type=="GP")))
Eh.GP.mod
summary(Eh.GP.mod)

mEhN<-aov(Eh~treatment+genotype, data=subset(topt_noneg,rate.type=="GP"))
TukeyHSD(mEhN)
summary(mEhN)


##### Statistical test for Eh R
#visualize groups R
EhR<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=treatment, Eh, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Eh")+
  labs(title="Mean Eh across Treatment (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("EhR")


EhsR<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=genotype, Eh, shape=genotype,colour=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Eh")+
  labs(title="Mean Eh across Genotype (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("EhsR")




#linear model
Eh.R.mod<-Anova(lm(Eh~treatment+genotype, data=subset(topt_noneg,rate.type=="R")))
Eh.R.mod
summary(Eh.R.mod)

mEhR<-aov(Eh~treatment+genotype, data=subset(topt_noneg,rate.type=="R"))
TukeyHSD(mEhR)
summary(mEhR)




##### Statistical test for Th GP

ThGP<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=treatment, Th, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Th")+
  labs(title="Mean Th across Treatment (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("ThGP")


ThsGP<-ggplot(subset(topt_noneg,rate.type=="GP"), aes(x=genotype, Th, shape=genotype,colour=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Th")+
  labs(title="Mean Th across Genotype (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("ThsGP")


#linear model
Th.GP.mod<-Anova(lm(Th~treatment+genotype, data=subset(topt_noneg,rate.type=="GP")))
Th.GP.mod
summary(Th.GP.mod)

mThN<-aov(Th~treatment+genotype, data=subset(topt_noneg,rate.type=="GP"))
TukeyHSD(mThN)
summary(mThN)


##### Statistical test for Th R
#visualize groups R
ThR<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=treatment, y=Th, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=7)+
  ylab("Th")+
  labs(title="Mean Th across Treatment (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("ThR")


ThsR<-ggplot(subset(topt_noneg,rate.type=="R"), aes(x=genotype, Th, colour=treatment,shape=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Th")+
  labs(title="Mean Th across Genotype (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("ThsR")


#linear model
Th.R.mod<-Anova(lm(Th~treatment+genotype, data=subset(topt_noneg,rate.type=="R")))
Th.R.mod
summary(Th.R.mod)

mThR<-aov(Th~treatment+genotype, data=subset(topt_noneg,rate.type=="R"))
TukeyHSD(mThR)
summary(mThR)

### PMax
#calculate Pmax values between sites
## have to make dataframe long format----- FINAL FINAL
## 
## 
Pmax1<- topt_nonegRep%>%
  select(replicas, individual.ID, rate.type, treatment, genotype, E, Eh, Th, lnc,topt, topt_C)%>%
  mutate(Pmax = schoolfield_high(lnc = lnc, E = E, Th = Th, Eh = Eh, temp = topt, Tc = 27)) %>% #add in factors that make up schoolfield function, reference topt to get pmax
  group_by(., replicas, individual.ID, rate.type, treatment, genotype)

Pmax<- topt_nonegRep%>%
  select(replicas, individual.ID, rate.type, E, Eh, Th, lnc,topt, topt_C)%>%
  mutate(Pmax = schoolfield_high(lnc = lnc, E = E, Th = Th, Eh = Eh, temp = topt, Tc = 27)) %>% #add in factors that make up schoolfield function, reference topt to get pmax
  group_by(., replicas, individual.ID, rate.type)


Pmax_datatest <- topt_noneg %>%
  #select(genotype,  term, estimate, treatment, rate.type) %>%
  #spread(term, estimate) %>%
  mutate(Pmax = schoolfield_high(lnc = lnc, E = E, Th = Th, Eh = Eh, temp = topt, Tc = 27)) %>% #add in factors that make up schoolfield function, reference topt to get pmax
  group_by(., genotype, rate.type, treatment)


write.csv(Pmax_datatest, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/June 2019/TPC_RawData/Pmax_final.geno.treat.csv")
write.csv(Pmax, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/June 2019/TPC_RawData/Pmax_final.replicas.csv")

##### PMAX------
##### Statistical test for Pmax in GP
#visualize groups R

pmaxtreat<-ggplot(subset(Pmax_datatest,rate.type=="GP"), aes(x=treatment, Pmax, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Pmax")+
  labs(title="Mean Pmax across Treatment (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("pmaxtreat")


pmaxsp<-ggplot(subset(Pmax_datatest,rate.type=="GP"), aes(x=genotype, Pmax, colour=treatment,shape=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Pmax")+
  labs(title="Mean Pmax across Genotype (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("pmaxsp")

pmaxtreat<-ggplot(subset(Pmax,rate.type=="GP"), aes(x=replicas, Pmax))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Pmax")+
  labs(title="Mean Pmax across Treatment (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("pmaxtreat")


pmaxever<-ggplot(subset(Pmax1,rate.type=="GP"), aes(x=replicas, Pmax, colour=treatment, shape= genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Pmax")+
  labs(title="Mean Pmax across Replicas (GP)")+
  coord_flip()
esquisse:::ggplot_to_ppt("pmaxever")


#linear model
Pmax.GP.mod<-Anova(lm(Pmax~treatment+genotype, data=subset(Pmax_datatest,rate.type=="GP")))
Pmax.GP.mod
summary(Pmax.GP.mod)

mPGP<-aov(Pmax~treatment+genotype, data=subset(Pmax_datatest,rate.type=="GP"))
TukeyHSD(mPGP)
summary(mPGP)



##### Statistical test for Pmax in R
#visualize groups R
#
#NEed too get rid of Negative values of Topt_
#
#
Pmax_datatest<-Pmax_datatest%>%
  filter(topt >0)


pmaxtrR<-ggplot(subset(Pmax_datatest,rate.type=="R"), aes(x=treatment, Pmax, shape=treatment, colour=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Pmax")+
  labs(title="Mean Pmax across Treatment(R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("pmaxtrR")


pmaxspR<-ggplot(subset(Pmax_datatest,rate.type=="R"), aes(x=genotype, Pmax, colour=treatment,shape=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Pmax")+
  labs(title="Mean Pmax across Genotype (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("pmaxspR")

### again filter out negatives for pmax calcuulations for Respiration
Pmax_neg<-Pmax1%>%
  filter(topt >0)

pmaxspRn<-ggplot(subset(Pmax_neg,rate.type=="R"), aes(x=replicas, Pmax, color=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Pmax")+
  labs(title="Mean Pmax across Replicas (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("pmaxspRn")


#linear model
Pmax.R.mod<-Anova(lm(Pmax~treatment+genotype, data=subset(Pmax_datatest,rate.type=="R")))
Pmax.R.mod
summary(Pmax.R.mod)

mPR<-aov(Pmax~treatment+genotype, data=subset(Pmax_datatest,rate.type=="R"))
TukeyHSD(mPR)
summary(mPR)


#anova function
Pmax.mod <- lm(Pmax~genotype*treatment, data=Pmax_datatest)
summary(Pmax.mod)
#check for normality, use normality plots

qqnorm(resid(Pmax.mod))
qqline(resid(Pmax.mod))

esquisse:::ggplot_to_ppt("NormPLot")

#check heteroscisity with boxplots

boxplot(resid(Pmax.mod)~Pmax_datatest$genotype*Pmax_datatest$treatment) 

#high R and low show inconsistent variances, may need to do weighted regression in the future

anova(Pmax.mod)
summary(Pmax.mod)
#intercept is significant- the intercept is different from zero
TukeyHSD(aov(Pmax.mod))
#compares genotype across 
anova(Pmax.mod)
summary(Pmax.mod)
TukeyHSD(aov(Pmax.mod))


library(magrittr)
library(plyr)
library(dplyr)

Pmax.means<-Pmax_datatest %>%
  group_by(genotype,treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(Pmax), se=sd(Pmax)/sqrt(n())) #calculates mean and s.e.


write.csv(Pmax.means,"C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/Pmax_negMeans.csv")



PMAX_P<-ggplot(Pmax.means, aes(x=rate.type, y=mean, col= genotype, group=factor(rate.type))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  labs(x="Rate Type", y="Maximum Rate of Performance (Pmax)", fill="treatment", color = "genotype") 

esquisse:::ggplot_to_ppt("PMAX_P")

PMAX_q<-ggplot(subset(Pmax.means, rate.type=="GP"), aes(x=genotype, y=mean, col= genotype, group=factor(rate.type))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=4) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  labs(x="genotype", y="Maximum Rate of Performance (Pmax) GP", fill="treatment", color = "genotype") 

esquisse:::ggplot_to_ppt("PMAX_q")


PMAX_R<-ggplot(subset(Pmax.means, rate.type=="R"), aes(x=genotype, y=mean, col= genotype, group=factor(rate.type))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=4) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=0.9), width=0.1) +
  labs(x="genotype", y="Maximum Rate of Performance (Pmax) R", fill="treatment", color = "genotype") 

esquisse:::ggplot_to_ppt("PMAX_R")

PMAX_P2<-ggplot(Pmax.means, aes(x=rate.type, y=mean, col= treatment, group=factor(rate.type))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=3) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position=position_dodge(width=1), width=0.5) +
  labs(x="Rate Type", y="Maximum Rate of Performance (Pmax)", fill="treatment", color = "treatment") 

esquisse:::ggplot_to_ppt("PMAX_P2")

ggsave(filename = "thermtol/Output/Pmax_graph.png", device = "png", width = 8, height = 5)

#subset data to compare across genotype- then genotype
#by genotype
DLAB.dfP<-Pmax_datatest%>%
  dplyr::filter(genotype=="DLAB")

MCAV.dfP<-Pmax_datatest%>%
  filter(genotype=="MCAV")

PAST.dfP<-Pmax_datatest%>%
  filter(genotype=="PAST")

OFRA.dfP<-Pmax_datatest%>%
  filter(genotype=="OFRA")

#by Rate.Type
#-------DLAB-------
DLAB.DF.GPP<-DLAB.dfP%>%
  filter(rate.type=="GP")

DLAB.DF.RP<-DLAB.dfP%>%
  filter(rate.type=="R")
#-------MCAV-------
MCAV.DF.GPP<-MCAV.dfP%>%
  filter(rate.type=="GP")

MCAV.DF.RP<-MCAV.dfP%>%
  filter(rate.type=="R")
#-------PAST-------
PAST.DF.GPP<-PAST.dfP%>%
  filter(rate.type=="GP")

PAST.DF.RP<-PAST.dfP%>%
  filter(rate.type=="R")
#-------OFRA-------
OFRA.DF.GPP<-OFRA.dfP%>%
  filter(rate.type=="GP")

OFRA.DF.RP<-OFRA.dfP%>%
  filter(rate.type=="R")

####### OFRA RESPIRATION -------
OFRA.DF.RP<-OFRA.dfP%>%
  filter(rate.type=="R")
#anova function
Pmax.modOFRA <- lm(Pmax~treatment, data=OFRA.DF.RP)
summary(Pmax.modOFRA)
#check for normality, use normality plots

qqnorm(resid(Pmax.modOFRA))
qqline(resid(Pmax.modOFRA))

#check heteroscisity with boxplots

boxplot(resid(Pmax.modOFRA)~OFRA.DF.RP$treatment)

#high R and low show inconsistent variances, may need to do weighted regression in the future

anova(Pmax.modOFRA)
summary(Pmax.modOFRA)
TukeyHSD(aov(Pmax.modOFRA))

data.summaryOFRA<-OFRA.DF.RP %>%
  group_by(treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(Pmax), se=sd(Pmax)/sqrt(n())) #calculates mean and s.e.

data.summaryOFRA



PMAX_POFRA<-ggplot(data.summaryOFRA, aes(x=treatment, y=mean, col= treatment, group=factor(treatment))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  labs(x="Depth", y="Maximum Rate of Performance (Pmax) OFRA (R)", fill="treatment", color = "treatment") 

PMAX_POFRA
esquisse:::ggplot_to_ppt("PMAX_POFRA")
### Topt OFRA
#anova function
Topt.modOFRA <- lm(topt_C~treatment, data=OFRA.DF.RP)
summary(Topt.modOFRA)
#check for normality, use normality plots

qqnorm(resid(Topt.modOFRA))
qqline(resid(Topt.modOFRA))

#check heteroscisity with boxplots

boxplot(resid(Topt.modOFRA)~OFRA.DF.RP$treatment)

#high R and low show inconsistent variances, may need to do weighted regression in the future

anova(Topt.modOFRA)
summary(Topt.modOFRA)
TukeyHSD(aov(Topt.modOFRA))

data.summaryTOFRA<-OFRA.DF.RP %>%
  group_by(treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(topt_C), se=sd(topt_C)/sqrt(n())) #calculates mean and s.e.
data.summaryTOFRA

Topt_POFRA<-ggplot(data.summaryTOFRA, aes(x=treatment, y=mean, col= treatment, group=factor(rate.type))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  labs(x="Depth", y="Thermal Optimum (Topt) OFRA (R)", fill="treatment", color = "treatment") 

Topt_POFRA
esquisse:::ggplot_to_ppt("Topt_POFRA")


### GP OFRA


#-------OFRA------- PHOTO------
OFRA.DF.GPP<-OFRA.dfP%>%
  filter(rate.type=="GP")
#anova function
Pmax.modOFRAP <- lm(Pmax~treatment, data=OFRA.DF.GPP)
summary(Pmax.modOFRAP)
#check for normality, use normality plots

qqnorm(resid(Pmax.modOFRAP))
qqline(resid(Pmax.modOFRAP))

#check heteroscisity with boxplots

boxplot(resid(Pmax.modOFRAP)~OFRA.DF.GPP$treatment)

#high R and low show inconsistent variances, may need to do weighted regression in the future

anova(Pmax.modOFRAP)
summary(Pmax.modOFRAP)
TukeyHSD(aov(Pmax.modOFRAP))

data.summaryOFRAGP<-OFRA.DF.GPP %>%
  group_by(treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(Pmax), se=sd(Pmax)/sqrt(n())) #calculates mean and s.e.
data.summaryOFRAGP



PMAX_POFRAGP<-ggplot(data.summaryOFRAGP, aes(x=treatment, y=mean, col= treatment, group=factor(treatment))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  labs(x="Depth", y="Maximum Rate of Performance (Pmax) OFRA (GP)", fill="treatment", color = "treatment") 

PMAX_POFRAGP
esquisse:::ggplot_to_ppt("PMAX_POFRAGP")
### Topt OFRA
#anova function
Topt.modOFRAGP <- lm(topt_C~treatment, data=OFRA.DF.GPP)
summary(Topt.modOFRAGP)
#check for normality, use normality plots

qqnorm(resid(Topt.modOFRAGP))
qqline(resid(Topt.modOFRAGP))

#check heteroscisity with boxplots

boxplot(resid(Topt.modOFRAGP)~OFRA.DF.GPP$treatment)

#high R and low show inconsistent variances, may need to do weighted regression in the future

anova(Topt.modOFRAGP)
summary(Topt.modOFRAGP)
TukeyHSD(aov(Topt.modOFRAGP))

data.summaryTOFRAGP<-OFRA.DF.GPP %>%
  group_by(treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(topt_C), se=sd(topt_C)/sqrt(n())) #calculates mean and s.e.
data.summaryTOFRAGP

Topt_POFRAGP<-ggplot(data.summaryTOFRAGP, aes(x=treatment, y=mean, col= treatment, group=factor(rate.type))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  labs(x="Depth", y="Thermal Optimum (Topt) OFRA (GP)", fill="treatment", color = "treatment") 

Topt_POFRAGP
esquisse:::ggplot_to_ppt("Topt_POFRAGP")



#-------PAST- PHOTOSYNTHESIS------
PAST.DF.GPP<-PAST.dfP%>%
  filter(rate.type=="GP")

#anova function
Pmax.modPASTP <- lm(Pmax~treatment, data=PAST.DF.GPP)
summary(Pmax.modPASTP)
#check for normality, use normality plots

qqnorm(resid(Pmax.modPASTP))
qqline(resid(Pmax.modPASTP))

#check heteroscisity with boxplots

boxplot(resid(Pmax.modPASTP)~PAST.DF.GPP$treatment)


data.summaryPASTGP<-PAST.DF.GPP %>%
  group_by(treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(Pmax), se=sd(Pmax)/sqrt(n())) #calculates mean and s.e.
data.summaryPASTGP

PMAX_PPASTGP<-ggplot(data.summaryPASTGP, aes(x=treatment, y=mean, col= treatment, group=factor(treatment))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  labs(x="Depth", y="Maximum Rate of Performance (Pmax) PAST (GP)", fill="treatment", color = "treatment") 

PMAX_PPASTGP
esquisse:::ggplot_to_ppt("PMAX_PPASTGP")
### Topt OFRA
#anova function
Topt.modPASTGP <- lm(topt_C~treatment, data=PAST.DF.GPP)
summary(Topt.modPASTGP)
#check for normality, use normality plots

qqnorm(resid(Topt.modPASTGP))
qqline(resid(Topt.modPASTGP))

#check heteroscisity with boxplots

boxplot(resid(Topt.modPASTGP)~PAST.DF.GPP$treatment)

#high R and low show inconsistent variances, may need to do weighted regression in the future


data.summaryTPASTGP<-PAST.DF.GPP %>%
  group_by(treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(topt_C), se=sd(topt_C)/sqrt(n())) #calculates mean and s.e.
data.summaryTPASTGP

Topt_PPASTGP<-ggplot(data.summaryTPASTGP, aes(x=treatment, y=mean, col= treatment, group=factor(rate.type))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  labs(x="Depth", y="Thermal Optimum (Topt) PAST (GP)", fill="treatment", color = "treatment") 

Topt_PPASTGP
esquisse:::ggplot_to_ppt("Topt_PPASTGP")

#####---- PAST RESPIRATION-----

PAST.DF.RP<-PAST.dfP%>%
  filter(rate.type=="R")

#anova function
Pmax.modPAST <- lm(Pmax~treatment, data=PAST.DF.RP)
summary(Pmax.modPAST)
#check for normality, use normality plots

qqnorm(resid(Pmax.modPAST))
qqline(resid(Pmax.modPAST))

#check heteroscisity with boxplots

boxplot(resid(Pmax.modPAST)~PAST.DF.RP$treatment)


data.summaryPAST<-PAST.DF.RP %>%
  group_by(treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(Pmax), se=sd(Pmax)/sqrt(n())) #calculates mean and s.e.
data.summaryPAST


PMAX_PPAST<-ggplot(data.summaryPAST, aes(x=treatment, y=mean, col= treatment, group=factor(treatment))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  labs(x="Depth", y="Maximum Rate of Performance (Pmax) PAST (R)", fill="treatment", color = "treatment") 

PMAX_PPAST
esquisse:::ggplot_to_ppt("PMAX_PPAST")
### Topt PAST
#anova function
Topt.modPAST <- lm(topt_C~treatment, data=PAST.DF.RP)
summary(Topt.modPAST)
#check for normality, use normality plots

qqnorm(resid(Topt.modPAST))
qqline(resid(Topt.modPAST))

#check heteroscisity with boxplots

boxplot(resid(Topt.modPAST)~PAST.DF.RP$treatment)

data.summaryTPAST<-PAST.DF.RP %>%
  group_by(treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(topt_C), se=sd(topt_C)/sqrt(n())) #calculates mean and s.e.
data.summaryTPAST

Topt_PPAST<-ggplot(data.summaryTPAST, aes(x=treatment, y=mean, col= treatment, group=factor(rate.type))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  labs(x="Depth", y="Thermal Optimum (Topt) PAST (R)", fill="treatment", color = "treatment") 

Topt_PPAST
esquisse:::ggplot_to_ppt("Topt_PPAST")
######------NEXT


#by Rate.Type
#-------DLAB PHOTOTSYNTHESIS-------
DLAB.DF.GPP<-DLAB.dfP%>%
  filter(rate.type=="GP")
#anova function
Pmax.modDLABP <- lm(Pmax~treatment, data=DLAB.DF.GPP)
summary(Pmax.modDLABP)
#check for normality, use normality plots

qqnorm(resid(Pmax.modDLABP))
qqline(resid(Pmax.modDLABP))

#check heteroscisity with boxplots

boxplot(resid(Pmax.modDLABP)~DLAB.DF.GPP$treatment)


data.summaryDLABGP<-DLAB.DF.GPP %>%
  group_by(treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(Pmax), se=sd(Pmax)/sqrt(n())) #calculates mean and s.e.
data.summaryDLABGP

PMAX_PDLABGP<-ggplot(data.summaryDLABGP, aes(x=treatment, y=mean, col= treatment, group=factor(treatment))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  labs(x="Depth", y="Maximum Rate of Performance (Pmax) DLAB (GP)", fill="treatment", color = "treatment") 

PMAX_PDLABGP
esquisse:::ggplot_to_ppt("PMAX_PDLABGP")
### Topt OFRA
#anova function
Topt.modDLABGP <- lm(topt_C~treatment, data=DLAB.DF.GPP)
summary(Topt.modDLABGP)
#check for normality, use normality plots

qqnorm(resid(Topt.modDLABGP))
qqline(resid(Topt.modDLABGP))

#check heteroscisity with boxplots

boxplot(resid(Topt.modDLABGP)~DLAB.DF.GPP$treatment)

#high R and low show inconsistent variances, may need to do weighted regression in the future


data.summaryTDLABGP<-DLAB.DF.GPP %>%
  group_by(treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(topt_C), se=sd(topt_C)/sqrt(n())) #calculates mean and s.e.
data.summaryTDLABGP

Topt_PDLABGP<-ggplot(data.summaryTDLABGP, aes(x=treatment, y=mean, col= treatment, group=factor(rate.type))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  labs(x="Depth", y="Thermal Optimum (Topt) DLAB (GP)", fill="treatment", color = "treatment") 

Topt_PDLABGP
esquisse:::ggplot_to_ppt("Topt_PDLABGP")
##### DLAB RESPIRATION------
DLAB.DF.RP<-DLAB.dfP%>%
  filter(rate.type=="R")

#anova function
Pmax.modDLAB <- lm(Pmax~treatment, data=DLAB.DF.RP)
summary(Pmax.modDLAB)
#check for normality, use normality plots

qqnorm(resid(Pmax.modDLAB))
qqline(resid(Pmax.modDLAB))

#check heteroscisity with boxplots

boxplot(resid(Pmax.modDLAB)~DLAB.DF.RP$treatment)


data.summaryDLAB<-DLAB.DF.RP %>%
  group_by(treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(Pmax), se=sd(Pmax)/sqrt(n())) #calculates mean and s.e.
data.summaryDLAB


PMAX_PDLAB<-ggplot(data.summaryDLAB, aes(x=treatment, y=mean, col= treatment, group=factor(treatment))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  labs(x="Depth", y="Maximum Rate of Performance (Pmax) DLAB (R)", fill="treatment", color = "treatment") 

PMAX_PDLAB
esquisse:::ggplot_to_ppt("PMAX_PDLAB")
### Topt DLAB
#anova function
Topt.modDLAB <- lm(topt_C~treatment, data=DLAB.DF.RP)
summary(Topt.modDLAB)
#check for normality, use normality plots

qqnorm(resid(Topt.modDLAB))
qqline(resid(Topt.modDLAB))

#check heteroscisity with boxplots

boxplot(resid(Topt.modDLAB)~DLAB.DF.RP$treatment)

data.summaryTDLAB<-DLAB.DF.RP %>%
  group_by(treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(topt_C), se=sd(topt_C)/sqrt(n())) #calculates mean and s.e.
data.summaryTDLAB

Topt_PDLAB<-ggplot(data.summaryTDLAB, aes(x=treatment, y=mean, col= treatment, group=factor(rate.type))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  labs(x="Depth", y="Thermal Optimum (Topt) DLAB (R)", fill="treatment", color = "treatment") 

Topt_PDLAB
esquisse:::ggplot_to_ppt("Topt_PDLAB")
#-------MCAV PHOTOSYNTHESIS-------
MCAV.DF.GPP<-MCAV.dfP%>%
  filter(rate.type=="GP")

#anova function
Pmax.modMCAVP <- lm(Pmax~treatment, data=MCAV.DF.GPP)
summary(Pmax.modMCAVP)
#check for normality, use normality plots

qqnorm(resid(Pmax.modMCAVP))
qqline(resid(Pmax.modMCAVP))

#check heteroscisity with boxplots

boxplot(resid(Pmax.modMCAVP)~MCAV.DF.GPP$treatment)


data.summaryMCAVGP<-MCAV.DF.GPP %>%
  group_by(treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(Pmax), se=sd(Pmax)/sqrt(n())) #calculates mean and s.e.
data.summaryMCAVGP

PMAX_PMCAVGP<-ggplot(data.summaryMCAVGP, aes(x=treatment, y=mean, col= treatment, group=factor(treatment))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  labs(x="Depth", y="Maximum Rate of Performance (Pmax) MCAV (GP)", fill="treatment", color = "treatment") 

PMAX_PMCAVGP
esquisse:::ggplot_to_ppt("PMAX_PMCAVGP")
### Topt OFRA
#anova function
Topt.modMCAVGP <- lm(topt_C~treatment, data=MCAV.DF.GPP)
summary(Topt.modMCAVGP)
#check for normality, use normality plots

qqnorm(resid(Topt.modMCAVGP))
qqline(resid(Topt.modMCAVGP))

#check heteroscisity with boxplots

boxplot(resid(Topt.modMCAVGP)~MCAV.DF.GPP$treatment)

#high R and low show inconsistent variances, may need to do weighted regression in the future


data.summaryTMCAVGP<-MCAV.DF.GPP %>%
  group_by(treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(topt_C), se=sd(topt_C)/sqrt(n())) #calculates mean and s.e.
data.summaryTMCAVGP

Topt_PMCAVGP<-ggplot(data.summaryTMCAVGP, aes(x=treatment, y=mean, col= treatment, group=factor(rate.type))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  labs(x="Depth", y="Thermal Optimum (Topt) MCAV (GP)", fill="treatment", color = "treatment") 

Topt_PMCAVGP
esquisse:::ggplot_to_ppt("Topt_PMCAVGP")


#### MCAV_) RESPIRATION-----
MCAV.DF.RP<-MCAV.dfP%>%
  filter(rate.type=="R")

#anova function
Pmax.modMCAV <- lm(Pmax~treatment, data=MCAV.DF.RP)
summary(Pmax.modMCAV)
#check for normality, use normality plots

qqnorm(resid(Pmax.modMCAV))
qqline(resid(Pmax.modMCAV))

#check heteroscisity with boxplots

boxplot(resid(Pmax.modMCAV)~MCAV.DF.RP$treatment)


data.summaryMCAV<-MCAV.DF.RP %>%
  group_by(treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(Pmax), se=sd(Pmax)/sqrt(n())) #calculates mean and s.e.
data.summaryMCAV


PMAX_PMCAV<-ggplot(data.summaryMCAV, aes(x=treatment, y=mean, col= treatment, group=factor(treatment))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  labs(x="Depth", y="Maximum Rate of Performance (Pmax) MCAV (R)", fill="treatment", color = "treatment") 

PMAX_PMCAV
esquisse:::ggplot_to_ppt("PMAX_PMCAV")
### Topt MCAV
#anova function
Topt.modMCAV <- lm(topt_C~treatment, data=MCAV.DF.RP)
summary(Topt.modMCAV)
#check for normality, use normality plots

qqnorm(resid(Topt.modMCAV))
qqline(resid(Topt.modMCAV))

#check heteroscisity with boxplots

boxplot(resid(Topt.modMCAV)~MCAV.DF.RP$treatment)

data.summaryTMCAV<-MCAV.DF.RP %>%
  group_by(treatment, rate.type) %>% #tells to group by these two factors
  dplyr::summarise(mean=mean(topt_C), se=sd(topt_C)/sqrt(n())) #calculates mean and s.e.
data.summaryTMCAV

Topt_PMCAV<-ggplot(data.summaryTMCAV, aes(x=treatment, y=mean, col= treatment, group=factor(rate.type))) +
  theme_bw()+
  theme(legend.title=element_text(colour="black", size=14), axis.text.x=element_text(face="bold", color="black", size=16), axis.text.y=element_text(face="bold", color="black", size=13), axis.title.x = element_text(color="black", size=18, face="bold"), axis.title.y = element_text(color="black", size=18, face="bold"),panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  geom_point(position="dodge", size=2) +
  theme(legend.text=element_text(size=rel(1))) + #makes legen elements larger
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  labs(x="Depth", y="Thermal Optimum (Topt) MCAV (R)", fill="treatment", color = "treatment") 

Topt_PMCAV
esquisse:::ggplot_to_ppt("Topt_PMCAV")


###### go to IndividualTPC's file to see TPC's for each genotype across each treatment and rate.type
###### 




