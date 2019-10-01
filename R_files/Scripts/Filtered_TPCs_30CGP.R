#filtering data more to compare TPC's
#
#Removing GP&R- Temp 30- M5
#removing GP- Temp 29- K2 & M6
#
#Remove ALL treatments 30C
#
##############################################-----

filtdata2 <- mydata %>%
  filter(!(temp.Cat == "30"& rate.type == "GP"))

#step one subset data: first make subsets for all of your treatment groups 

# by genotype, K2, M5, M6# 

K2.DF<-filtdata2%>%
  filter(genotype=="K2")
M5.DF<-filtdata2%>%
  filter(genotype=="M5")
M6.DF<-filtdata2%>%
  filter(genotype=="M6")

#by Rate.Type GP / R#

#---- K2------
K2.DF.GP<-K2.DF%>%
  filter(rate.type=="GP")
K2.DF.R<-K2.DF%>%
  filter(rate.type=="R")

#----M5------
M5.DF.GP<-M5.DF%>%
  filter(rate.type=="GP")
M5.DF.R<-M5.DF%>%
  filter(rate.type=="R")

#------M6------
M6.DF.GP<-M6.DF%>%
  filter(rate.type=="GP")
M6.DF.R<-M6.DF%>%
  filter(rate.type=="R")


#-finally by treatment T1, T2, C & T3

K2.DF.GP.T1<-K2.DF.GP%>%filter(treatment=="T1")
M5.DF.GP.T1<-M5.DF.GP%>%filter(treatment=="T1")
M6.DF.GP.T1<-M6.DF.GP%>%filter(treatment=="T1")
K2.DF.GP.T2<-K2.DF.GP%>%filter(treatment=="T2")
M5.DF.GP.T2<-M5.DF.GP%>%filter(treatment=="T2")
M6.DF.GP.T2<-M6.DF.GP%>%filter(treatment=="T2")
K2.DF.GP.C<-K2.DF.GP%>%filter(treatment=="C")
M5.DF.GP.C<-M5.DF.GP%>%filter(treatment=="C")
M6.DF.GP.C<-M6.DF.GP%>%filter(treatment=="C")
K2.DF.GP.T3<-K2.DF.GP%>%filter(treatment=="T3")
M5.DF.GP.T3<-M5.DF.GP%>%filter(treatment=="T3")
M6.DF.GP.T3<-M6.DF.GP%>%filter(treatment=="T3")


K2.DF.R.T1<-K2.DF.R%>%filter(treatment=="T1")
M5.DF.R.T1<-M5.DF.R%>%filter(treatment=="T1")
M6.DF.R.T1<-M6.DF.R%>%filter(treatment=="T1")
K2.DF.R.T2<-K2.DF.R%>%filter(treatment=="T2")
M5.DF.R.T2<-M5.DF.R%>%filter(treatment=="T2")
M6.DF.R.T2<-M6.DF.R%>%filter(treatment=="T2")
K2.DF.R.C<-K2.DF.R%>%filter(treatment=="C")
M5.DF.R.C<-M5.DF.R%>%filter(treatment=="C")
M6.DF.R.C<-M6.DF.R%>%filter(treatment=="C")
K2.DF.R.T3<-K2.DF.R%>%filter(treatment=="T3")
M5.DF.R.T3<-M5.DF.R%>%filter(treatment=="T3")
M6.DF.R.T3<-M6.DF.R%>%filter(treatment=="T3")



#-------calculate parameter estimates per genotype-------
#Step one: make empty DF to fill with fitted values and run function
##"""it is very important to empty the dataframes before running the function multiple times"""

All.fitted<-data.frame()#create an empty df to fill with fitted values
All.CI<-data.frame()#create an empty df. to fill with confidence intervals

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
  
  All.fitted<<-rbind(All.fitted,Data)
  
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
    unnest(fit %>% map(augment)) %>%
    ungroup()
  
  new_preds <- Data %>%
    do(., data.frame(K = seq(min(.$K), max(.$K), length.out = 250), stringsAsFactors = FALSE))
  
  preds <- augment(fit2, newdata = new_preds)
  
  df.genotype<-(Data$genotype[1])
  df.rate.type<-(Data$rate.type[1])
  df.treatment<-(Data$treatment[1])
  
  
  preds <- fit_boots %>%
    unnest(fit %>% map(augment, newdata = new_preds)) %>%
    # group by each value of K and get quantiles
    group_by(., K) %>%
    summarise(lwr_CI = quantile(.fitted, 0.025),
              upr_CI = quantile(.fitted, 0.975)) %>%
    ungroup() %>%
    merge(., preds, by = 'K')%>%
    mutate(genotype=factor(df.genotype),
           treatment=factor(df.treatment),
           rate.type=factor(df.rate.type))
  All.CI<<-rbind(All.CI, preds)
  
}

#"""it is very important to empty the dataframes before running the function multiple times"""


#-----put all the subsetted data through the function--------
#First run all of these and get All.fitted. and All. CI


mult.fit.curves(K2.df.GP.T1)
mult.fit.curves(M5.DF.GP.T1)
mult.fit.curves(M6.DF.GP.T1)
mult.fit.curves(K2.DF.GP.T2)
mult.fit.curves(M5.DF.GP.T2)
mult.fit.curves(M6.DF.GP.T2)
mult.fit.curves(K2.DF.GP.C)
mult.fit.curves(M5.DF.GP.C)
mult.fit.curves(M6.DF.GP.C)
mult.fit.curves(K2.DF.GP.T3)
mult.fit.curves(M5.DF.GP.T3)
mult.fit.curves(M6.DF.GP.T3)

mult.fit.curves(K2.DF.R.T1)
mult.fit.curves(M5.DF.R.T1)
mult.fit.curves(M6.DF.R.T1)
mult.fit.curves(K2.DF.R.T2)
mult.fit.curves(M5.DF.R.T2)
mult.fit.curves(M6.DF.R.T2)
mult.fit.curves(K2.DF.R.C)
mult.fit.curves(M5.DF.R.C)
mult.fit.curves(M6.DF.R.C)
mult.fit.curves(K2.DF.R.T3)
mult.fit.curves(M5.DF.R.T3)
mult.fit.curves(M6.DF.R.T3)




#----Graph the outputs--------
#graph with your wobbly predictions
All.fitted%<>%
  mutate(group=factor(paste(rate.type,genotype)))

All.CI%<>%
  mutate(group=factor(paste(rate.type,genotype)))


#-----Save plot in PP-------##
library(esquisse)
library(rvg)

TPCall<-ggplot() +
  geom_ribbon(data=subset(All.CI, lwr_CI>0),
              aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fitted, aes(x=(K - 273.15), y=log.rate, shape=genotype)) +
  geom_line(data=All.fitted, aes(x=(K - 273.15), y=fitted, colour=group)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~treatment)
esquisse:::ggplot_to_ppt("TPCall")



#----Graph the outputs--------
#graph with your wobbly predictions
All.fittedg2<-All.fitted%<>%
  mutate(group=factor(paste(treatment,rate.type)))

All.CIg2<-All.CI%<>%
  mutate(group=factor(paste(treatment,rate.type)))
####
TPCallg2<-ggplot() +
  geom_ribbon(data=subset(All.CIg2, lwr_CI>0),
              aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedg2, aes(x=(K - 273.15), y=log.rate, shape=treatment)) +
  geom_line(data=All.fittedg2, aes(x=(K - 273.15), y=fitted, colour=group)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~genotype)
esquisse:::ggplot_to_ppt("TPCallg2")




#graph using a smoothing function
STPCall2<-ggplot() +
  geom_ribbon(data=subset(All.CI, lwr_CI>0),aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fitted, aes(x=(K - 273.15), y=log.rate, shape=genotype)) +
  geom_smooth(data=All.fitted, aes(x=(K - 273.15), y=fitted, colour=group, se=F)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~treatment)

esquisse:::ggplot_to_ppt("STPCall2")

#graph using a smoothing function
STPCallg2<-ggplot() +
  geom_ribbon(data=subset(All.CIg2, lwr_CI>0),aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedg2, aes(x=(K - 273.15), y=log.rate, shape=treatment)) +
  geom_smooth(data=All.fittedg2, aes(x=(K - 273.15), y=fitted, colour=group,se=F)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~genotype)

esquisse:::ggplot_to_ppt("STPCallg2")


#graph using a smoothing function and fancy up this graphical representation of biological processes
TPC2<-ggplot() +
  geom_ribbon(data=subset(All.CI, lwr_CI>0),aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fitted, aes(x=(K - 273.15), y=log.rate, shape=genotype))+ 
  geom_smooth(data=All.fitted, aes(x=(K - 273.15), y=fitted, colour=group,se=F)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  facet_wrap(~treatment)+
  theme_classic()+
  labs(title="Thermal Performance Curves Across Treatments")
esquisse:::ggplot_to_ppt("TPC2")

TPCg2<-ggplot() +
  geom_ribbon(data=subset(All.CIg2, lwr_CI>0),aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedg2, aes(x=(K - 273.15), y=log.rate, shape=treatment))+ 
  geom_smooth(data=All.fittedg2, aes(x=(K - 273.15), y=fitted, colour=group,se=F)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  facet_wrap(~genotype)+
  theme_classic()+
  labs(title="Thermal Performance Curves Across Treatments")
esquisse:::ggplot_to_ppt("TPCg2")




##### more subsets
##### 
##### 
##### 
  
K201a<-(subset(K2.df.GP.T1, individual.ID== "K201a"))
M501a<-(subset(M5.df.GP.T1, individual.ID== "M501a"))
M601a<-(subset(M6.df.GP.T1, individual.ID== "M601a"))
K202a<-(subset(K2.df.GP.T2, individual.ID== "K202a"))
M502a<-(subset(M5.df.GP.T2, individual.ID== "M502a"))
M602a<-(subset(M6.df.GP.T2, individual.ID== "M602a"))
K2Ca<-(subset(K2.df.GP.C, individual.ID== "K2Ca"))
M5Ca<-(subset(M5.df.GP.C, individual.ID== "M5Ca"))
M6Ca<-(subset(M6.df.GP.C, individual.ID== "M6Ca"))
K203a<-(subset(K2.df.GP.T3, individual.ID== "K203a"))
M503a<-(subset(M5.df.GP.T3, individual.ID== "M503a"))
M603a<-(subset(M6.df.GP.T3, individual.ID== "M603a"))
K201a<-(subset(K2.df.R.T1, individual.ID==  "K201a"))
M501a<-(subset(M5.df.R.T1, individual.ID==  "M501a"))
M601a<-(subset(M6.df.R.T1, individual.ID==  "M601a"))
K202a<-(subset(K2.df.R.T2, individual.ID==  "K202a"))
M502a<-(subset(M5.df.R.T2, individual.ID==  "M502a"))
M602a<-(subset(M6.df.R.T2, individual.ID==  "M602a"))
K20Ca<-(subset(K2.df.R.C, individual.ID==  "K20Ca"))
M50Ca<-(subset(M5.df.R.C, individual.ID==  "M50Ca"))
M60Ca<-(subset(M6.df.R.C, individual.ID==  "M60Ca"))
K203a<-(subset(K2.df.R.T3, individual.ID==  "K203a"))
M503a<-(subset(M5.df.R.T3, individual.ID==  "M503a"))
M603a<-(subset(M6.df.R.T3, individual.ID==  "M603a"))
K201b<-(subset(K2.df.GP.T1, individual.ID== "K201b"))
M501b<-(subset(M5.df.GP.T1, individual.ID== "M501b"))
M601b<-(subset(M6.df.GP.T1, individual.ID== "M601b"))
K202b<-(subset(K2.df.GP.T2, individual.ID== "K202b"))
M502b<-(subset(M5.df.GP.T2, individual.ID== "M502b"))
M602b<-(subset(M6.df.GP.T2, individual.ID== "M602b"))
K20Cb<-(subset(K2.df.GP.C, individual.ID== "K20Cb"))
M50Cb<-(subset(M5.df.GP.C, individual.ID== "M50Cb"))
M60Cb<-(subset(M6.df.GP.C, individual.ID== "M60Cb"))
K203b<-(subset(K2.df.GP.T3, individual.ID== "K203b"))
M503b<-(subset(M5.df.GP.T3, individual.ID== "M503b"))
M603b<-(subset(M6.df.GP.T3, individual.ID== "M603b"))
K201b<-(subset(K2.df.R.T1, individual.ID== "K201b"))
M501b<-(subset(M5.df.R.T1, individual.ID== "M501b"))
M601b<-(subset(M6.df.R.T1, individual.ID== "M601b"))
K202b<-(subset(K2.df.R.T2, individual.ID== "K202b"))
M502b<-(subset(M5.df.R.T2, individual.ID== "M502b"))
M602b<-(subset(M6.df.R.T2, individual.ID== "M602b"))
K20Cb<-(subset(K2.df.R.C, individual.ID== "K20Cb"))
M50Cb<-(subset(M5.df.R.C, individual.ID== "M50Cb"))
M60Cb<-(subset(M6.df.R.C, individual.ID== "M60Cb"))
K203b<-(subset(K2.df.R.T3, individual.ID==  "K203b"))
M503b<-(subset(M5.df.R.T3, individual.ID==  "M503b"))
M603b<-(subset(M6.df.R.T3, individual.ID==  "M603b"))
K201c<-(subset(K2.df.GP.T1, individual.ID== "K201c"))
M501c<-(subset(M5.df.GP.T1, individual.ID== "M501c"))
M601c<-(subset(M6.df.GP.T1, individual.ID== "M601c"))
K202c<-(subset(K2.df.GP.T2, individual.ID== "K202c"))
M502c<-(subset(M5.df.GP.T2, individual.ID== "M502c"))
M602c<-(subset(M6.df.GP.T2, individual.ID== "M602c"))
K20Cc<-(subset(K2.df.GP.C, individual.ID== "K20Cc"))
M50Cc<-(subset(M5.df.GP.C, individual.ID== "M50Cc"))
M60Cc<-(subset(M6.df.GP.C, individual.ID== "M60Cc"))
K203c<-(subset(K2.df.GP.T3, individual.ID== "K203c"))
M503c<-(subset(M5.df.GP.T3, individual.ID== "M503c"))
M603c<-(subset(M6.df.GP.T3, individual.ID== "M603c"))
K201c<-(subset(K2.df.R.T1, individual.ID==  "K201c"))
M501c<-(subset(M5.df.R.T1, individual.ID==  "M501c"))
M601c<-(subset(M6.df.R.T1, individual.ID==  "M601c"))
K202c<-(subset(K2.df.R.T2, individual.ID==  "K202c"))
M502c<-(subset(M5.df.R.T2, individual.ID==  "M502c"))
M602c<-(subset(M6.df.R.T2, individual.ID==  "M602c"))
K20Cc<-(subset(K2.df.R.C, individual.ID== "K20Cc"))
M50Cc<-(subset(M5.df.R.C, individual.ID== "M50Cc"))
M60Cc<-(subset(M6.df.R.C, individual.ID== "M60Cc"))
K203c<-(subset(K2.df.R.T3, individual.ID==  "K203c"))
M503c<-(subset(M5.df.R.T3, individual.ID==  "M503c"))
M603c<-(subset(M6.df.R.T3, individual.ID==  "M603c"))
K201d<-(subset(K2.df.GP.T1, individual.ID== "K201d"))
M501d<-(subset(M5.df.GP.T1, individual.ID== "M501d"))
M601d<-(subset(M6.df.GP.T1, individual.ID== "M601d"))
K202d<-(subset(K2.df.GP.T2, individual.ID== "K202d"))
M502d<-(subset(M5.df.GP.T2, individual.ID== "M502d"))
M602d<-(subset(M6.df.GP.T2, individual.ID== "M602d"))
K20Cd<-(subset(K2.df.GP.C, individual.ID== "K20Cd"))
M50Cd<-(subset(M5.df.GP.C, individual.ID== "M50Cd"))
M60Cd<-(subset(M6.df.GP.C, individual.ID== "M60Cd"))
K203d<-(subset(K2.df.GP.T3, individual.ID== "K203d"))
M503d<-(subset(M5.df.GP.T3, individual.ID== "M503d"))
M603d<-(subset(M6.df.GP.T3, individual.ID== "M603d"))
K201d<-(subset(K2.df.R.T1, individual.ID==  "K201d"))
M501d<-(subset(M5.df.R.T1, individual.ID==  "M501d"))
M601d<-(subset(M6.df.R.T1, individual.ID==  "M601d"))
K202d<-(subset(K2.df.R.T2, individual.ID==  "K202d"))
M502d<-(subset(M5.df.R.T2, individual.ID==  "M502d"))
M602d<-(subset(M6.df.R.T2, individual.ID==  "M602d"))
K20Cd<-(subset(K2.df.R.C, individual.ID== "K20Cd"))
M50Cd<-(subset(M5.df.R.C, individual.ID== "M50Cd"))
M60Cd<-(subset(M6.df.R.C, individual.ID== "M60Cd"))
K203d<-(subset(K2.df.R.T3, individual.ID==  "K203d"))
M503d<-(subset(M5.df.R.T3, individual.ID==  "M503d"))
M603d<-(subset(M6.df.R.T3, individual.ID==  "M603d"))



#### Extra code from Parameter means
#### ####
####
#-------Visualize and fit to linear model----
#library(car)
#library(lme4)
#library(stats)
##library(Rmisc)
#library(Hmisc)


#visualize groups GP
#GPbT<-ggplot(subset(topt_data,rate.type=="GP"), aes(x=Treatment, topt_C, colour=species, shape=Treatment))+
#  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
#  geom_point(size=4)+
#  ylab("Topt (ºC)")+
#  labs(title="Photosynthesis Topt by Treatment")+
#  coord_flip()
#
#esquisse:::ggplot_to_ppt("GPbT")
#
#GPbS<-ggplot(subset(topt_data,rate.type=="GP"), aes(x=species, topt_C, colour=species, shape=Treatment))+
#  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
#  geom_point(size=4)+
#  ylab("Topt (ºC)")+
#  labs(title="Photosynthesis Topt by species")+
#  coord_flip()
#
#esquisse:::ggplot_to_ppt("GPbS")
#
#ggplot(subset(topt_data,rate.type=="GP"), aes(x=Treatment, topt_C, colour=Treatment))+
#  geom_boxplot()+
#  geom_point()
#
#ggplot(subset(topt_data,rate.type=="GP"), aes(x=species, topt_C, colour=species))+
#  geom_boxplot()+
#  geom_point()
#
#
#
##linear model
#Topt.modGP<-Anova(lm(topt_C~Treatment+species, data=subset(topt_data,rate.type=="GP")))
#summary(Topt.modGP)
#
##check for normality, use normality plots
#Topt.mod2 <- lm(topt_C~Treatment+species, data=subset(topt_data,rate.type=="GP"))
#
#
##do I need to do this?
##qqnorm(resid(All.fitted$fitted))
##qqline(resid(All.fitted$fitted))
#
#
##check heteroscisity with boxplots
#
##boxplot(resid(Topt.mod)~Topt_data$treatment*Topt_data$rate.type)
#
############
#Anova(lm(topt_C~Treatment+species, data=subset(topt_data,rate.type=="GP")))
#
#m1<-aov(topt_C~Treatment+species, data=subset(topt_data,rate.type=="GP"))
#TukeyHSD(m1, "Treatment")
#TukeyHSD(m1, "species")
#
##visualize groups R
#ggplot(subset(topt_data,rate.type=="R"), aes(x=Treatment, topt_C, colour=Treatment))+
#  geom_boxplot()+
#  geom_point()
#
#ggplot(subset(topt_data,rate.type=="R"), aes(x=species, topt_C, colour=species))+
#  geom_boxplot()+
#  geom_point()
#
#
#RbTB<-ggplot(subset(topt_data,rate.type=="R"), aes(x=Treatment, topt_C, colour=Treatment))+
#  geom_boxplot()+
#  geom_point()+
#  ylab("Topt (ºC)")+
#  labs(title="Repiration Topt by Treatment")
#esquisse:::ggplot_to_ppt("RbTB")
#
#RbT<-ggplot(subset(topt_data,rate.type=="R"), aes(x=Treatment, topt_C, colour=species, shape=Treatment))+
#  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
#  geom_point(size=4)+
#  ylab("Topt (ºC)")+
#  labs(title="Photosynthesis Topt by species")+
#  coord_flip()
#
#esquisse:::ggplot_to_ppt("GPbS")
#
#
#
#ggplot(subset(topt_data,rate.type=="R"), aes(x=species, topt_C, colour=species, shape=Treatment))+
#  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
#  geom_point(size=4)
#
##linear model
#Anova(lm(topt_C~Treatment+species, data=subset(topt_data,rate.type=="R")))
#
#m2<-aov(topt_C~Treatment+species, data=subset(topt_data,rate.type=="R"))
#TukeyHSD(m2, "Treatment")
#TukeyHSD(m2, "species")
#
##-------linear models does treatment or species lnc----
##visualize groups GP
#ggplot(subset(topt_data,rate.type=="GP"), aes(x=Treatment, lnc, colour=Treatment))+
#  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
#  geom_point(size=4)
#
#ggplot(subset(topt_data,rate.type=="GP"), aes(x=species, lnc, colour=species))+
#  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
#  geom_point(size=4)
#
#
#
##linear model
#Anova(lm(lnc~Treatment+species, data=subset(topt_data,rate.type=="GP")))
#
#m3<-aov(lnc~Treatment+species, data=subset(topt_data,rate.type=="GP"))
#TukeyHSD(m3, "Treatment")
#TukeyHSD(m3, "species")
#
##visualize groups R
#ggplot(subset(topt_data,rate.type=="R"), aes(x=Treatment, lnc, colour=Treatment))+
#  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
#  geom_point(size=4)
#
#
#ggplot(subset(topt_data,rate.type=="R"), aes(x=species, lnc, colour=species))+
#  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
#  geom_point(size=4)
#
#
#
#
##linear model
#Anova(lm(lnc~Treatment+species, data=subset(topt_data,rate.type=="R")))
#
#m4<-aov(lnc~Treatment+species, data=subset(topt_data,rate.type=="R"))
#TukeyHSD(m4, "Treatment")
#TukeyHSD(m4, "species")
############3right here#########3-------------