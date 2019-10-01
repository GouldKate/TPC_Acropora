

GPdata<-filtdata%>%
  filter(rate.type=="GP")
  
Redata<-filtdata%>%
  filter(rate.type=="R")

GPT1<-GPdata%>%
  filter(treatment=="T1")
GPT2<-GPdata%>%
  filter(treatment=="T2")
GPT3<-GPdata%>%
  filter(treatment=="T3")
GP0C<-GPdata%>%
  filter(treatment=="C")
ReT1<-Redata%>%
  filter(treatment=="T1")
ReT2<-Redata%>%
  filter(treatment=="T2")
ReT3<-Redata%>%
  filter(treatment=="T3")
Re0C<-Redata%>%
  filter(treatment=="C")






#----run this second to get parameters for each genotype
All.parametersbT<-data.frame()#create an empty df to fill with calculated parameters


mult.fit.curvesbT<-function(Data){
  fit2bT <- nls_multstart(log.rate ~ schoolfield_high(lnc, E, Eh, Th, temp = K, Tc = 27),
                        data = Data,
                        iter = 500,
                        start_lower = c(lnc = -10, E = 0.1, Eh = 0.2, Th = 285),
                        start_upper = c(lnc = 10, E = 2, Eh = 5, Th = 330),
                        supp_errors = 'Y',
                        na.action = na.omit,
                        lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))
  #print(fit2)
  params <- tidy(fit2bT)
  
  params%<>%
    mutate(treatment=Data$treatment[1],
           rate.type=Data$rate.type[1])
  print(params)
  All.parametersbT<<-rbind(params, All.parametersbT)
}


#  
mult.fit.curvesbT(GPT1)
mult.fit.curvesbT(GPT2)
mult.fit.curvesbT(GPT3)
mult.fit.curvesbT(GP0C)
mult.fit.curvesbT(ReT1)
mult.fit.curvesbT(ReT2)
mult.fit.curvesbT(ReT3)
mult.fit.curvesbT(Re0C)

#------T opt from parameters-------
get_topt<-function(E, Th, Eh){
  return((Eh*Th)/(Eh+(8.62e-05*Th*(log((Eh/E))-1))))
}



Topt_bT<-All.parametersbT%>%
  select(rate.type, treatment, estimate, term)%>%
  group_by(rate.type, treatment)%>%
  spread(term,estimate)

write.csv(Topt_bTbT, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/TPC_param_byTreatment.csv")
write.csv(All.parametersbT, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/All.parametersbT.csv")
#bytreat<-read.csv("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/All.parameters_bytreat.csv")

Topt_bT$topt<-get_topt(E=Topt_bT$E, Th=Topt_bT$Th, Eh=Topt_bT$Eh)

Topt_bT$topt_C<-Topt_bT$topt-273.15

write.csv(Topt_bT, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/R_files/R_Output/Topt_estimatesbT.csv")


#-------linear Models does treatment or species affect Topt----
library(car)
library(lme4)
library(stats)
library(Rmisc)
library(Hmisc)

# plot distribution of Topt ##### ------ Parameter graphs--- 
Toptmeans<-ggplot(Topt_bT, aes(topt_C)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~rate.type)+
  xlab('Optimum Temperature (ºC)') +
  ggtitle('Distribution of optimum temperatures')

esquisse:::ggplot_to_ppt("Toptmeans")

GPtbT<-ggplot(subset(Topt_bT,rate.type=="GP"),
              aes(x=treatment, topt_C, colour=treatment, shape=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 5)+
  geom_point(size=4)+
  ylab("Topt (ºC)")+
  labs(title="Mean Topt across Treatment (GP)")+
  coord_flip()
GPtbT
esquisse:::ggplot_to_ppt("GPtbT")






