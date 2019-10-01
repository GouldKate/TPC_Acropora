topt_Cever<-ggplot(subset(Pmax1,rate.type=="GP"), aes(x=replicas, topt_C, colour=treatment, shape= genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Topt (ºC)")+
  xlab("Replicas")+
  labs(title="Mean Topt (ºC) across Replicas (GP)")+
    coord_flip()
esquisse:::ggplot_to_ppt("topt_Cever")



topt_RCever<-ggplot(subset(Pmax1,rate.type=="R"), aes(x=replicas, topt_C, colour=treatment, shape= genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=4)+
  ylab("Topt (ºC)")+
  xlab("Replicas")+
  labs(title="Mean Topt (ºC) across Replicas (R)")+
  coord_flip()
esquisse:::ggplot_to_ppt("topt_RCever")


mydataP$log.rate <- mydata$umol.cm2.hr + )  #logging and adding 0.3(-2 was smallest in data set) because a log of zero does not exist
mydataR$log.rate <- log(mydata$umol.cm2.hr + 3)  #logging and adding 0.3(-2 was smallest in data set) because a log of zero does not exist
