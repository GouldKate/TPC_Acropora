#genotype
TPCallg<-ggplot() +
  geom_ribbon(data=subset(All.CIg, lwr_CI>0),
              aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedg, aes(x=(K - 273.15), y=log.rate, shape=genotype)) +
  geom_line(data=All.fittedg, aes(x=(K - 273.15), y=fitted, colour=group)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~treatment)
esquisse:::ggplot_to_ppt("TPCallg")

# grouped by treatment and rate.type- facetwrapped by genotype- sad

TPCall<-ggplot() +
  geom_ribbon(data=subset(All.CIT, lwr_CI>0),
              aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedT, aes(x=(K - 273.15), y=log.rate, shape=treatment)) +
  theme(legend.position = "none")+
  geom_line(data=All.fittedT, aes(x=(K - 273.15), y=fitted, colour=group)) +
  theme(legend.position = "none")+
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~genotype)
esquisse:::ggplot_to_ppt("TPCall")





TPCallg<-ggplot() +
  geom_ribbon(data=subset(All.CIr, lwr_CI>0),
              aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedr, aes(x=(K - 273.15), y=log.rate, shape=treatment)) +
  geom_line(data=All.fittedr, aes(x=(K - 273.15), y=fitted, color=group)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~genotype)
esquisse:::ggplot_to_ppt("TPCallg")


#individual.ID
TPCallID<-ggplot() +
  geom_ribbon(data=subset(All.CIr, lwr_CI>0),
              aes(x=K - 273.15, ymin = lwr_Cr, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedr, aes(x=(K - 273.15), y=log.rate, shape=rate.type)) +
  geom_line(data=All.fittedr, aes(x=(K - 273.15), y=fitted, colour=group)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  theme_bw(base_size = 12, base_family = 'Helvetica') +
  facet_wrap(~genotype)
esquisse:::ggplot_to_ppt("TPCallID")

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


TPCTm<-ggplot() +
  geom_ribbon(data=subset(All.CIi, lwr_CI>0),aes(x=K - 273.15, ymin = lwr_CI, ymax = upr_CI, group=group),fill = 'grey', alpha = .4) +
  geom_point(data=All.fittedID, aes(x=(K - 273.15), y=log.rate, shape=genotype))+ 
  geom_smooth(data=All.fittedID, aes(x=(K - 273.15), y=fitted, colour=group,se=F)) +
  ylab('log Metabolic rate') +
  xlab('Temperature (ºC)') +
  facet_wrap(~genotype)+
  theme_classic()+
  labs(title="Thermal Performance Curves Across treatments")
esquisse:::ggplot_to_ppt("TPCTm")