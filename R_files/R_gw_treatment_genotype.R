
###############################
install.packages("ggplot2")
library(ggplot2)
##################################
gw_acropora_04252018
head(GW_acropora_04252018)

gw<-(GW_acropora_04252018)

head(gw)

ggplot(gw)+geom_boxplot(aes(x=factor(treatment),y=chnge_bt))


ggplot(gw)+geom_boxplot(aes(x=factor(genotype),y=chnge_bt))

ggplot(gw)+geom_boxplot(mapping=aes(x=genotype,y=perc_wt_chng))

ggplot(gw)+geom_point(aes(x=factor(coralid),y=perc_wt_chng))

ggplot(gw)+geom_point(aes(x=genotype,y=perc_wt_chng))

ggplot(gw)+geom_boxplot(aes(x=genotype,y=perc_wt_chng))+ 
  labs(x="Genotype",y="Percent Weight Change")
##########################################
meanperchan <- mean(gw$perc_wt_chng)
meanperchan
21.32177
sd(gw$perc_wt_chng)
10.45504
seperchan <- sd(gw$perc_wt_chng)/length(gw$perc_wt_chng)
seperchan
0.217813

dnorm(0)
0.3989423

pnorm(0)
0.5

qnorm(.975)
1.959964
pnorm(qnorm(.975))
0.975

serich*qnorm(.975)
0.4269064
qnorm(.975, sd=serich)
0.4269064

c(meanperchan+seperchan*qnorm(.975),meanperchan-seperchan*qnorm(.975))
21.74868 20.89486

c(meanperchan+seperchan*qt(.975,length(gw$perc_wt_chng)-1),
  meanperchan-seperchan*qt(.975,length(gw$perc_wt_chng-1))
  
  ## didn't work
  ######################3
  

  
  
  meantreatment <- tapply(gw$perc_wt_chng,gw$treatment,mean)
  meantreatment
  
  1        2        3        C 
  22.41404 20.85684 19.02169 22.99452 
  1        2        3        4 
  22.41404 20.85684 19.02169 22.99452
  meangeno <- tapply(gw$perc_wt_chng,gw$genotype,mean)
  meangeno
  
  K2       M5       M6 
  29.69014 16.76906 17.50612
  2        5        6 
  29.69014 16.76906 17.50612
  
  
  
  vartreat <- tapply(gw$perc_wt_chng,gw$treatment,var)
  vartreat
  1         2         3         C 
  87.80492 168.52800  28.77259 171.57621 
  sdn <- tapply(gw$perc_wt_chng,gw$treatment,length)
  sdn
  1  2  3  C 
  12 12 12 12 
  
  vargeno <- tapply(gw$perc_wt_chng,gw$genotype,var)
  vargeno
  
  K2        M5        M6 
  173.31448  43.91084  12.93573
  
  sdn <- tapply(gw$perc_wt_chng,gw$genotype,length)
  sdn
  
  K2 M5 M6 
  16 16 16 
  
  ############################
bw2$treatment.f<- factor(bw2$treatment)
  ### warnings
  
  ggplot(bw2,aes(x=treatment,y=chng_bw, color=treatment))+ geom_point()
  
  
  
  ##############################
  
  t.test(rikz$richness[rikz$week==2],
         rikz$richness[rikz$week==3],var.equal=T)
  
  
  t.test(gw$perc_wt_chng[gw$treatment==1],
         gw$perc_wt_chng[gw$treatment==2],var.equal=T)
  
  Two Sample t-test
  
  data:  gw$perc_wt_chng[gw$treatment == 1] and gw$perc_wt_chng[gw$treatment == 2]
  t = 0.33693, df = 22, p-value = 0.7394
  alternative hypothesis: true difference in means is not equal to 0
  95 percent confidence interval:
    -8.027831 11.142241
  sample estimates:
    mean of x mean of y 
  22.41404  20.85684 
  
  t.test(gw$perc_wt_chng[gw$treatment==1],
         gw$perc_wt_chng[gw$treatment==3],var.equal=T)
  
  Two Sample t-test
  
  data:  gw$perc_wt_chng[gw$treatment == 1] and gw$perc_wt_chng[gw$treatment == 3]
  t = 1.0884, df = 22, p-value = 0.2882
  alternative hypothesis: true difference in means is not equal to 0
  95 percent confidence interval:
    -3.071611  9.856318
  sample estimates:
    mean of x mean of y 
  22.41404  19.02169 
  
  t.test(gw$perc_wt_chng[gw$treatment==1],
         gw$perc_wt_chng[gw$treatment==4],var.equal=T)
  
  Two Sample t-test
  
  data:  gw$perc_wt_chng[gw$treatment == 1] and gw$perc_wt_chng[gw$treatment == 4]
  t = -0.12486, df = 22, p-value = 0.9018
  alternative hypothesis: true difference in means is not equal to 0
  95 percent confidence interval:
    -10.222336   9.061381
  sample estimates:
    mean of x mean of y 
  22.41404  22.99452 
  
  
  t.test(gw$perc_wt_chng[gw$genotype==2],
         gw$perc_wt_chng[gw$genotype==5],var.equal=T)
  
  Two Sample t-test
  
  data:  gw$perc_wt_chng[gw$genotype == 2] and gw$perc_wt_chng[gw$genotype == 5]
  t = 3.5067, df = 30, p-value = 0.001451
  alternative hypothesis: true difference in means is not equal to 0
  95 percent confidence interval:
    5.396038 20.446122
  sample estimates:
    mean of x mean of y 
  29.69014  16.76906 
  
  t.test(gw$perc_wt_chng[gw$genotype==2],
         gw$perc_wt_chng[gw$genotype==6],var.equal=T)
  
  Two Sample t-test
  
  data:  gw$perc_wt_chng[gw$genotype == 2] and gw$perc_wt_chng[gw$genotype == 6]
  t = 3.5711, df = 30, p-value = 0.001222
  alternative hypothesis: true difference in means is not equal to 0
  95 percent confidence interval:
    5.216117 19.151922
  sample estimates:
    mean of x mean of y 
  29.69014  17.50612 
  
  
  
  t.test(gw$perc_wt_chng[gw$genotype==5],
         gw$perc_wt_chng[gw$genotype==6],var.equal=T)
  
  Two Sample t-test
  
  data:  gw$perc_wt_chng[gw$genotype == 5] and gw$perc_wt_chng[gw$genotype == 6]
  t = -0.39103, df = 30, p-value = 0.6985
  alternative hypothesis: true difference in means is not equal to 0
  95 percent confidence interval:
    -4.586574  3.112453
  sample estimates:
    mean of x mean of y 
  16.76906  17.50612 
  
  
aov(gw$perc_wt_chng[gw$genotype==5],
    gw$perc_wt_chng[gw$genotype==6],var.equal=T)
  
  ##################################

genmod2<- lm(chng_bw~factor(genotype), data=gw) 
  
summary(genmod2)

Call:
  lm(formula = chng_bw ~ factor(genotype), data = bw2)

Residuals:
  Min      1Q  Median      3Q     Max 
-2.2638 -1.0125 -0.4463  0.3791  8.4700 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)    
(Intercept)         1.6300     0.4550   3.583 0.000831 ***
  factor(genotype)5  -0.2663     0.6434  -0.414 0.680979    
factor(genotype)6   0.2787     0.6434   0.433 0.666912    
---
  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 1.82 on 45 degrees of freedom
Multiple R-squared:  0.0157,	Adjusted R-squared:  -0.02805 
F-statistic: 0.3588 on 2 and 45 DF,  p-value: 0.7005

angen<-anova(genmod)  

angen

treatmod<- lm(chng_bw~factor(treatment), data=bw2)  

Analysis of Variance Table

Response: chng_bw
Df  Sum Sq Mean Sq F value Pr(>F)
factor(genotype)  2   2.377  1.1883  0.3588 0.7005
Residuals        45 149.033  3.3118

  summary(treatmod)
  
  Call:
    lm(formula = chng_bw ~ factor(treatment), data = bw2)
  
  Residuals:
    Min      1Q  Median      3Q     Max 
  -2.9117 -0.9992 -0.4767  0.4996  8.0883 
  
  Coefficients:
    Estimate Std. Error t value Pr(>|t|)   
  (Intercept)          1.7642     0.5274   3.345  0.00169 **
  factor(treatment)2   0.2475     0.7459   0.332  0.74160   
  factor(treatment)3  -0.1725     0.7459  -0.231  0.81818   
  factor(treatment)4  -0.5950     0.7459  -0.798  0.42933   
  ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
  
  Residual standard error: 1.827 on 44 degrees of freedom
  Multiple R-squared:  0.02991,	Adjusted R-squared:  -0.03623 
  F-statistic: 0.4523 on 3 and 44 DF,  p-value: 0.717
  
  
  
antreat<-anova(treatmod)  
  
antreat

Analysis of Variance Table

Response: chng_bw
                  Df  Sum Sq Mean Sq F value Pr(>F)
factor(treatment)  3   4.529  1.5097  0.4523  0.717
Residuals         44 146.880  3.3382 




treano<-anova(treatmod)


  
  
summary(treano)  
head(bw2)


  
treat.aov<- aov(bw2$chng_bw ~ bw2$treatment)

treat.aov
summary(treat.aov)


anova(bw2$chng_bw ~ bw2$treatment==)
  
  wilcox.test(rikz$richness[rikz$week==2],rikz$richness[rikz$week==3])
  ##########################
  
  
  fit<- aov()