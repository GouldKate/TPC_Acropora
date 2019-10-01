
###############################
install.packages("ggplot2")
library(ggplot2)
##################################
bw2_acropora_04252018
head(bw2_acropora_04252018)

bw2<-(bw2_acropora_04252018)

head(bw2)

ggplot(bw2)+geom_boxplot(aes(x=factor(treatment),y=perc_wt_chng))

ggplot(bw2)+geom_boxplot(aes(x=treatment,y=perc_wt_chng))+ 
  labs(x="Treatment",y="Percent Weight Change")

ggplot(bw2)+geom_boxplot(aes(x=factor(genotype),y=perc_wt_chng))

ggplot(bw2)+geom_boxplot(mapping=aes(x=genotype,y=perc_wt_chng))

ggplot(bw2)+geom_point(aes(x=factor(coralid),y=perc_wt_chng))

ggplot(bw2)+geom_point(aes(x=genotype,y=perc_wt_chng))

ggplot(bw2)+geom_boxplot(aes(x=genotype,y=perc_wt_chng))+ 
  labs(x="Genotype",y="Percent Weight Change")
##########################################
meanperchan <- mean(bw2$perc_wt_chng)
meanperchan
21.32177
sd(bw2$perc_wt_chng)
10.45504
seperchan <- sd(bw2$perc_wt_chng)/length(bw2$perc_wt_chng)
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

c(meanperchan+seperchan*qt(.975,length(bw2$perc_wt_chng)-1),
  meanperchan-seperchan*qt(.975,length(bw2$perc_wt_chng-1))
  
  ## didn't work
  
  
  meantreatment <- tapply(bw2$perc_wt_chng,bw2$treatment,mean)
  meantreatment
  
  1        2        3        C 
  22.41404 20.85684 19.02169 22.99452 
  1        2        3        4 
  22.41404 20.85684 19.02169 22.99452
  meangeno <- tapply(bw2$perc_wt_chng,bw2$genotype,mean)
  meangeno
  
  K2       M5       M6 
  29.69014 16.76906 17.50612
  2        5        6 
  29.69014 16.76906 17.50612
  
  
  
  vartreat <- tapply(bw2$perc_wt_chng,bw2$treatment,var)
  vartreat
  1         2         3         C 
  87.80492 168.52800  28.77259 171.57621 
  sdn <- tapply(bw2$perc_wt_chng,bw2$treatment,length)
  sdn
  1  2  3  C 
  12 12 12 12 
  
  vargeno <- tapply(bw2$perc_wt_chng,bw2$genotype,var)
  vargeno
  
  K2        M5        M6 
  173.31448  43.91084  12.93573
  
  sdn <- tapply(bw2$perc_wt_chng,bw2$genotype,length)
  sdn
  
  K2 M5 M6 
  16 16 16 
  
  ##############################
  
  t.test(rikz$richness[rikz$week==2],
         rikz$richness[rikz$week==3],var.equal=T)
  
  
  t.test(bw2$perc_wt_chng[bw2$treatment==1],
         bw2$perc_wt_chng[bw2$treatment==2],var.equal=T)
  
  Two Sample t-test
  
  data:  bw2$perc_wt_chng[bw2$treatment == 1] and bw2$perc_wt_chng[bw2$treatment == 2]
  t = 0.33693, df = 22, p-value = 0.7394
  alternative hypothesis: true difference in means is not equal to 0
  95 percent confidence interval:
    -8.027831 11.142241
  sample estimates:
    mean of x mean of y 
  22.41404  20.85684 
  
  t.test(bw2$perc_wt_chng[bw2$treatment==1],
         bw2$perc_wt_chng[bw2$treatment==3],var.equal=T)
  
  Two Sample t-test
  
  data:  bw2$perc_wt_chng[bw2$treatment == 1] and bw2$perc_wt_chng[bw2$treatment == 3]
  t = 1.0884, df = 22, p-value = 0.2882
  alternative hypothesis: true difference in means is not equal to 0
  95 percent confidence interval:
    -3.071611  9.856318
  sample estimates:
    mean of x mean of y 
  22.41404  19.02169 
  
  t.test(bw2$perc_wt_chng[bw2$treatment==1],
         bw2$perc_wt_chng[bw2$treatment==4],var.equal=T)
  
  Two Sample t-test
  
  data:  bw2$perc_wt_chng[bw2$treatment == 1] and bw2$perc_wt_chng[bw2$treatment == 4]
  t = -0.12486, df = 22, p-value = 0.9018
  alternative hypothesis: true difference in means is not equal to 0
  95 percent confidence interval:
    -10.222336   9.061381
  sample estimates:
    mean of x mean of y 
  22.41404  22.99452 
  
  
  t.test(bw2$perc_wt_chng[bw2$genotype==2],
         bw2$perc_wt_chng[bw2$genotype==5],var.equal=T)
  
  Two Sample t-test
  
  data:  bw2$perc_wt_chng[bw2$genotype == 2] and bw2$perc_wt_chng[bw2$genotype == 5]
  t = 3.5067, df = 30, p-value = 0.001451
  alternative hypothesis: true difference in means is not equal to 0
  95 percent confidence interval:
    5.396038 20.446122
  sample estimates:
    mean of x mean of y 
  29.69014  16.76906 
  
  t.test(bw2$perc_wt_chng[bw2$genotype==2],
         bw2$perc_wt_chng[bw2$genotype==6],var.equal=T)
  
  Two Sample t-test
  
  data:  bw2$perc_wt_chng[bw2$genotype == 2] and bw2$perc_wt_chng[bw2$genotype == 6]
  t = 3.5711, df = 30, p-value = 0.001222
  alternative hypothesis: true difference in means is not equal to 0
  95 percent confidence interval:
    5.216117 19.151922
  sample estimates:
    mean of x mean of y 
  29.69014  17.50612 
  

  
  t.test(bw2$perc_wt_chng[bw2$genotype==5],
         bw2$perc_wt_chng[bw2$genotype==6],var.equal=T)
  
  Two Sample t-test
  
  data:  bw2$perc_wt_chng[bw2$genotype == 5] and bw2$perc_wt_chng[bw2$genotype == 6]
  t = -0.39103, df = 30, p-value = 0.6985
  alternative hypothesis: true difference in means is not equal to 0
  95 percent confidence interval:
    -4.586574  3.112453
  sample estimates:
    mean of x mean of y 
  16.76906  17.50612 
  
  ##################################
 
  
  
  wilcox.test(rikz$richness[rikz$week==2],rikz$richness[rikz$week==3])
  
  
  
  ############################
  
  
  ggplot(bw2)+geom_boxplot(aes(x=factor(treatment),y=perc_wt_chng))
  
  
  
  
  genmod2<- lm(perc_wt_chng~factor(genotype), data=BW_acropora_04252018)
  
  
  summary(genmod2)
  
  
  Call:
    lm(formula = perc_wt_chng ~ factor(genotype), data = BW_acropora_04252018)
  
  Residuals:
    Min      1Q  Median      3Q     Max 
  -15.933  -5.113  -0.546   2.166  26.210 
  
  Coefficients:
                         Estimate Std. Error t value Pr(>|t|)    
  (Intercept)          29.690      2.190  13.559  < 2e-16 ***
    factor(genotype)M5  -12.921      3.097  -4.172 0.000136 ***
    factor(genotype)M6  -12.184      3.097  -3.934 0.000286 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
  
  Residual standard error: 8.759 on 45 degrees of freedom
  Multiple R-squared:  0.328,	Adjusted R-squared:  0.2981 
  F-statistic: 10.98 on 2 and 45 DF,  p-value: 0.0001306
  
  
  treatmod2<- lm(perc_wt_chng~factor(treatment), data=BW_acropora_04252018)
  
  summary(treatmod2)
  Call:
    lm(formula = perc_wt_chng ~ factor(genotype), data = BW_acropora_04252018)
  
  Residuals:
    Min      1Q  Median      3Q     Max 
  -15.933  -5.113  -0.546   2.166  26.210 
  
  Coefficients:
    Estimate Std. Error t value Pr(>|t|)    
  (Intercept)          29.690      2.190  13.559  < 2e-16 ***
    factor(genotype)M5  -12.921      3.097  -4.172 0.000136 ***
    factor(genotype)M6  -12.184      3.097  -3.934 0.000286 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
  
  Residual standard error: 8.759 on 45 degrees of freedom
  Multiple R-squared:  0.328,	Adjusted R-squared:  0.2981 
  F-statistic: 10.98 on 2 and 45 DF,  p-value: 0.0001306
  
  > treatmod2<- lm(perc_wt_chng~factor(treatment), data=BW_acropora_04252018)
  There were 50 or more warnings (use warnings() to see the first 50)
  > summary(treatmod2)
  
  Call:
    lm(formula = perc_wt_chng ~ factor(treatment), data = BW_acropora_04252018)
  
  Residuals:
    Min      1Q  Median      3Q     Max 
  -15.258  -6.757  -3.626   2.857  35.043 
  
  Coefficients:
    Estimate Std. Error t value Pr(>|t|)    
  (Intercept)         22.4140     3.0845   7.267 4.66e-09 ***
    factor(treatment)2  -1.5572     4.3622  -0.357    0.723    
  factor(treatment)3  -3.3924     4.3622  -0.778    0.441    
  factor(treatment)C   0.5805     4.3622   0.133    0.895    
  ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
  
  Residual standard error: 10.69 on 44 degrees of freedom
  Multiple R-squared:  0.02218,	Adjusted R-squared:  -0.04448 
  F-statistic: 0.3328 on 3 and 44 DF,  p-value: 0.8017
  
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
  
  
  ggplot(bw2)+geom_boxplot(aes(x=factor(treatment),y=perc_wt_chng))
  
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
  
  
  
  
  treatmod<- lm(chng_bw~factor(treatment), data=bw2)
  
  genmod2<- lm(perc_wt_chng~factor(genotype), data=bw2)
  
  
  interactmod <- lm(richness~week.f*NAP,data=rikz) interactmod
  
  
  
  
  
  
  One Way Anova (Completely Randomized Design)
  fit <- aov(y ~ A, data=mydataframe)
  
fit1<- aov(perc_wt_chng~ treatment*genotype, data= BW_acropora_04252018)
  
fit1

summary(fit1)




#############################
boxplot(perc_wt_chng ~ treatment, data= BW_acropora_04252018,
        ylim = c(0, 100), ylab="Coral cover increase/year"))
stripchart(perc_wt_chng ~ treatment, data= BW_acropora_04252018, r5,
           vertical = TRUE, method = "jitter",
           pch = 21, col = "maroon", bg = "bisque",
           add = TRUE)
  

boxplot(perc_wt_chng ~ treatment, data= BW_acropora_04252018, outpch = NA,
        ylim = c(0, 50), ylab="Change in Coral Weight", xlab="Treatment")
stripchart(perc_wt_chng ~ treatment, data= BW_acropora_04252018,
           vertical = TRUE, method = "jitter",
           pch = 21, col = "maroon", bg = "bisque",
           add = TRUE)


boxplot(perc_wt_chng ~ genotype, data= BW_acropora_04252018, outpch = NA,
        ylim = c(0,60), ylab="Change in Coral Weight", xlab="Genotype")
stripchart(perc_wt_chng ~ genotype, data= BW_acropora_04252018,
           vertical = TRUE, method = "jitter",
           pch = 21, col = "maroon", bg = "bisque",
           add = TRUE)

install.packages("ggpubr")
library(ggpubr)

compare_means(perc_wt_chng ~ genotype, data= BW_acropora_04252018, method = "anova")

# A tibble: 1 x 6
.y.                 p    p.adj p.format p.signif method
<chr>           <dbl>    <dbl> <chr>    <chr>    <chr> 
  1 perc_wt_chng 0.000131 0.000131 0.00013  ***      Anova

compare_means(perc_wt_chng ~ treatment, data= BW_acropora_04252018, method = "anova")

# A tibble: 1 x 6
.y.              p p.adj p.format p.signif method
<chr>        <dbl> <dbl> <chr>    <chr>    <chr> 
  1 perc_wt_chng 0.802 0.802 0.8      ns       Anova


ggboxplot(BW_acropora_04252018, x = "genotype", y = "perc_wt_chng",
          palette = "jco")+
  stat_compare_means(method = "anova")

############################################



boxplot(perc_wt_chng ~ genotype, data= BW_acropora_04252018, outpch = NA,
        ylim = c(0,60), ylab="Change in Coral Weight", xlab="Genotype")
stripchart(perc_wt_chng ~ genotype, data= BW_acropora_04252018,
           vertical = TRUE, method = "jitter",
           pch = 21, col = "maroon", bg = "bisque",
           add = TRUE)+ geom_line(aes(group = interaction(index, variable)),
            alpha = 0.5, colour = "darkgrey") +
  facet_grid(sp~variable,scales="free_x") +
  scale_x_discrete(labels = "")

