
###############################
install.packages("ggplot2")
library(ggplot2)
##################################
BW_acropora_04252018
head(BW_acropora_04252018)

bw<-(BW_acropora_04252018)

head(bw)

ggplot(bw)+geom_boxplot(aes(x=factor(treatment),y=perc_wt_chng))

ggplot(bw)+geom_boxplot(aes(x=treatment,y=perc_wt_chng))+ 
  labs(x="Treatment",y="Percent Weight Change")

ggplot(bw)+geom_boxplot(aes(x=factor(genotype),y=perc_wt_chng))

ggplot(bw)+geom_boxplot(mapping=aes(x=genotype,y=perc_wt_chng))

ggplot(bw)+geom_point(aes(x=factor(coralid),y=perc_wt_chng))

ggplot(bw)+geom_point(aes(x=genotype,y=perc_wt_chng))

ggplot(bw)+geom_boxplot(aes(x=genotype,y=perc_wt_chng))+ 
  labs(x="Genotype",y="Percent Weight Change")
##########################################
meanperchan <- mean(bw$perc_wt_chng)
meanperchan
sd(bw$perc_wt_chng)
seperchan <- sd(bw$perc_wt_chng)/length(bw$perc_wt_chng)
seperchan

dnorm(0)
pnorm(0)


qnorm(.975)
pnorm(qnorm(.975))

serich*qnorm(.975)
qnorm(.975, sd=serich)


c(meanperchan+seperchan*qnorm(.975),meanperchan-seperchan*qnorm(.975))

c(meanperchan+seperchan*qt(.975,length(bw$perc_wt_chng)-1),
  meanperchan-seperchan*qt(.975,length(bw$perc_wt_chng-1))

## didn't work


meantreatment <- tapply(bw$perc_wt_chng,bw$treatment,mean)
meantreatment

1        2        3        C 
22.41404 20.85684 19.02169 22.99452 

meangeno <- tapply(bw$perc_wt_chng,bw$genotype,mean)
meangeno

K2       M5       M6 
29.69014 16.76906 17.50612

  

vartreat <- tapply(bw$perc_wt_chng,bw$treatment,var)
vartreat
1         2         3         C 
87.80492 168.52800  28.77259 171.57621 
sdn <- tapply(bw$perc_wt_chng,bw$treatment,length)
sdn
1  2  3  C 
12 12 12 12 

vargeno <- tapply(bw$perc_wt_chng,bw$genotype,var)
vargeno

K2        M5        M6 
173.31448  43.91084  12.93573

sdn <- tapply(bw$perc_wt_chng,bw$genotype,length)
sdn

K2 M5 M6 
16 16 16 

##############################

t.test(rikz$richness[rikz$week==2],
       rikz$richness[rikz$week==3],var.equal=T)


t.test(bw$perc_wt_chng[bw$treatment==1],
       bw$perc_wt_chng[bw$treatment==2],var.equal=T)

Two Sample t-test

data:  bw$perc_wt_chng[bw$treatment == 1] and bw$perc_wt_chng[bw$treatment == 2]
t = 0.33693, df = 22, p-value = 0.7394
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -8.027831 11.142241
sample estimates:
  mean of x mean of y 
22.41404  20.85684 

t.test(bw$perc_wt_chng[bw$treatment==1],
       bw$perc_wt_chng[bw$treatment==3],var.equal=T)

Two Sample t-test

data:  bw$perc_wt_chng[bw$treatment == 1] and bw$perc_wt_chng[bw$treatment == 3]
t = 1.0884, df = 22, p-value = 0.2882
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -3.071611  9.856318
sample estimates:
  mean of x mean of y 
22.41404  19.02169 

t.test(bw$perc_wt_chng[bw$treatment==1],
       bw$perc_wt_chng[bw$treatment==C],var.equal=T)

t.test(bw$perc_wt_chng[bw$treatment==2],
       bw$perc_wt_chng[bw$treatment==3],var.equal=T)

t.test(bw$perc_wt_chng[bw$treatment==2],
       bw$perc_wt_chng[bw$treatment==C],var.equal=T)

t.test(bw$perc_wt_chng[bw$treatment==3],
       bw$perc_wt_chng[bw$treatment==C],var.equal=T)

##################################
t.test(bw$perc_wt_chng[bw$genotype==K2],
       bw$perc_wt_chng[bw$genotype==M5],var.equal=T)

t.test(bw$perc_wt_chng[bw$genotype==K2],
       bw$perc_wt_chng[bw$genotype==M6],var.equal=T)

t.test(bw$perc_wt_chng[bw$genotype==M6],
       bw$perc_wt_chng[bw$genotype==M5],var.equal=T)


    ################
boxplot(Recovery ~ Recovery_group, data=recovery_data, outpch = NA,
        ylim = c(0, 100), ylab="Coral cover increase/year")
stripchart(Recovery ~ Recovery_group, data=recovery_data,
           vertical = TRUE, method = "jitter",
           pch = 21, col = "maroon", bg = "bisque",
           add = TRUE)



