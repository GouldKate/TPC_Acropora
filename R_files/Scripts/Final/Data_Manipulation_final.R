##Photosynthesis and Respiration code

rm(list=ls())


##Install packages
if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 
library(devtools)
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') install_github('colin-olito/LoLinR') 
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate') 
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron') 
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 

#Read in required libraries

##### Include Versions of libraries
#install_github('colin-olito/LoLinR')
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library('plyr')
library('dplyr')
library('tidyverse')
library('xlsx')
library('dplyr')
library('tidyr')


# get the file path

setwd("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/June 2019")
getwd()

path.p<-"C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/June 2019/Raw_R" #the location of all your respirometry files


#Change all Excel files to CSV- if you need to here is the code for it
#kate_files <- list.files("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/Presens Measurements/All", pattern="*.xlsx", full.names=TRUE)
#
#for (all in kate_files){
#  
#  n <- strsplit(all,"C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/Presens Measurements/All")[[1]][2]
#  n <- paste(strsplit(n,"//.")[[1]][1],".csv", sep = "")
#  df<-read.xlsx(all, sheetIndex = 6)
#  write.csv(df, paste("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Acropora/Presens Measurements/All/csv",n, sep = ""), 
#            row.names = FALSE, quote=FALSE)
#}

path.p<-"C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/June 2019/Raw_R" #the location of all your respirometry files

# bring in the respiration files
file.names<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders

#basename above removes the subdirectory name from the file
file.names.full<-list.files(path = path.p, pattern = "csv$", recursive = TRUE) 

#generate a 3 column dataframe with specific column names
Photo.R<- data.frame(matrix(NA, nrow=length(file.names)*2, ncol=4))
colnames(Photo.R) <- c("fragment.ID.full","Intercept", "umol.L.sec","Temp.C")
View(Photo.R)

#Load Sample Info 
library(readxl)
Sample.Info <- read_excel("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/Raw_data/Sample.Info.xlsx")#read in sample.info data
View(Sample.Info)

# load surface area data
SA <- read_excel("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/Raw_data/SA.xlsx") #read sample.info data
View(SA)

#Calculate the volume of water
as.numeric(SA$volume.mls)


# joint the sample info and surface area and volume measurements
Sample.Info<-left_join(Sample.Info, SA)
Sample.Info$date<-NULL

View(Sample.Info)

# for every file in list calculate O2 uptake or release rate and add the data to the Photo.R dataframe
for(i in 1:length(file.names.full)) { # for every file in list calculate O2 uptake or release rate and add the data to the Photo.R dataframe
  
  #find the lines in sample info that have the same file name that is being brought it
  FRow<-which(Sample.Info$fragment.ID.full==strsplit(file.names[i],'.csv'))
  # read in the O2 data one by one
  Photo.Data1 <-read.csv(file.path(path.p,file.names.full[i]), header=T) # skips the first line
  Photo.Data1 <-plyr::rename(Photo.Data1, c("Date"="Time","Oxygen"= "Value","Temperature"="Temp"))
  #rename columns 
  Photo.Data1  <- Photo.Data1[,c("Time","Value","Temp")] #subset columns of interest
  # Photo.Data1$Time <-strptime(Photo.Data1$Time, format = "%H:%M:%S") 
  Photo.Data1$Time <- strftime(Photo.Data1$Time, format = "%H:%M:%S") 
  #Photo.Data1$Time <- as.POSIXct(Photo.Data1$Time,format="%H:%M:%S", tz = "") #convert time from character to time
  Photo.Data1 <- na.omit(Photo.Data1)
  
  
  # clean up some of the data
  n<-dim(Photo.Data1)[1] # length of full data
  Photo.Data1 <-Photo.Data1[240:(n-3),] #start at data point ~4 minute in to avoid excess noise from start of run and remove last 3 lines containing text
  n<-dim(Photo.Data1)[1] #list length of trimmed data
  Photo.Data1$sec <- (1:n) #set seconds by one from start to finish of run in a new column
  
  
  #Save plot prior to and after data thinning to make sure thinning is not too extreme
  rename <- sub(".csv","", file.names[i]) # remove all the extra stuff in the file name
  
  pdf(paste0("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/Thin",rename,"thinning.pdf")) # open the graphics device
  
  par(omi=rep(0.3, 4)) #set size of the outer margins in inches
  par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
  plot(Value ~ sec, data=Photo.Data1 , xlab='Time (seconds)', ylab=expression(paste(' O'[2],' (',mu,'mol/L)')), type='n', axes=FALSE) #plot (empty plot to fill) data as a function of time
  usr  <-  par('usr') # extract the size of the figure margins
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA) # put a grey background on the plot
  whiteGrid() # make a grid
  box() # add a box around the plot
  points(Photo.Data1 $Value ~ Photo.Data1 $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1) # add the x axis
  axis(2, las=1) # add the y-axis
  
  # Thin the data to make the code run faster
  Photo.Data.orig<-Photo.Data1#save original unthinned data
  Photo.Data1 <-  thinData(Photo.Data1 ,by=20)$newData1 #thin data by every 20 points for all the O2 values
  Photo.Data1$sec <- as.numeric(rownames(Photo.Data1 )) #maintain numeric values for time
  Photo.Data1$Temp<-NA # add a new column to fill with the thinned data
  Photo.Data1$Temp <-  thinData(Photo.Data.orig,xy = c(1,3),by=20)$newData1[,2] #thin data by every 20 points for the temp values
  
  # plot the thinned data
  plot(Value ~ sec, data=Photo.Data1 , xlab='Time (seconds)', ylab=expression(paste(' O'[2],' (',mu,'mol/L)')), type='n', axes=FALSE) #plot thinned data
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  box()
  points(Photo.Data1 $Value ~ Photo.Data1 $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1)
  axis(2, las=1)
  ##Olito et al. 2017: It is running a bootstrapping technique and calculating the rate based on density
  #option to add multiple outputs method= c("z", "eq", "pc")
  Regs  <-  rankLocReg(xall=Photo.Data1$sec, yall=Photo.Data1$Value, alpha=0.5, method="pc", verbose=TRUE)  
  
  # add the regression data
  plot(Regs)
  dev.off()
  
  # fill in all the O2 consumption and rate data
  Photo.R[i,2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
  Photo.R[i,1] <- rename #stores the file name in the Date column
  Photo.R[i,4] <- mean(Photo.Data1$Temp, na.rm=T)  #stores the Temperature in the Temp.C column
  #Photo.R[i,5] <- PR[j] #stores stores whether it is photosynthesis or respiration
  
  
  # rewrite the file everytime... I know this is slow, but it will save the data that is already run

  
View(Photo.R)
  
Photo.R <- na.omit(Photo.R)
  
write.csv(Photo.R,"C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/Photo.R.csv") 
  

Photo.R<-read.csv("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/Photo.R.csv")
Photo.R$X<-NULL 


#join dataframes
Photo.R<-left_join(Photo.R, Sample.Info)
View(Photo.R)

# Calculate P and R rate
#Convert sample volume to mL
Photo.R$volume <- Photo.R$volume.mls/1000 #calculate volume

#Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
Photo.R$umol.sec <- Photo.R$umol.L.sec*Photo.R$volume

#Account for blank rate by temperature
#convert character columns to factors
Photo.R <- Photo.R %>%
  mutate_if(sapply(., is.character), as.factor)
View(Photo.R)


#make the blank column a factor
Photo.R$BLANK<-ifelse(Photo.R$treatment=='blank', 1,0)
Photo.R$BLANK<-as.factor(Photo.R$BLANK)
View(Photo.R)



photo.blnk <- aggregate(umol.sec ~ species*temp.Cat*light_dark*run*BLANK, data=Photo.R, mean)
View(photo.blnk)
# pull out only the blanks
#photo.blnk<-photo.blnk[photo.blnk$Species=='BK',]
photo.blnk<-photo.blnk[photo.blnk$BLANK==1,]
# remove the species column and join with the full data set
photo.blnk$species<-NULL
# remove the blank column
photo.blnk$BLANK<-NULL


colnames(photo.blnk)[4]<-'blank.rate' # rename the blank rate 

write.csv(photo.blnk, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/photo.blk.csv") 

# join the blank data with the rest of the data
Photo.R<-left_join(Photo.R, photo.blnk)
View(Photo.R)




# subtract the blanks######################
Photo.R$umol.sec.corr<-Photo.R$umol.sec-Photo.R$blank.rate

View(Photo.R)

#### Normalize to SA (surface area)#####

Photo.R$umol.cm2.hr <- (Photo.R$umol.sec.corr*3600)/Photo.R$surf.area.cm2 #mmol cm-2 hr-1

#Photo.R<-Photo.R[complete.cases(Photo.R),] # remove NAs and blanks
Photo.R<-Photo.R[Photo.R$BLANK==0,]

Photo.R <- na.omit(Photo.R)
#Photo.Rtest<-Photo.Rtest[complete.cases(Photo.Rtest),] # # remove NAs and blanks

#Calculate net P and R
Photo.Rtest<-Photo.Rtest[Photo.Rtest$BLANK==0,]
Photo.Rfinal<-Photo.R
write.csv(Photo.Rfinal,"C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/Final_Photo.R.csv") 

### Raw data visualization- Extensively NP and Respiration
photo.Rrrr<-Photo.Rfinal

photo.Rrrr$umol.cm2.hr[photo.Rrrr$light_dark=="dark"]<-photo.Rrrr$umol.cm2.hr[photo.Rrrr$light_dark=="dark"]*-1 
lessthan <- which(photo.Rrrr$light_dark=="dark" & photo.Rrrr$umol.cm2.hr < 0)
photo.Rrrr$umol.cm2.hr[lessthan] <- 0
#####
raw1<-ggplot(photo.Rrrr, aes(x=Temp.C, y=umol.cm2.hr,group = light_dark, col = genotype))+
  geom_point()+  
  facet_wrap(~ treatment*light_dark, labeller = labeller(.multi_line = FALSE))

esquisse:::ggplot_to_ppt("raw1")

raw2<-ggplot(photo.Rrrr, aes(x=Temp.C, y=umol.cm2.hr,group = light_dark, col = genotype))+
  geom_point()+  
  facet_wrap(~ treatment*light_dark*genotype, labeller = labeller(.multi_line = FALSE))

esquisse:::ggplot_to_ppt("raw2")

raw3<-ggplot(photo.Rrrr, aes(x=Temp.C, y=umol.cm2.hr, group = light_dark, col = genotype))+
  geom_point()+  
  facet_wrap(~ light_dark*genotype, labeller = labeller(.multi_line = FALSE))

esquisse:::ggplot_to_ppt("raw3")


######## SPECIES
raw1<-ggplot(photo.Rrrr, aes(x=Temp.C, y=umol.cm2.hr,group = light_dark, col = species))+
geom_point()+  
  facet_wrap(~ treatment*light_dark, labeller = labeller(.multi_line = FALSE))

esquisse:::ggplot_to_ppt("raw1")

raw2<-ggplot(photo.Rrrr, aes(x=Temp.C, y=umol.cm2.hr,group = light_dark, col = species))+
  geom_point()+  
  facet_wrap(~ treatment*light_dark*species, labeller = labeller(.multi_line = FALSE))

esquisse:::ggplot_to_ppt("raw2")


raw3<-ggplot(photo.Rrrr, aes(x=Temp.C, y=umol.cm2.hr, group = light_dark, col = species))+
  geom_point()+  
  facet_wrap(~ light_dark*species, labeller = labeller(.multi_line = FALSE))

esquisse:::ggplot_to_ppt("raw3")
##########################


write.csv(photo.Rrrr, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/PhotoPrrr.csv") 


#Restart and look at FINAL data frame and calculate GROSS PHOTOSYNTHESIS
#
Photo.R2<-Photo.Rfinal

###calculating gross photosynthesis from data frame###

#make ifelse statements to assign light treatments as NP and dark treatments as resp
#light will be assigned NP for net photosynthesis 

Photo.R2$rate.type <-ifelse(Photo.R2$light_dark=='light', "NP", "R")
Photo.R2$rate.type<-as.factor(Photo.R2$rate.type)
View(Photo.R2)


#rename fragment ID 
Photo.R2$individual.ID <- str_split(Photo.R2$fragment.ID, "D", n = Inf, simplify = TRUE)[,1]
Photo.R2$individual.ID <- as.factor(Photo.R2$individual.ID)

library(stringr)
Photo.R2$individual.ID<-str_sub(Photo.R2$individual.ID , 1, str_length(Photo.R2$individual.ID)-1)

write.csv(Photo.R2, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/Photo.R2NP_R.csv") 

Pnet2<- subset(Photo.R2, rate.type=="NP")
Rdark2<- subset(Photo.R2, rate.type == "R")



### to get umol.cm2.hr in ONE row (subset rows then re-combine)
pgoss<-rbind(Pnet2,Rdark2)

#Data was not measured correctly for the 2602 dark Trial in the BLANK CHamber- 
#remove trial 2602D 
pgoss%<>%
  dplyr::filter(run!="2602 dark")
View(pgoss)  

Photo.Res$filter.column<-NULL

### pgoss is Pnet- Rdark when Rdark is negative

write.csv(pgoss, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/np_R_GP.csv")
##### made Grossphoto, by Pnet-Rdark in pgoss in excel and saving it as Pgross_Pnet-Rdark.csv
##### ### (respiration is negative so Pnet-Res=GP)- if positive Np+Res= GP

# re-read file into make GP graphs and compare NP and R rates

resp.datafin<-read.csv("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/Pgross_Pnet-Rdark.csv")
#### This data frame has NPO, R, and GP ALL in one row for each Fragment.ID






GP<-ggplot(resp.datafin, aes(x=Temp.C, y=GP_umol.cm2.hr, group =fragment.ID, col = genotype))+
  geom_point()+   
  theme_bw () +
  ylim(-2,4)+  
  facet_wrap(~ genotype*treatment, labeller = labeller(.multi_line = FALSE))


esquisse:::ggplot_to_ppt("GP")

#### combined Pgross_Pnet-Rdark.csv, ad np_R_GP.csv to make one column on umol.cm.hr and distinguish via rate.type
resp.datafin2<-read.csv("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/final_NP_R_GP_edited.csv")


GP<-ggplot(resp.datafin2, aes(x=Temp.C, y=umol.cm2.hr, group = fragment.ID, col = genotype))+
  geom_point()+   
  theme_bw () +
  #ylim(0,1.5)+  
  facet_wrap(~ rate.type, labeller = labeller(.multi_line = FALSE))

esquisse:::ggplot_to_ppt("GP")


GP<-ggplot(resp.datafin2, aes(x=Temp.C, y=umol.cm2.hr, group = rate.type, col = genotype))+
  geom_point()+   
  theme_bw () +
  #ylim(0,1.5)+  
  facet_wrap(~ treatment*rate.type, labeller = labeller(.multi_line = FALSE))


esquisse:::ggplot_to_ppt("GP")

GP<-ggplot(resp.datafin2, aes(x=Temp.C, y=umol.cm2.hr, group = rate.type, col = genotype))+
  geom_point()+   
  theme_bw () +
  #ylim(0,1.5)+  
  facet_wrap(~ treatment*species, labeller = labeller(.multi_line = FALSE))


esquisse:::ggplot_to_ppt("GP")

esquisse:::ggplot_to_ppt("GP")

GP<-ggplot(resp.datafin2, aes(x=Temp.C, y=umol.cm2.hr, group = rate.type, col = genotype))+
  geom_point()+   
  theme_bw () +
  #ylim(0,1.5)+  
  facet_wrap(~ treatment*genotype, labeller = labeller(.multi_line = FALSE))


esquisse:::ggplot_to_ppt("GP")

GP<-ggplot(resp.datafin2, aes(x=Temp.C, y=umol.cm2.hr, group = rate.type, col = genotype))+
  geom_point()+   
  theme_bw () +
  #ylim(0,1.5)+  
  facet_wrap(~ treatment*species*rate.type, labeller = labeller(.multi_line = FALSE))


esquisse:::ggplot_to_ppt("GP")



#Calculate means from resp.datafin with ALL NP, R, and GP umo.cm2.hr on it
AllMeansG <- ddply(resp.datafin, c('genotype'), summarise,
                   #pnet
                   Pnet.mean= mean(Pnet_umol.cm2.hr, na.rm=T), #mean pnet
                   N = sum(!is.na(Pnet_umol.cm2.hr)), # sample size
                   Pnet.se = sd(Pnet_umol.cm2.hr, na.rm=T)/sqrt(N), #SE
                   #Rdark
                   Rdark.mean= mean(Rdark_umol.cm2.hr, na.rm=T), #mean rdark
                   Rdark.se = sd(Rdark_umol.cm2.hr, na.rm=T)/sqrt(N), #SE
                   #Pgross
                   Pgross.mean  = mean(GP_umol.cm2.hr, na.rm=TRUE),
                   Pgross.se = sd(GP_umol.cm2.hr, na.rm=TRUE)/sqrt(N),
                   Temp.mean = mean(Temp.C, na.rm=TRUE),
                   Temp.se = sd(Temp.C, na.rm=TRUE)/sqrt(N)
)

View(AllMeansG)

write.csv(AllMeansG, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/AllMeansGfin.csv") 


AllMeansS <- ddply(resp.datafin, c('species','temp.Cat', 'treatment'), summarise,
                   #pnet
                   Pnet.mean= mean(Pnet_umol.cm2.hr, na.rm=T), #mean pnet
                   N = sum(!is.na(Pnet_umol.cm2.hr)), # sample size
                   Pnet.se = sd(Pnet_umol.cm2.hr, na.rm=T)/sqrt(N), #SE
                   #Rdark
                   Rdark.mean= mean(Rdark_umol.cm2.hr, na.rm=T), #mean rdark
                   Rdark.se = sd(Rdark_umol.cm2.hr, na.rm=T)/sqrt(N), #SE
                   #Pgross
                   Pgross.mean  = mean(GP_umol.cm2.hr, na.rm=TRUE),
                   Pgross.se = sd(GP_umol.cm2.hr, na.rm=TRUE)/sqrt(N),
                   Temp.mean = mean(Temp.C, na.rm=TRUE),
                   Temp.se = sd(Temp.C, na.rm=TRUE)/sqrt(N)
)

View(AllMeansS)

write.csv(AllMeansS, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/AllMeansSfin.csv") 


AllMeansdepth <- ddply(resp.datafin, c('treatment','temp.Cat'), summarise,
                       #pnet
                       Pnet.mean= mean(Pnet_umol.cm2.hr, na.rm=T), #mean pnet
                       N = sum(!is.na(Pnet_umol.cm2.hr)), # sample size
                       Pnet.se = sd(Pnet_umol.cm2.hr, na.rm=T)/sqrt(N), #SE
                       #Rdark
                       Rdark.mean= mean(Rdark_umol.cm2.hr, na.rm=T), #mean rdark
                       Rdark.se = sd(Rdark_umol.cm2.hr, na.rm=T)/sqrt(N), #SE
                       #Pgross
                       Pgross.mean  = mean(GP_umol.cm2.hr, na.rm=TRUE),
                       Pgross.se = sd(GP_umol.cm2.hr, na.rm=TRUE)/sqrt(N),
                       Temp.mean = mean(Temp.C, na.rm=TRUE),
                       Temp.se = sd(Temp.C, na.rm=TRUE)/sqrt(N)
)



View(AllMeansdepth)

write.csv(AllMeansdepth, "C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/AllMeansdepthfin.csv") 




### GRaph Pgross means by species 
GPplot<-ggplot(AllMeansS, aes(x=temp.Cat, y=Pgross.mean, color=species))+
  geom_point()+
  facet_wrap(~species)

esquisse:::ggplot_to_ppt("GPplot")

Rrplot<-ggplot(AllMeansS, aes(x=temp.Cat, y=Rdark.mean, color=species))+
  geom_point()+
  facet_wrap(~species)

esquisse:::ggplot_to_ppt("Rrplot")

Nplot<-ggplot(AllMeansS, aes(x=temp.Cat, y=Pnet.mean, color=species))+
  geom_point()+
  facet_wrap(~species)

esquisse:::ggplot_to_ppt("Nplot")

#### By treatment

Rplot<-ggplot(AllMeansdepth, aes(x=temp.Cat, y=Rdark.mean, color=treatment))+
  geom_point()+
  facet_wrap(~treatment)

esquisse:::ggplot_to_ppt("Rplot")

Gplot<-ggplot(AllMeansdepth, aes(x=temp.Cat, y=Pgross.mean, color=temp.Cat))+
  geom_point()+
  facet_wrap(~treatment)

esquisse:::ggplot_to_ppt("Gplot")

NPplot<-ggplot(AllMeansdepth, aes(x=temp.Cat, y=Pnet.mean, color=temp.Cat))+
  geom_point()+
  facet_wrap(~treatment)

esquisse:::ggplot_to_ppt("NPplot")

####
# by genotype
# 
#AllMeansG$treatment<-str_sub(AllMeansG$genotype, 1, str_length(AllMeansG$genotype)-2)
#AllMeansG$species<-str_sub(AllMeansG$genotype, 1, str_length(AllMeansG$genotype)-1)
  
AllMeansG<-read.csv("C:/Users/Kate/OneDrive - University of North Carolina at Chapel Hill/2019laptop/Documents/TPC_Mesop_BDA/R_output/Final/AllMeansGfin.csv")

Rplot<-ggplot(AllMeansG, aes(x=X.1, y=Rdark.mean, color=genotype, shape=X.1))+
  geom_point()+
  facet_wrap(~treatment)

esquisse:::ggplot_to_ppt("Rplot")

Gplot<-ggplot(AllMeansG, aes(x=X.1, y=Pgross.mean, color=genotype, shape=X.1))+
  geom_point()+
  facet_wrap(~treatment)

esquisse:::ggplot_to_ppt("Gplot")

NPplot<-ggplot(AllMeansG, aes(x=X.1, y=Pnet.mean, color=genotype, shape=X.1))+
  geom_point()+
  facet_wrap(~treatment)

esquisse:::ggplot_to_ppt("NPplot")



#### 
#### 

GPtbT<-ggplot(subset(resp.datafin2,rate.type=="GP"), aes(x=treatment, umol.cm2.hr, colour=genotype, shape=treatment))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=2)+
  geom_point(position="dodge", size=2)+
  ylab("GP(umol.cm2.hr)")+
  labs(title="Mean GP across Treatment")+
  coord_flip()
GPtbT
esquisse:::ggplot_to_ppt("GPtbT")

GPtb<-ggplot(subset(resp.datafin2,rate.type=="GP"), aes(x=genotype, umol.cm2.hr, colour=treatment, shape=genotype))+
  stat_summary(fun.data = "mean_cl_boot", colour = "black", size = 2)+
  geom_point(size=2)+
  #geom_errorbar(AllMeansS,aes(ymax=mean+se, ymin=mean-se), size=1.2, width=0.4) +
  ylab("GP(umol.cm2.hr)")+
  labs(title="Mean GP across genotype")+
  coord_flip()
GPtb
esquisse:::ggplot_to_ppt("GPtb")


