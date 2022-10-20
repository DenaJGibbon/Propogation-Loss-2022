# Notes to Dena: Need to add time categories

# Noise plots
library(ggpubr)
library(stringr)
library(plyr)
library(tidyverse)
library(seewave)
library(data.table)
library(seewave)
library(nlme)
library(ggplot2)
library(MASS)
library(multcomp)
library(grid)
library(reshape2)
library(plyr)
library(lme4)

#standard error function
se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))


RunganNoise <-
  read.csv("/Users/denaclink/Desktop/RStudio Projects/Propogation-Loss-2022/ThirdOctaveBandDFRunganAugust1.csv")

# Match .wav file and habitat
RunganMetaData <- read.csv('PropLoss_test_7Jul22.csv')

RunganNoiseAddHabitat <- data.frame()
unique.wav <- unique(RunganNoise$wav.file)
for(a in 1:length(unique.wav)){
  unique.wav.single <-  unique.wav[a]
  unique.wav.single <- str_split_fixed(unique.wav.single,pattern = '.wav',n=2)[,1]
  wav.index <- which(str_detect(RunganMetaData$Sel_Table,unique.wav.single))
  RunganMatchRow <- RunganMetaData[wav.index[1],]
  RunganNoiseSub <- subset(RunganNoise, wav.file==unique.wav[a])
  RunganNoiseSub$habitat <- rep(RunganMatchRow$Habitat,nrow(RunganNoiseSub))
  RunganNoiseAddHabitat <- rbind.data.frame(RunganNoiseAddHabitat,RunganNoiseSub)
  }

RunganNoiseAddHabitat$time <- str_split_fixed(RunganNoiseAddHabitat$wav.file,pattern='_',n=3)[,3]
RunganNoiseAddHabitat$time <- as.numeric(substr(RunganNoiseAddHabitat$time,1,2))

RunganNoiseAddHabitat <- subset(RunganNoiseAddHabitat, time > 5 & time < 17)

RunganNoiseMedian <- RunganNoiseAddHabitat %>%
  group_by(habitat,time,center.freq) %>%
  dplyr::summarise(ambient.noise.high=meandB(noise.valuedb)+se(noise.valuedb),
            ambient.noise.low=meandB(noise.valuedb)-se(noise.valuedb),
            ambient.noise=meandB(noise.valuedb))

# Drop missing habitat
RunganNoiseMedian <- subset(RunganNoiseMedian,habitat=='K'| habitat=='LP' | habitat=='MS')

RunganNoiseMedian$time <- as.factor(RunganNoiseMedian$time)
RunganNoiseMedian$time <- as.factor(RunganNoiseMedian$time)

RunganNoisePlot <- ggplot(RunganNoiseMedian,aes(center.freq,ambient.noise,
                       group=time,colour=time,linetype=time))+
  stat_summary(data=RunganNoiseMedian,fun.y=meandB,geom="line",alpha=0.2,aes(group=time))+
  stat_summary(data=RunganNoiseMedian,fun.y=meandB,geom="line",aes(group=time))+
  geom_ribbon(data=RunganNoiseMedian,aes(ymin=ambient.noise.low,ymax=ambient.noise.high
                                         ,fill=time,group=time,color=NULL),alpha=0.25)+ 
  scale_fill_manual(values= matlab::jet.colors(length(unique(RunganNoiseMedian$time))) )+ theme_bw()+
  scale_color_manual(values= matlab::jet.colors(length(unique(RunganNoiseMedian$time))) )+ theme_bw()#+ylim(25,70)
#+ylim(25,70)

table(RunganNoiseMedian$time)

MaliauNoise <-
  read.csv("/Users/denaclink/Desktop/RStudio Projects/Propogation-Loss-2022/ThirdOctaveBandDFMaliauAddtimeswith20.csv")

MaliauNoise$time <- str_split_fixed(MaliauNoise$wav.file,pattern='_',n=3)[,3]
MaliauNoise$time <- as.numeric(substr(MaliauNoise$time,1,2))

recorders <- str_split_fixed(MaliauNoise$wav.file,pattern = '_',n=3)[,1]

# Remove close weird ones
MaliauNoise <-
  MaliauNoise[which(recorders=='M6' |recorders=='M7' | recorders=='M8' ),]

head(MaliauNoise)

MaliauNoiseMedian <- MaliauNoise %>%
  group_by(time,center.freq) %>%
  dplyr::summarise(ambient.noise.high=meandB(noise.valuedb)+se(noise.valuedb),
            ambient.noise.low=meandB(noise.valuedb)-se(noise.valuedb),
            ambient.noise=meandB(noise.valuedb))

MaliauNoise$time <- as.factor(MaliauNoise$time)
MaliauNoiseMedian$time <- as.factor(MaliauNoiseMedian$time)

MaliauNoisePlot <- ggplot(MaliauNoiseMedian,aes(center.freq,ambient.noise,
                          group=time,colour=time,linetype=time))+
  stat_summary(data=MaliauNoiseMedian,fun.y=meandB,geom="line",alpha=0.2,aes(group=time))+
  stat_summary(data=MaliauNoiseMedian,fun.y=meandB,geom="line",aes(group=time))+
  geom_ribbon(data=MaliauNoiseMedian,aes(ymin=ambient.noise.low,ymax=ambient.noise.high
                                  ,fill=time,group=time,color=NULL),alpha=0.25)+ 
  scale_fill_manual(values= matlab::jet.colors(length(unique(MaliauNoise$time))) )+ 
  scale_color_manual(values= matlab::jet.colors(length(unique(MaliauNoise$time))) )+ theme_bw()#+ylim(25,70)

cowplot::plot_grid(RunganNoisePlot,MaliauNoisePlot)

## Just do by habitat type

MaliauNoise$habitat <- rep('D',nrow(MaliauNoise))

CombinedNoise <- 
  rbind.data.frame(RunganNoiseAddHabitat,MaliauNoise)

CombinedNoiseMedian <- CombinedNoise %>%
  group_by(habitat,center.freq) %>%
  dplyr::summarise(ambient.noise.high=meandB(noise.valuedb)+se(noise.valuedb),
            ambient.noise.low=meandB(noise.valuedb)-se(noise.valuedb),
            ambient.noise=meandB(noise.valuedb))



# Drop missing habitat
CombinedNoiseMedian <- subset(CombinedNoiseMedian,habitat=='D'|habitat=='K'| habitat=='LP' | habitat=='MS')


ggplot(CombinedNoiseMedian,aes(center.freq,ambient.noise,
                                                group=habitat,colour=habitat,linetype=habitat))+
  stat_summary(data=CombinedNoiseMedian,fun.y=meandB,geom="line",aes(group=habitat),lwd=1)+
   stat_summary(data=CombinedNoiseMedian,fun.y=meandB,geom="line",aes(group=habitat))+
  # geom_ribbon(data=CombinedNoiseMedian,aes(ymin=ambient.noise.low,ymax=ambient.noise.high
  #                                        ,fill=habitat,group=habitat,color=NULL),alpha=0.5)+ 
  # scale_fill_manual(values= matlab::jet.colors(length(unique(CombinedNoise$habitat))) )+ 
  scale_color_manual(values= matlab::jet.colors(length(unique(CombinedNoise$habitat))) )+ 
  theme_bw()+ylab(expression(paste('Ambient sound level (dB re 20', mu,'Pa)',sep=' ')))+xlab('Center Frequency (Hz)') #+ylim(25,70)

ggplot(CombinedNoiseMedian,aes(center.freq,ambient.noise,
                               group=habitat,colour=habitat,linetype=habitat))+
  stat_summary(data=CombinedNoiseMedian,fun.y=meandB,geom="line",aes(group=habitat),lwd=1)+
  stat_summary(data=CombinedNoiseMedian,fun.y=meandB,geom="line",aes(group=habitat))+
   geom_ribbon(data=CombinedNoiseMedian,aes(ymin=ambient.noise.low,ymax=ambient.noise.high
                                           ,fill=habitat,group=habitat,color=NULL),alpha=0.5)+ 
  scale_fill_manual(values= matlab::jet.colors(length(unique(CombinedNoise$habitat))) )+ 
  scale_color_manual(values= matlab::jet.colors(length(unique(CombinedNoise$habitat))) )+ 
  theme_bw()+ylab(expression(paste('Ambient sound level (dB re 20', mu,'Pa)',sep=' ')))+
  xlab('Center Frequency (Hz)')+xlim(150,2000)+ylim(25,45)

      
                  

#CombinedNoiseMedian$time <- as.numeric(CombinedNoiseMedian$time)

ggplot(CombinedNoiseMedian,aes(center.freq,ambient.noise,
                               group=habitat,colour=habitat,linetype=habitat))+
  stat_summary(data=CombinedNoiseMedian,fun.y=meandB,geom="line",alpha=0.2,aes(group=habitat))+
  stat_summary(data=CombinedNoiseMedian,fun.y=meandB,geom="line",aes(group=habitat))+
  geom_ribbon(data=CombinedNoiseMedian,aes(ymin=ambient.noise.low,ymax=ambient.noise.high
                                           ,fill=habitat,group=habitat,color=NULL),alpha=0.25)+ 
  scale_fill_manual(values= matlab::jet.colors(length(unique(CombinedNoiseMedian$habitat))) )+ 
  scale_color_manual(values= matlab::jet.colors(length(unique(CombinedNoiseMedian$habitat))) )+ 
  theme_bw()+ylab(expression(paste('Ambient sound level (dB re 20', mu,'Pa)',sep=' ')))+
  xlab('Center Frequency (Hz)')

ggplot(CombinedNoiseMedian,aes(center.freq,ambient.noise,
                               group=habitat,colour=habitat,linetype=habitat))+
  stat_summary(data=CombinedNoiseMedian,fun.y=meandB,geom="line",alpha=0.2,aes(group=habitat))+
  stat_summary(data=CombinedNoiseMedian,fun.y=meandB,geom="line",aes(group=habitat))+
  geom_ribbon(data=CombinedNoiseMedian,aes(ymin=ambient.noise.low,ymax=ambient.noise.high
                                           ,fill=habitat,group=habitat,color=NULL),alpha=0.25)+ 
  scale_fill_manual(values= matlab::jet.colors(length(unique(CombinedNoiseMedian$habitat))) )+ 
  scale_color_manual(values= matlab::jet.colors(length(unique(CombinedNoiseMedian$habitat))) )+ 
  theme_bw()+ylab(expression(paste('Ambient sound level (dB re 20', mu,'Pa)',sep=' ')))+
  xlab('Center Frequency (Hz)')+xlim(0,2000)+ylim(25,45)



# Model selection ---------------------------------------------------------
CombinedNoiseSub <- subset(CombinedNoise,center.freq > 50 & center.freq < 2000)
CombinedNoiseSub <- subset(CombinedNoiseSub,habitat=='D'|habitat=='K'| habitat=='LP' | habitat=='MS')

CombinedNoiseMedian <- CombinedNoiseSub %>%
  group_by(habitat,center.freq,time) %>%
  summarise(ambient.noise.high=meandB(noise.valuedb)+se(noise.valuedb),
            ambient.noise.low=meandB(noise.valuedb)-se(noise.valuedb),
            ambient.noise=meandB(noise.valuedb))


CombinedNoiseMedian$TimeCat <- recode(CombinedNoiseMedian$time, '16' = "Afternoon",
                                      '10'="Morning",'11'='Morning','13'="Afternoon",
                                      '14'="Afternoon",'15'="Afternoon",'6'='Morning',
                                      '7'='Morning','8'='Morning','9'='Morning')

CombinedNoiseMedian$Site <- recode(CombinedNoiseMedian$habitat, 'K' = "Rungan",
                                      'LP'="Rungan",'MS'='Rungan','D'="Maliau")

CombinedNoiseMedian$center.freq <- as.factor(CombinedNoiseMedian$center.freq)

Combined.lmerm.prop.loss.null <- lmer(ambient.noise ~  (1|Site), data=CombinedNoiseMedian) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmerm.prop.loss.center.freq <- lmer(ambient.noise ~  center.freq + (1|Site), data=CombinedNoiseMedian) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmerm.prop.habitat <- lmer(ambient.noise ~  habitat+ (1|Site), data=CombinedNoiseMedian) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmerm.prop.time <- lmer(ambient.noise ~  TimeCat+ (1|Site), data=CombinedNoiseMedian) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmerm.prop.loss.full <- lmer(ambient.noise ~  center.freq + habitat+ TimeCat+ (1|Site), data=CombinedNoiseMedian) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmerm.prop.loss.nocenterfreq <- lmer(ambient.noise ~  habitat+ TimeCat+ (1|Site), data=CombinedNoiseMedian) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmerm.prop.loss.interaction <- lmer(ambient.noise ~  habitat*TimeCat+ (1|Site), data=CombinedNoiseMedian) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmerm.prop.loss.interaction.cf <- lmer(ambient.noise ~ center.freq + habitat*TimeCat+ (1|Site), data=CombinedNoiseMedian) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmerm.prop.loss.interaction.habitat <- lmer(ambient.noise ~ center.freq *habitat + TimeCat+ (1|Site), data=CombinedNoiseMedian) # + (Call.Type|recorder.ID + Call.Type|recorder.location)

sjPlot::plot_model(Combined.lmerm.prop.loss.interaction.habitat,intercept=F)+ggtitle('Ambient noise model coefficents')+theme_bw()+
  geom_hline(yintercept = 0,linetype='dashed')

bbmle::AICctab(Combined.lmerm.prop.loss.null,Combined.lmerm.prop.loss.center.freq,
               Combined.lmerm.prop.habitat,Combined.lmerm.prop.time,Combined.lmerm.prop.loss.full,Combined.lmerm.prop.loss.nocenterfreq,
               Combined.lmerm.prop.loss.interaction,Combined.lmerm.prop.loss.interaction.cf,Combined.lmerm.prop.loss.interaction.habitat,weights=T)

hist(resid(Combined.lmerm.prop.loss.interaction))

CombinedNoiseMedian500 <- subset(CombinedNoiseMedian,center.freq==250 |center.freq==500 |center.freq==500 | center.freq==1250 |
                                   center.freq==1600)

ggboxplot(data=CombinedNoiseMedian,x='habitat',y='ambient.noise',
          fill  = 'center.freq',outlier.shape = NA)+xlab('Habitat')+ylab(expression(paste('Ambient sound level (dB re 20', mu,'Pa)',sep=' ')))

