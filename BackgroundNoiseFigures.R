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

#standard error function
se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))


RunganNoise <-
  read.csv("ThirdOctaveBandDFRunganRungan.csv")

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
  summarise(ambient.noise.high=meandB(noise.valuedb)+se(noise.valuedb),
            ambient.noise.low=meandB(noise.valuedb)-se(noise.valuedb),
            ambient.noise=meandB(noise.valuedb))

# Drop missing habitat
RunganNoiseMedian <- subset(RunganNoiseMedian,habitat=='K'| habitat=='LP' | habitat=='MS')

RunganNoiseMedian$time <- as.factor(RunganNoiseMedian$time)
RunganNoiseMedian$time <- as.factor(RunganNoiseMedian$time)

ggplot(RunganNoiseMedian,aes(center.freq,ambient.noise,
                       group=time,colour=time,linetype=time))+
  stat_summary(data=RunganNoiseMedian,fun.y=meandB,geom="line",alpha=0.2,aes(group=time))+
  stat_summary(data=RunganNoiseMedian,fun.y=meandB,geom="line",aes(group=time))+
  geom_ribbon(data=RunganNoiseMedian,aes(ymin=ambient.noise.low,ymax=ambient.noise.high
                                         ,fill=time,group=time,color=NULL),alpha=0.5)+ 
  scale_fill_manual(values= matlab::jet.colors(length(unique(RunganNoiseMedian$time))) )+ theme_bw()+
  scale_color_manual(values= matlab::jet.colors(length(unique(MaliauNoise$time))) )+ theme_bw()#+ylim(25,70)
#+ylim(25,70)

table(RunganNoiseMedian$time)

MaliauNoise <-
  read.csv("ThirdOctaveBandDFMaliau.csv")

MaliauNoise$time <- str_split_fixed(MaliauNoise$wav.file,pattern='_',n=3)[,3]
MaliauNoise$time <- as.numeric(substr(MaliauNoise$time,1,2))


MaliauNoiseMedian <- MaliauNoise %>%
  group_by(time,center.freq) %>%
  summarise(ambient.noise.high=meandB(noise.valuedb)+se(noise.valuedb),
            ambient.noise.low=meandB(noise.valuedb)-se(noise.valuedb),
            ambient.noise=meandB(noise.valuedb))

MaliauNoise$time <- as.factor(MaliauNoise$time)
MaliauNoiseMedian$time <- as.factor(MaliauNoiseMedian$time)

MaliauNoisePlot <- ggplot(MaliauNoiseMedian,aes(center.freq,ambient.noise,
                          group=time,colour=time,linetype=time))+
  stat_summary(data=MaliauNoiseMedian,fun.y=meandB,geom="line",alpha=0.2,aes(group=time))+
  stat_summary(data=MaliauNoiseMedian,fun.y=meandB,geom="line",aes(group=time))+
  geom_ribbon(data=MaliauNoiseMedian,aes(ymin=ambient.noise.low,ymax=ambient.noise.high
                                  ,fill=time,group=time,color=NULL),alpha=0.5)+ 
  scale_fill_manual(values= matlab::jet.colors(length(unique(MaliauNoise$time))) )+ 
  scale_color_manual(values= matlab::jet.colors(length(unique(MaliauNoise$time))) )+ theme_bw()#+ylim(25,70)

cowplot::plot_grid(RunganNoisePlot,MaliauNoisePlot)

## Just do by habitat type

MaliauNoise$habitat <- rep('D',nrow(MaliauNoise))

CombinedNoise <- 
  rbind.data.frame(RunganNoiseAddHabitat,MaliauNoise)

CombinedNoiseMedian <- CombinedNoise %>%
  group_by(habitat,center.freq) %>%
  summarise(ambient.noise.high=meandB(noise.valuedb)+se(noise.valuedb),
            ambient.noise.low=meandB(noise.valuedb)-se(noise.valuedb),
            ambient.noise=meandB(noise.valuedb))



# Drop missing habitat
CombinedNoiseMedian <- subset(CombinedNoiseMedian,habitat=='D'|habitat=='K'| habitat=='LP' | habitat=='MS')


ggplot(CombinedNoiseMedian,aes(center.freq,ambient.noise,
                                                group=habitat,colour=habitat,linetype=habitat))+
  stat_summary(data=CombinedNoiseMedian,fun.y=meandB,geom="line",alpha=0.2,aes(group=habitat))+
  stat_summary(data=CombinedNoiseMedian,fun.y=meandB,geom="line",aes(group=habitat))+
  geom_ribbon(data=CombinedNoiseMedian,aes(ymin=ambient.noise.low,ymax=ambient.noise.high
                                         ,fill=habitat,group=habitat,color=NULL),alpha=0.5)+ 
  scale_fill_manual(values= matlab::jet.colors(length(unique(CombinedNoise$habitat))) )+ 
  scale_color_manual(values= matlab::jet.colors(length(unique(CombinedNoise$habitat))) )+ 
  theme_bw()+ylab(expression(paste('Ambient sound level (dB re 20', mu,'Pa)',sep=' ')))+xlab('Center Frequency (Hz)') #+ylim(25,70)

ggplot(CombinedNoiseMedian,aes(center.freq,ambient.noise,
                               group=habitat,colour=habitat,linetype=habitat))+
  stat_summary(data=CombinedNoiseMedian,fun.y=meandB,geom="line",alpha=0.2,aes(group=habitat))+
  stat_summary(data=CombinedNoiseMedian,fun.y=meandB,geom="line",aes(group=habitat))+
  geom_ribbon(data=CombinedNoiseMedian,aes(ymin=ambient.noise.low,ymax=ambient.noise.high
                                           ,fill=habitat,group=habitat,color=NULL),alpha=0.5)+ 
  scale_fill_manual(values= matlab::jet.colors(length(unique(CombinedNoise$habitat))) )+ 
  scale_color_manual(values= matlab::jet.colors(length(unique(CombinedNoise$habitat))) )+ 
  theme_bw()+ylab(expression(paste('Ambient sound level (dB re 20', mu,'Pa)',sep=' ')))+
  xlab('Center Frequency (Hz)')+xlim(0,2000)+ylim(25,45)

      
                  

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
