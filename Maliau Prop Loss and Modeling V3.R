# Load packages
library(XML)
library(dplyr)
library(stringr)
library(ggplot2)
library(seewave)
library(ggpubr)
library(stringr)
library(tidyverse)
library(geosphere)
library(lme4)

# Playback template table
SelectionIDsMaliau <- 
  read.delim("/Users/denaclink/Desktop/RStudio Projects/Propagation-Loss-2020-2021/SelectionLabels_S00974_20190811_101922_updated.txt")

# Remove replicate number
TempSoundType <- 
  str_split_fixed(SelectionIDsMaliau$Sound.Type, pattern = '_',n=3)[,2]

TempSoundType <- substr(TempSoundType,start = 1,stop=2)

# Remove pulses
PulsesToRemove <- which(TempSoundType!="Hf" & TempSoundType!="Ha")

PlaybackSeq <- seq(1,nrow(SelectionIDsMaliau),1)

PlaybackSeqUpdated <- PlaybackSeq[-PulsesToRemove]
SelectionIDsMaliau <- SelectionIDsMaliau[-PulsesToRemove,]


#MaliauDF <- read.csv("/Users/denaclink/Desktop/RStudio Projects/Propogation-Loss-2022/BackgroundNoiseRemovedMaliauAugust9AdaptiveNoiseMinNoise.csv")
MaliauDF <- read.csv("/Users/denaclink/Desktop/RStudio Projects/Propogation-Loss-2022/BackgroundNoiseRemovedMaliauFeb2023.csv")
#MaliauDF <- read.csv("/Users/denaclink/Desktop/RStudio Projects/Propogation-Loss-2022/BackgroundNoiseRemovedMaliauJuly2022.csv")
PredictedSpreading <- read.csv("Predicted_dB_Spherical.csv")
PredictedSpreadingMaliau <- subset(PredictedSpreading,Site=='Maliau')


MaliauDF <- na.omit(MaliauDF)
MaliauDF <- droplevels(subset(MaliauDF, date != '20190825'))

# Remove pulses
# Remove replicate number
TempSoundType <- 
  str_split_fixed(MaliauDF$Sound.Type, pattern = '_',n=3)[,2]

TempSoundType <- substr(TempSoundType,start = 1,stop=2)

PulsesToRemove <- which(TempSoundType!="Hf" & TempSoundType!="Ha")

MaliauDF <- MaliauDF[-PulsesToRemove,]


# Read in GPS data
source('readGPX.R')
recorder.gps <- readGPX("/Users/denaclink/Downloads/MB Playbacks 50 m.GPX") 


# Convert name so that it matches dataframe
recorder.gps$waypoints$name <- str_remove(recorder.gps$waypoints$name, '0')

# Subset only necessary columns from GPS data
small.gps.df <- recorder.gps$waypoints[,c('lon','lat','name')]
colnames(small.gps.df) <- c('lon','lat','recorder')

# Combine into a distance matrix
xy.coords <- cbind(c(small.gps.df$lon), 
                   c(small.gps.df$lat))

dist.mat <- distm( xy.coords, fun = distHaversine)

# Add recorder names to distance matrix
colnames(dist.mat) <- c(as.character(small.gps.df$recorder))
rownames(dist.mat) <- c(as.character(small.gps.df$recorder))


dist.to.playback.maliau <- 26.4    #17.1

# Check output
dist.source.vector <- ((dist.mat+dist.to.playback.maliau)[,1])

# Calculate hypotenuse distances
dist.source.vector <- sqrt(20^2 + dist.source.vector^2) 

#dist.source.vector[2:9] <- PredictedSpreadingMaliau$Dist2

# Create an index with unique date/time combinations
date.time.combo <- paste(MaliauDF$date,MaliauDF$time,sep='_')
unique.date.time.combo <- unique(date.time.combo)


# Create empty dataframe for propagation loss
observed.prop.lossMaliau <- data.frame()

# Loop to calculate propagation loss
for(z in 1:length(unique.date.time.combo)) { #tryCatch({ 
  
  # Subset by unique date/time index
  temp.date.time.subset <- 
    str_split_fixed(unique.date.time.combo[z],pattern = '_',n=2) 
  
  # Subset data frame to focus on unique date/time
  temp.playback <- subset(MaliauDF, date==temp.date.time.subset[,1] & time==temp.date.time.subset[,2])
  
  # See how many unique playbacks
  unique(temp.playback$file.name)
  
  # Create an index for each unique file in the playback
  file.index <- unique(temp.playback$file.name)
  SelectionIndex <- unique(temp.playback$Sound.Type)
  
  # This isolates each selection in the original template one by one
  for(a in 1:length(SelectionIndex)){
    
    # Subset the same selection from each of the recorders
    small.sample.playback.test <- data.frame()
    for(b in 1:length(file.index) ){
      temp.table <- subset(temp.playback,file.name==file.index[b])
      #temp.table$Sound.Type <- SelectionIDsMaliau$Sound.Type
      temp.table <- temp.table[a,]
      small.sample.playback.test <- rbind.data.frame(small.sample.playback.test,temp.table )
    }
    
    small.sample.playback.test <- na.omit(small.sample.playback.test[order(small.sample.playback.test$recorder),])
    
    # Create an index for each unique recorder in the new subset dataset
    recorder.index.test <- unique(small.sample.playback.test$recorder)
    
    # Create a new column with receive levels standardized so the closest recorder is 0
    small.sample.playback.test$PowerDb.zero <- 
      small.sample.playback.test$PowerDb-small.sample.playback.test$PowerDb[1]
    

    # Loop to calculate propagation loss; note the index starts at 2 since we use the closest one as the reference
    if(length(recorder.index.test)>0){
    if(recorder.index.test[1]=='M1'){
    for(c in 2:length(recorder.index.test)){tryCatch({ 
      print(c)
      # Isolate the recorder that we will use to estimate receive levels
      temp.recorder.received <- subset(small.sample.playback.test,recorder==recorder.index.test[c])
      
      # Isolate the recorder we consider as the 'source' for our relative calculations
      temp.recorder.source <- subset(small.sample.playback.test,recorder==recorder.index.test[1])
      
      # Based on our distance matrix above calculate the distance between the two recorders
      receive.dist <- dist.source.vector[c(temp.recorder.received$recorder)]
      
      source.dist <- dist.source.vector[c(temp.recorder.source$recorder)]
       
      distance.from.source <- receive.dist - source.dist
     
       # Assign the actual receive level (not zeroed) to new variable
      actual.receive.level <- temp.recorder.received$PowerDb
      
      # Assign zeroed receive level to new variable 
      zero.receive.level <- temp.recorder.received$PowerDb.zero
      
      # Assign 'source' level to new variable
      source.level <- temp.recorder.source$PowerDb.zero
      
      # Assign distance to new variable
      distance <- distance.from.source
      print(distance)
      isolate.distance <- which.min(abs(PredictedSpreadingMaliau$Dist2 - distance))
      
      PredictedSpreadingMaliauTemp <- PredictedSpreadingMaliau[isolate.distance,]
      
      ActualDbDifference <- temp.recorder.source$PowerDb  - temp.recorder.received$PowerDb  
      
      ExcessAttenuation <-  ActualDbDifference -PredictedSpreadingMaliauTemp$dBLoss_Spherical
      
      # Assign noise level estimate to new variable
      noise.level <- temp.recorder.received$NoisevalueDb
      
      # Calculate the distance ratio for propagation loss equation
      dist.ratio <- log10(receive.dist/source.dist)
      
      # Calculate the 'magic x'
      magic.x <-  zero.receive.level/dist.ratio
      
      # dB per doubling distance
      dBdoubledist <- magic.x*log10(10/5)
      print(dBdoubledist)
      
      # Assign sound type to new variable
      Sound.type <- SelectionIndex[a]
      
      # Assign time  to new variable
      time <- temp.recorder.received$time
      
      # Assign date to new variable
      date <- temp.recorder.received$date
      
      noise.level <- temp.recorder.received$NoisevalueDb
      
      # Combine all into a new temp dataframe
      temp.df <- cbind.data.frame(zero.receive.level,actual.receive.level,source.level,distance,Sound.type,time,date,magic.x,ExcessAttenuation,dBdoubledist,noise.level)
      print(temp.df)
      # Combine all observations into one data frame
      observed.prop.lossMaliau <- rbind.data.frame(observed.prop.lossMaliau,temp.df)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    #  DoublingDistanceDF <- rbind.data.frame(DoublingDistanceDF,DoublingDistanceDFtemp)
    }
    }
    }
  }
#}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

observed.prop.lossMaliauSubset<- observed.prop.lossMaliau
observed.prop.lossMaliauSubset<- droplevels(subset(observed.prop.lossMaliauSubset,date!='20190825'))

observed.prop.lossMaliauSubset$distance <- round(observed.prop.lossMaliauSubset$distance,0)


unique(observed.prop.lossMaliauSubset$distance)
#observed.prop.lossMaliauSubset <- read.csv('observed.prop.lossMaliauSubset.csv')
observed.prop.lossMaliauSubset$time <- as.factor(observed.prop.lossMaliauSubset$time)
observed.prop.lossMaliauSubset$Call.category <- str_split_fixed(observed.prop.lossMaliauSubset$Sound.type,pattern = '_',n=3)[,2]
# observed.prop.lossMaliauSubset <-droplevels( subset(observed.prop.lossMaliauSubset,
#                                          Call.category!="Pmor"))
observed.prop.lossMaliauSubset$distance <- round(observed.prop.lossMaliauSubset$distance,0)


# Prep Maliau Data

observed.prop.lossMaliauSubset$Species <- 
  recode(observed.prop.lossMaliauSubset$Call.category, Hfunstart = "NGreyGibbon",
         Hfuntrill = "NGreyGibbon",
         Halbstart='WhiteBeardGibbon',Halbend='WhiteBeardGibbon',Halbpeak='WhiteBeardGibbon')

# Add time category
observed.prop.lossMaliauSubset$time <- as.numeric(as.character(observed.prop.lossMaliauSubset$time))
 
observed.prop.lossMaliauSubset<-  observed.prop.lossMaliauSubset%>%
  mutate(TimeCat = case_when(
    time <= 700  ~ 'Dawn',
    time >= 800 & time <= 1200 ~ 'Morning',
    TRUE ~ 'Afternoon'
  ))


unique(observed.prop.lossMaliauSubset$TimeCat)

observed.prop.lossMaliauSubset$playback.num <- paste(observed.prop.lossMaliauSubset$date,observed.prop.lossMaliauSubset$time,sep='_')
observed.prop.lossMaliauSubset$habitat <- rep('D',nrow(observed.prop.lossMaliauSubset))
observed.prop.lossMaliauSubset$site <- rep('Maliau',nrow(observed.prop.lossMaliauSubset)) 
observed.prop.lossMaliauSubset$distance <- as.factor(observed.prop.lossMaliauSubset$distance)

observed.prop.lossMaliauSubset$distance <- as.numeric(as.character(observed.prop.lossMaliauSubset$distance))

# Remove outliers
hist(observed.prop.lossMaliauSubset$magic.x)

Call_dist <- unique(paste(observed.prop.lossMaliauSubset$Call.category,
      observed.prop.lossMaliauSubset$distance, sep='_'))

observed.prop.lossMaliauOutRM <- data.frame()
for(b in 1:length(Call_dist)){
  print(b)
  TempSubVals <- str_split_fixed(Call_dist[b],pattern = '_',n=2)
  TempDF <- subset(observed.prop.lossMaliauSubset,
                   Call.category==TempSubVals[1] & distance==TempSubVals[2])
  outliers <- boxplot.stats(TempDF$magic.x)$out
  IndexRM <- which(TempDF$magic.x %in% outliers)
  TempDF <- TempDF[-IndexRM,]

  observed.prop.lossMaliauOutRM <- rbind.data.frame(observed.prop.lossMaliauOutRM,TempDF)
}


hist(observed.prop.lossMaliauOutRM$magic.x)


observed.prop.lossMaliauOutRMdf  <- observed.prop.lossMaliauOutRM

# Remove this playback due to presence of calling gibbons
observed.prop.lossMaliauOutRMdf  <- 
  subset(observed.prop.lossMaliauOutRMdf, time!='720' )


round(as.numeric(as.character(observed.prop.lossMaliauOutRMdf$time)),0)

table(observed.prop.lossMaliauOutRMdf$time,observed.prop.lossMaliauOutRMdf$date)
hist(observed.prop.lossMaliauOutRMdf$magic.x)

observed.prop.lossMaliauOutRMdf$log.distance <- log10(observed.prop.lossMaliauOutRMdf$distance)

unique(observed.prop.lossMaliauOutRMdf$time)                                                                                  

ggboxplot(data=observed.prop.lossMaliauOutRMdf,
          y='actual.receive.level', x='Call.category')

ggboxplot(data=observed.prop.lossMaliauOutRMdf,
          y='magic.x', x='Call.category')


observed.prop.lossMaliauOutRMdf$date <- as.factor(observed.prop.lossMaliauOutRMdf$date)
observed.prop.lossMaliauOutRMdf$TimeCat <- as.factor(observed.prop.lossMaliauOutRMdf$TimeCat)
observed.prop.lossMaliauOutRMdf$Call.category <- as.factor(observed.prop.lossMaliauOutRMdf$Call.category )



observed.prop.lossMaliauOutRMdf$time <- as.factor(observed.prop.lossMaliauOutRMdf$time)
hist((observed.prop.lossMaliauOutRMdf$dBdoubledist))
hist((observed.prop.lossMaliauOutRMdf$magic.x))

unique(observed.prop.lossMaliauOutRMdf$time)

gghistogram(data=observed.prop.lossMaliauOutRMdf,
            x='dBdoubledist',facet.by = 'time',fill='date')

gghistogram(data=observed.prop.lossMaliauOutRMdf,
            x="noise.level",facet.by = 'time',fill='Species')

observed.prop.lossMaliauOutRMdf <-subset(observed.prop.lossMaliauOutRMdf,
       dBdoubledist <0 & dBdoubledist >-20)


observed.prop.lossMaliauOutRMdf$time <-
  revalue(observed.prop.lossMaliauOutRMdf$time, c('640' = "600", '840' = "800",
                                                  '920' = "900",
                                                  '1120' = "1100",
                                                '1040' = "1000",'1440' = "1400",
                                                '1320' = "1300",
                                                '1520' = "1500",
                                                '1640' = "1600"))


Maliau.lmm.prop.loss.null <- lmer(dBdoubledist ~  (1|date), data=observed.prop.lossMaliauOutRMdf) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Maliau.lmm.prop.loss.full <- lmer(dBdoubledist ~ log.distance  + Species + TimeCat + (1|date), data=observed.prop.lossMaliauOutRMdf) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Maliau.lmm.prop.loss.notime <- lmer(dBdoubledist ~log.distance  + Species  + (1|date), data=observed.prop.lossMaliauOutRMdf) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Maliau.lmm.prop.loss.nolog.distance <- lmer(dBdoubledist ~  Species + TimeCat + (1|date), data=observed.prop.lossMaliauOutRMdf) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Maliau.lmm.prop.loss.time.hour <- lmer(dBdoubledist ~ log.distance + Species + time + (1|date), data=observed.prop.lossMaliauOutRMdf) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Maliau.lmm.prop.loss.call.cat <- lmer(dBdoubledist ~ log.distance + Call.category  + TimeCat + (1|date), data=observed.prop.lossMaliauOutRMdf) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Maliau.lmm.prop.loss.distonly <- lmer(dBdoubledist ~ log.distance + (1|date), data=observed.prop.lossMaliauOutRMdf) # + (Call.Type|recorder.ID + Call.Type|recorder.location)


bbmle::AICctab(Maliau.lmm.prop.loss.null,Maliau.lmm.prop.loss.full,
               Maliau.lmm.prop.loss.notime,Maliau.lmm.prop.loss.nolog.distance,
              Maliau.lmm.prop.loss.time.hour,
               Maliau.lmm.prop.loss.distonly,
               weights=T)


hist(resid(Maliau.lmm.prop.loss.time.hour))

simulationOutputMaliau <- simulateResiduals(fittedModel = Maliau.lmm.prop.loss.time.hour, plot = F)
plotQQunif(simulationOutputMaliau)

MaliauPropLoss <- sjPlot::plot_model(Maliau.lmm.prop.loss.time.hour  ,intercept=F,sort.est = F)+
  ggtitle('Maliau Propagation loss')+ theme_bw()+geom_hline(yintercept = 0)+ scale_color_manual(values=c('black','red'))

MaliauPropLoss

sjPlot::plot_model(Maliau.lmm.prop.loss.time.hour,type='pred',
                   intercept=F,sort.est = TRUE)


PlotlistMaliau <-sjPlot::plot_model(Maliau.lmm.prop.loss.time.hour, type='pred')

# Marginal effects tells us how a dependent variable (outcome) changes when a specific independent variable (explanatory variable) changes.
MaliauMarginalEffectsDist <- PlotlistMaliau[[1]]+ylab('Propagation loss \n (dB decrease per doubling distance)')+
  ggtitle('')+xlab('Log Distance')+theme_bw()

MaliauMarginalEffectsSpecies <- PlotlistMaliau[[2]]+ylab('Propagation loss \n (dB decrease per doubling distance)')+
  ggtitle('')+theme_bw()

MaliauMarginalEffectsTime <- PlotlistMaliau[[3]]+ylab('Propagation loss \n (dB decrease per doubling distance)')+
  ggtitle('')+xlab('Local time')+theme_bw()


cowplot::plot_grid(MaliauMarginalEffectsDist,
                   MaliauMarginalEffectsSpecies,
                   MaliauMarginalEffectsTime,
                   labels=c('A','B','C'),label_x = 0.9)


noise.model<- lmer(noise.level ~ log.distance + Call.category  + time + (1|date), data=observed.prop.lossMaliauOutRMdf) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
noise.model.plot <-sjPlot::plot_model(noise.model  ,intercept=F,sort.est = F)+ggtitle('Maliau Noise loss')+ theme_bw()+geom_hline(yintercept = 0)
cowplot::plot_grid( MaliauPropLoss,noise.model.plot)

ggpubr::ggboxplot(data=observed.prop.lossMaliauOutRMdf,
                  fill='Species',y='noise.level',
                  x ='distance') +ylab('Receive levels')


ggpubr::ggboxplot(data=observed.prop.lossMaliauOutRMdf,
                  fill='Call.category',y='dBdoubledist',
                  x ='time') +ylab('Receive levels')




ggpubr::ggboxplot(data=observed.prop.lossMaliauOutRMdf,
                  fill='Species',y='dBdoubledist',#shape='TimeCat',
                  x ='Species') +ylab('dB loss per doubling distance')

ggpubr::gghistogram(data=observed.prop.lossMaliauOutRMdf,
                  facet.by  ='Call.category',x='dBdoubledist',fill='log.distance',
                  ) +ylab('dB loss per doubling distance')


Trillonly <- subset(observed.prop.lossMaliauOutRMdf,Call.category=='Hfuntrill')
ggpubr::gghistogram(data=Trillonly,
                    facet.by  ='distance',x='dBdoubledist',
) +ylab('Num Obs')



observed.prop.lossRunganSubset$Species <- factor(observed.prop.lossRunganSubset$Species, levels = c("NGreyGibbon", "WhiteBeardGibbon", "OrangCKalimantan",
                                                                                                    "OrangSabah"))

ggpubr::ggboxplot(data=observed.prop.lossRunganSubset,
                  fill='Species',y='dBdoubledist',#shape='TimeCat',
                  x ='habitat',
                  facet.by = 'Species') +ylab('Sound Transmission \n (dB loss per doubling distance)')+
   xlab('Habitat type')+ guides(fill=F)

# Combine into a single plot
CombinedDF <- rbind.data.frame(observed.prop.lossRunganSubset[,c('magic.x',"dBdoubledist","habitat", "Species",'dBdoubledist',"noise.level","time","TimeCat" )],
observed.prop.lossMaliauOutRMdf[,c('magic.x',"dBdoubledist","habitat", "Species",'dBdoubledist',"noise.level","time" ,"TimeCat" )])


ggpubr::ggboxplot(data=CombinedDF,
                  fill='Species',y='dBdoubledist',#shape='TimeCat',
                  x ='habitat',
                  facet.by = 'Species',outlier.shape = NA) +ylab('Sound Transmission \n (dB loss per doubling distance)')+
  xlab('Habitat type')+ guides(fill=F)+ geom_hline(yintercept=-6, linetype="dashed", color='black')

ggpubr::ggboxplot(data=CombinedDF,
                  fill='Species',y='magic.x',#shape='TimeCat',
                  x ='habitat',
                  facet.by = 'Species') +ylab('Propagation loss')+
  xlab('Habitat type')+ guides(fill=F)+ geom_hline(yintercept=-20, linetype="dashed", color='grey')

# Background noise --------------------------------------------------------
ggboxplot(data=CombinedDF,x='Species',y="noise.level",fill='habitat',
          outlier.shape = NA)


