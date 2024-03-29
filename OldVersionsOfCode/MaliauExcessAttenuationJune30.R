library(XML)
library(dplyr)
library(stringr)
library(ggplot2)
library(seewave)
library(ggpubr)
library(stringr)
library(tidyverse)
library(geosphere)

# Playback template table
SelectionIDsMaliau <- 
  read.delim("/Users/denaclink/Desktop/RStudio Projects/Propagation-Loss-2020-2021/SelectionLabels_S00974_20190811_101922_updated.txt")

TempSoundType <- 
  str_split_fixed(SelectionIDsMaliau$Sound.Type, pattern = '_',n=3)[,2]

TempSoundType <- substr(TempSoundType,start = 1,stop=2)

# Remove pulses
PulsesToRemove <- which(TempSoundType!="Hf" & TempSoundType!="Ha")

PlaybackSeq <- seq(1,nrow(SelectionIDsMaliau),1)

PlaybackSeqUpdated <- PlaybackSeq[-PulsesToRemove]
SelectionIDsMaliau <- SelectionIDsMaliau[-PulsesToRemove,]


MaliauDF <- read.csv('/Users/denaclink/Desktop/RStudio Projects/Propagation-Loss-2020-2021/BackgroundNoiseRemovedMaliauJuly2022add20times.csv')

PredictedSpreading <- read.csv("/Users/denaclink/Desktop/RStudio Projects/Propagation-Loss-2020-2021/Predicted_dB_Spherical.csv")
PredictedSpreadingMaliau <- subset(PredictedSpreading,Site=='Maliau')

head(MaliauDF)
table(MaliauDF$date)
unique(MaliauDF$time)

MaliauDF <-
  subset(MaliauDF,time==600 |time==640 |time==800|time==840|time==1000|time==1040|
         time==1400 |time==1440 |time==1600|time==1640)

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


dist.to.playback.maliau <- 26.4

# Check output
dist.source.vector <- ((dist.mat+dist.to.playback.maliau)[,1])


# Create an index with unique date/time combinations
date.time.combo <- paste(MaliauDF$date,MaliauDF$time,sep='_')
unique.date.time.combo <- unique(date.time.combo)


# Create empty dataframe for propagation loss
observed.prop.lossMaliau <- data.frame()

# Loop to calculate propagation loss
for(z in 1:length(unique.date.time.combo)) { tryCatch({ 
  
  # Subset by unique date/time index
  temp.date.time.subset <- 
    str_split_fixed(unique.date.time.combo[z],pattern = '_',n=2) 
  
  # Subset data frame to focus on unique date/time
  temp.playback <- subset(MaliauDF, date==temp.date.time.subset[,1] & time==temp.date.time.subset[,2])
  
  # See how many unique playbacks
  unique(temp.playback$file.name)
  
  # Create an index for each unique file in the playback
  file.index <- unique(temp.playback$file.name)
  SelectionIndex <- (SelectionIDsMaliau$Sound.Type)
  
  # This isolates each selection in the original template one by one
  for(a in 1:length(SelectionIndex)){
    
    # Subset the same selection from each of the recorders
    small.sample.playback.test <- data.frame()
    for(b in 1:length(file.index) ){
      temp.table <- subset(temp.playback,file.name==file.index[b])
      temp.table$Sound.Type <- SelectionIDsMaliau$Sound.Type
      temp.table <- temp.table[a,]
      small.sample.playback.test <- rbind.data.frame(small.sample.playback.test,temp.table )
    }
    
    # Create an index for each unique recorder in the new subset dataset
    recorder.index.test <- unique(small.sample.playback.test$recorder)
    
    # Create a new column with receive levels standardized so the closest recorder is 0
    small.sample.playback.test$PowerDb.zero <- 
      small.sample.playback.test$PowerDb-small.sample.playback.test$PowerDb[1]
    

  
     
    
    # Loop to calculate propagation loss; note the index starts at 2 since we use the closest one as the reference
    for(c in 2:length(recorder.index.test)){tryCatch({ 
      print(c)
      # Isolate the recorder that we will use to estimate receive levels
      temp.recorder.received <- subset(small.sample.playback.test,recorder==recorder.index.test[c])
      
      # Isolate the recorder we consider as the 'source' for our relative calculations
      temp.recorder.source <- subset(small.sample.playback.test,recorder==recorder.index.test[1])
      
      # Based on our distance matrix above calculate the distance between the two recorders
      receive.dist <- dist.source.vector[c(temp.recorder.received$recorder)]
      
      source.dist <- dist.source.vector[1]
       
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
      dist.ratio <- log10(distance/dist.to.playback.maliau)
      
      # Calculate the 'magic x'
      magic.x <-  zero.receive.level/dist.ratio
      
      # dB per doubling distance
      dBdoubledist <- magic.x*log10(10/5)
      print(dBdoubledist)
      
      # Assign sound type to new variable
      Sound.type <- temp.recorder.received$Sound.Type
      
      # Assign time  to new variable
      time <- temp.recorder.received$time
      
      # Assign date to new variable
      date <- temp.recorder.received$date
      
      # Combine all into a new temp dataframe
      temp.df <- cbind.data.frame(zero.receive.level,actual.receive.level,source.level,distance,Sound.type,time,date,magic.x,noise.level,ExcessAttenuation,dBdoubledist)
      
      # Combine all observations into one data frame
      observed.prop.lossMaliau <- rbind.data.frame(observed.prop.lossMaliau,temp.df)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    #  DoublingDistanceDF <- rbind.data.frame(DoublingDistanceDF,DoublingDistanceDFtemp)
    }
    
  }
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

write.csv(observed.prop.lossMaliau,'observed.prop.lossMaliau.csv',row.names = F)

observed.prop.lossMaliau <- read.csv('observed.prop.lossMaliau.csv')
observed.prop.lossMaliau$ExcessAttenuation <- round(as.numeric(observed.prop.lossMaliau$ExcessAttenuation),2)
observed.prop.lossMaliau$time <- as.factor(observed.prop.lossMaliau$time)
observed.prop.lossMaliau$Call.category <- str_split_fixed(observed.prop.lossMaliau$Sound.type,pattern = '_',n=3)[,2]
observed.prop.lossMaliau$distance <- round(observed.prop.lossMaliau$distance,0)
observed.prop.lossMaliauSubset <- observed.prop.lossMaliau #subset(observed.prop.lossMaliau,Call.category=="Pmor" | Call.category=="Hfuntrill" |Call.category== "Hfunstart")

#observed.prop.lossMaliauSubset <- droplevels(subset(observed.prop.lossMaliauSubset,ExcessAttenuation>0))
ggboxplot(data=observed.prop.lossMaliauSubset, x='distance',
            y='ExcessAttenuation',facet.by = 'time')#+ylim(0,500)


ggboxplot(data=observed.prop.lossMaliauSubset, x='distance',
          y='magic.x',color = 'time',facet.by = 'Call.category')#+ylim(0,500)

ggboxplot(data=observed.prop.lossMaliauSubset, y='dBdoubledist',
          x='distance',fill = 'Call.category')#+ylim(0,500)

MaliauPlot <- ggboxplot(data=observed.prop.lossMaliauSubset, y='magic.x',
          x='distance',fill = 'Call.category')#+ylim(0,500)


observed.prop.lossMaliauSubsetPmor <- subset(observed.prop.lossMaliau,Call.category=="Pmor")
ggboxplot(data=observed.prop.lossMaliauSubsetPmor, x='distance',
          y='ExcessAttenuation',fill = 'time',facet.by = 'Sound.type')#+ylim(0,500)

#subset(observed.prop.lossMaliau,distance== 0)
MaliauPlot <- ggscatter(data=observed.prop.lossMaliauSubset, y='distance',
                        x='ExcessAttenuation',color = 'Sound.type',facet.by = 'Sound.type')+ylim(0,500)+ geom_jitter(width = 1.5, height = 1)

MaliauPlot

# Median values -----------------------------------------------------------
# Take the median value for replicate selections
RecorderDateTime.index <- unique(paste(observed.prop.lossMaliau$distance,
                                       observed.prop.lossMaliau$date,
                                       observed.prop.lossMaliau$time,sep='_'))

observed.prop.lossMaliau$RecorderDateTime <- unique(paste(#observed.prop.lossMaliau$distance,
  observed.prop.lossMaliau$date,
  observed.prop.lossMaliau$time,sep='_'))

observed.prop.lossMaliau$sound.type.for.median <- str_split_fixed (observed.prop.lossMaliau$Sound.type, pattern = '_',n=2)[,2]

# Isolate sound category names
observed.prop.lossMaliau$Sound.category <- 
  str_split_fixed(observed.prop.lossMaliau$Sound.type, pattern = '_',n=3)[,2]


observed.prop.lossMaliau.median.df <- data.frame()

for(a in 1:length(RecorderDateTime.index)){
  temp.observed <-  subset(observed.prop.lossMaliau,RecorderDateTime==RecorderDateTime[a])
  
  sound.type.index <- 
    unique(temp.observed$sound.type.for.median  )
  
  distance.index <- unique(temp.observed$distance  )
  
  for(b in 1:length(sound.type.index)){
    
    temp.subset <- subset(temp.observed,sound.type.for.median==sound.type.index[b])
    
    for(c in 1:length(distance.index)){
      
      temp.subset.distance <-  subset(temp.subset,distance==distance.index[c])
      median.x <- median(temp.subset.distance$magic.x) 
      
      temp.subset.distance <- temp.subset.distance[1,c("distance", "actual.receive.level","Sound.type", "time", "date", "magic.x", 
                                                       "RecorderDateTime", "sound.type.for.median", "Sound.category")]
      
      temp.subset.distance$median.x <- median.x
      observed.prop.lossMaliau.median.df <- rbind.data.frame(observed.prop.lossMaliau.median.df,temp.subset.distance)
    }
  }
  
}

observed.prop.lossMaliau.median.df <- na.omit(observed.prop.lossMaliau.median.df)

# Subset so we focus on primates
observed.prop.lossMaliau.subset <-  subset(observed.prop.lossMaliau.median.df,
                                           Sound.category=="Hfunstart" |Sound.category=="Hfuntrill" 
                                           |Sound.category=="Halbpeak"
                                           |Sound.category=="Halbend"
                                           |Sound.category=="Pmor"|Sound.category=="PwurP"|Sound.category=="PwurS") #

# Give more informative names
observed.prop.lossMaliau.subset$Sound.category <- 
  plyr::revalue(observed.prop.lossMaliau.subset$Sound.category,
                c(Halbpeak= "Gibbon C. Kali Peak",Halbend= "Gibbon C. Kali End",
                  Hfunstart = "Gibbon Sabah Intro",
                  Hfuntrill = "Gibbon Sabah Trill",Pmor="Orangutan Sabah", 
                  PwurP="Orangutan C. Kali Pulse",
                  PwurS="Orangutan C. Kali Sigh"))

# Plot observed change in receive level by distance
MaliauPlot <- ggpubr::ggscatter(data = observed.prop.lossMaliau.subset,
                                x='distance', y='ExcessAttenuation',
                                color='time',
                                facet.by = 'Sound.category',xlim=c(0,400), #ylim=c(-70,-20),
                                add = c("loess"),title='')+
  xlab('Distance (m)')+ ylab('Receive Level (dB)')+ theme(legend.position = "none")+ggtitle('Maliau')

MaliauPlot
# Plot observed change in receive level by distance for all signals
ggpubr::ggscatter(data = observed.prop.lossMaliauSubset,x='distance', y='actual.receive.level',
                  color='time',
                  facet.by = 'Call.category',#xlim=c(0,300), #ylim=c(-70,-20),
                  add = c("loess"),title='')+
  xlab('Distance (m)')+ ylab('Receive Level (dB)')+ theme(legend.position = "none")

ggpubr::ggscatter(data = observed.prop.lossMaliauSubset,x='distance', y='magic.x',
                  color='date',
                  facet.by = 'Call.category',#xlim=c(0,300), #ylim=c(-70,-20),
                  add = c("loess"),title='')+
  xlab('Distance (m)')+ ylab('Magic x')#+ theme(legend.position = "none")

observed.prop.lossMaliauSubset[which(observed.prop.lossMaliauSubset$magic.x >0),]

head(observed.prop.lossMaliauSubset)
hist(observed.prop.lossMaliauSubset$magic.x)
