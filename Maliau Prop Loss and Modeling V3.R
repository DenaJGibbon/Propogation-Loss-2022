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

TempSoundType <- 
  str_split_fixed(SelectionIDsMaliau$Sound.Type, pattern = '_',n=3)[,2]

TempSoundType <- substr(TempSoundType,start = 1,stop=2)

# Remove pulses
PulsesToRemove <- which(TempSoundType!="Hf" & TempSoundType!="Ha")

PlaybackSeq <- seq(1,nrow(SelectionIDsMaliau),1)

PlaybackSeqUpdated <- PlaybackSeq[-PulsesToRemove]
SelectionIDsMaliau <- SelectionIDsMaliau[-PulsesToRemove,]


#MaliauDF <- read.csv("/Users/denaclink/Desktop/RStudio Projects/Propogation-Loss-2022/BackgroundNoiseRemovedMaliauAugust9AdaptiveNoiseMinNoise.csv")
MaliauDF <- read.csv("/Users/denaclink/Desktop/RStudio Projects/Propogation-Loss-2022/BackgroundNoiseRemovedMaliauAugust24SubtractMoreNoise.csv")
#MaliauDF <- read.csv("/Users/denaclink/Desktop/RStudio Projects/Propogation-Loss-2022/BackgroundNoiseRemovedMaliauJuly2022.csv")
PredictedSpreading <- read.csv("Predicted_dB_Spherical.csv")
PredictedSpreadingMaliau <- subset(PredictedSpreading,Site=='Maliau')

head(MaliauDF)
table(MaliauDF$date)
unique(MaliauDF$time)
nrow(MaliauDF)

MaliauDF <- na.omit(MaliauDF)

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
      
      # Combine all into a new temp dataframe
      temp.df <- cbind.data.frame(zero.receive.level,actual.receive.level,source.level,distance,Sound.type,time,date,magic.x,ExcessAttenuation,dBdoubledist)
      print(temp.df)
      # Combine all observations into one data frame
      observed.prop.lossMaliau <- rbind.data.frame(observed.prop.lossMaliau,temp.df)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    #  DoublingDistanceDF <- rbind.data.frame(DoublingDistanceDF,DoublingDistanceDFtemp)
    }
    
  }
#}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}



#observed.prop.lossMaliau <- read.csv('observed.prop.lossMaliau.csv')
observed.prop.lossMaliau$ExcessAttenuation <- round(as.numeric(observed.prop.lossMaliau$ExcessAttenuation),2)
observed.prop.lossMaliau$time <- as.factor(observed.prop.lossMaliau$time)
observed.prop.lossMaliau$Call.category <- str_split_fixed(observed.prop.lossMaliau$Sound.type,pattern = '_',n=3)[,2]
observed.prop.lossMaliau$distance <- round(observed.prop.lossMaliau$distance,0)
observed.prop.lossMaliauSubset <- observed.prop.lossMaliau #subset(observed.prop.lossMaliau,Call.category=="Pmor" | Call.category=="Hfuntrill" |Call.category== "Hfunstart")

#write.csv(observed.prop.lossMaliau,'observed.prop.lossMaliauAugust14adaptive.csv',row.names = F)

# Prep Maliau Data

observed.prop.lossMaliauSubset$Species <- 
  recode(observed.prop.lossMaliauSubset$Call.category, Hfunstart = "NGreyGibbon",
         Hfuntrill = "NGreyGibbon",
         Halbstart='WhiteBeardGibbon',Halbend='WhiteBeardGibbon',Halbpeak='WhiteBeardGibbon')

# Add time category
observed.prop.lossMaliauSubset$time <- as.numeric(as.character(observed.prop.lossMaliauSubset$time))

observed.prop.lossMaliauSubset <-  observed.prop.lossMaliauSubset %>%
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
# magic.x.mean <- mean(observed.prop.lossMaliauSubset$magic.x)
#  magic.x.sd <- sd((observed.prop.lossMaliauSubset$magic.x))
#  magic.x.sd <- magic.x.sd*3
# 
outliers <- boxplot.stats(observed.prop.lossMaliauSubset$magic.x)$out
# 
out_ind <- which(observed.prop.lossMaliauSubset$magic.x %in% c(outliers))
out_ind
# 
# # Subset removing outliers
observed.prop.lossMaliauSubset <- observed.prop.lossMaliauSubset[-out_ind,]

#observed.prop.lossMaliauSubset <- subset(observed.prop.lossMaliauSubset, magic.x > -60 & magic.x < -15)

observed.prop.lossMaliauSubset <- subset(observed.prop.lossMaliauSubset,
                                             distance <=315)

# Remove this playback due to presence of calling gibbons
observed.prop.lossMaliauSubset <- subset(observed.prop.lossMaliauSubset, time!='720')

hist(observed.prop.lossMaliauSubset$magic.x)

observed.prop.lossMaliauSubset$distance <- log10(observed.prop.lossMaliauSubset$distance)

unique(observed.prop.lossMaliauSubset$time)                                                                                  

ggboxplot(data=observed.prop.lossMaliauSubset,
          y='actual.receive.level', x='Call.category')

ggboxplot(data=observed.prop.lossMaliauSubset,
          y='magic.x', x='Call.category')

#observed.prop.lossMaliauSubset <- droplevels(subset(observed.prop.lossMaliauSubset, time!=920 & time!=1120  &time!=1320 & time!=1520))
observed.prop.lossMaliauSubset$date <- as.factor(observed.prop.lossMaliauSubset$date)
observed.prop.lossMaliauSubset$TimeCat <- as.factor(observed.prop.lossMaliauSubset$TimeCat)
observed.prop.lossMaliauSubset$Call.category <- as.factor(observed.prop.lossMaliauSubset$Call.category )

levels(observed.prop.lossMaliauSubset$Call.category)

Maliau.lmm.prop.loss.null <- lmer(magic.x ~  (1|date), data=observed.prop.lossMaliauSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Maliau.lmm.prop.loss.full <- lmer(magic.x ~ distance  + Species + TimeCat + (1|date), data=observed.prop.lossMaliauSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Maliau.lmm.prop.loss.notime <- lmer(magic.x ~ distance  + Species  + (1|date), data=observed.prop.lossMaliauSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Maliau.lmm.prop.loss.nodistance <- lmer(magic.x ~  Species + TimeCat + (1|date), data=observed.prop.lossMaliauSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Maliau.lmm.prop.loss.full.cat <- lmer(magic.x ~ distance*Species + TimeCat + (1|date), data=observed.prop.lossMaliauSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Maliau.lmm.prop.loss.full.species <- lmer(magic.x ~ distance  + Species + TimeCat + (1|date), data=observed.prop.lossMaliauSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Maliau.lmm.prop.loss.full.species <- lmer(magic.x ~ distance  + Species + TimeCat + (1|date), data=observed.prop.lossMaliauSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)


bbmle::AICctab(Maliau.lmm.prop.loss.null,Maliau.lmm.prop.loss.full,Maliau.lmm.prop.loss.full.species,
               Maliau.lmm.prop.loss.notime,Maliau.lmm.prop.loss.nodistance, weights=T)

summary(Maliau.lmm.prop.loss.full)
hist(resid(Maliau.lmm.prop.loss.full))
MaliauPropLoss <- sjPlot::plot_model(Maliau.lmm.prop.loss.full.species ,intercept=F,sort.est = TRUE)+ggtitle('Maliau Propogation loss')+ theme_bw()+geom_hline(yintercept = 0)
MaliauPropLoss

sjPlot::plot_model(Maliau.lmm.prop.loss.full,type='eff',intercept=F,sort.est = TRUE)



ggpubr::ggboxplot(data=observed.prop.lossMaliauSubset,
                  fill='Call.category',y='actual.receive.level',
                  x ='distance') +ylab('Receive levels')



ggpubr::ggboxplot(data=observed.prop.lossMaliauSubset,
                  fill='Call.category',y='magic.x',
                  x ='distance') +ylab('Magic x')

