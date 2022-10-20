library(seewave)
library(tuneR)
library(stringr)
library(plyr)
library(ggpubr)
library(geosphere)
library(plotKML)
library(lme4)
library(lmerTest)
library(ggplot2)
library(MASS)
library(multcomp)
library(FSA)
library(DHARMa)
library(ggpubr)
library(dplyr)

dist.to.playback <- 10

SelectionIDsRungan <- 
  read.delim("/Users/denaclink/Desktop/RStudio Projects/Propagation-Loss-2020-2021/SelectionLabels_S00974_20190811_101922_updated_fundamental.txt")

TempSoundType <- 
  str_split_fixed(SelectionIDsRungan$Sound.Type, pattern = '_',n=3)[,2]

TempSoundType <- substr(TempSoundType,start = 1,stop=2)

# Remove pulses
PulsesToRemove <- which(TempSoundType!="Hf" & TempSoundType!="Ha"
                        & TempSoundType!="Pm" & TempSoundType!="Pw" )

PlaybackSeq <- seq(1,nrow(SelectionIDsRungan),1)

PlaybackSeqUpdated <- PlaybackSeq[-PulsesToRemove]
SelectionIDsRungan <- SelectionIDsRungan[-PulsesToRemove,]

#NOTE that there are two Pwur call types
RunganDF <- read.csv("BackgroundNoiseRemovedDFRunganAugust102022adaptiveCombinedDF.csv")
PredictedSpreading <- read.csv("/Users/denaclink/Desktop/RStudio Projects/Propagation-Loss-2020-2021/Predicted_dB_Spherical.csv")
PredictedSpreadingRungan <- subset(PredictedSpreading,Site=='Munkgu')

# Read in data file
# Each row corresponds to a playback
rungan_data <- read.csv("PropLoss_test_7Jul22.csv")

# Create an index with unique date/time combinations
date.time.combo <- paste(RunganDF$date,RunganDF$time,sep='_')
unique.date.time.combo <- unique(date.time.combo)

Loc_Name.index <- 
  unique(RunganDF$Loc_Name)

Loc_Name.index  <- Loc_Name.index[- which(Loc_Name.index %in% c('char1','char2','char3'))]

# Create dataframe to match characterization units  
char.matching <- data.frame(
  char = c('char1','char2','char3'),
  rec=c('S00976', 'S01143', 'S00974')
)


# Create empty dataframe for propagation loss
observed.prop.lossRungan <- data.frame()

# Loop to calculate propagation loss
for(z in 1:length(Loc_Name.index)) { #tryCatch({ 
  
  # Subset data frame to focus on unique date/time
  temp.playback <- subset(RunganDF,Loc_Name==Loc_Name.index[z])
  
  
  # Create an index for each unique file in the playback
  SelectionIndex <- (SelectionIDsRungan$Sound.Type)
  
  playback.index <-  which(rungan_data$Loc_Name == unique(temp.playback$Loc_Name))
  
  temp.playback.data <- rungan_data[playback.index[1],]
  
  temp.char.info <- char.matching[which(char.matching$rec== temp.playback.data$ARU_ID),]
  
  TempReferenceDF <- subset(RunganDF,Loc_Name== temp.char.info$char)
  
  # Create an index for each unique file in the playback
  distance.index <- unique(temp.playback$distance.from.source)
  SelectionIndex <- (SelectionIDsRungan$Sound.Type)
  
  temp.playback$distance.from.source <- as.character(temp.playback$distance.from.source)
  
  
  # This isolates each selection in the original template one by one
  for(a in 1:length(SelectionIndex)){
    
    # Subset the same selection from each of the recorders
    small.sample.playback.test <- data.frame()
    for(b in 1:length(distance.index) ){
      temp.table <- subset(temp.playback,distance.from.source==distance.index[b])
      temp.table$Sound.Type <- SelectionIDsRungan$Sound.Type
      temp.table <- temp.table[a,]
      small.sample.playback.test <- rbind.data.frame(small.sample.playback.test,temp.table )
    }
    
    
    small.sample.playback.test <- rbind.data.frame(small.sample.playback.test,
                                                   TempReferenceDF[a,])
    
    # Reorder based on distance from speaker
    small.sample.playback.test <-  arrange(small.sample.playback.test, distance.from.source)  
    
    # Create a new column with receive levels standardized so the closest recorder is 0
    small.sample.playback.test$PowerDb.zero <- 
      small.sample.playback.test$PowerDb-small.sample.playback.test$PowerDb[1]
    
    
    distance.index.test <- unique(small.sample.playback.test$distance.from.source)
    
    # Loop to calculate propagation loss; note the index starts at 2 since we use the closest one as the reference
    for(c in 2:length(distance.index.test)){#tryCatch({ 
     
      # Isolate the recorder that we will use to estimate receive levels
      temp.recorder.received <- subset(small.sample.playback.test,distance.from.source==distance.index.test[c])
      
      # Isolate the recorder we consider as the 'source' for our relative calculations
      temp.recorder.source <- subset(small.sample.playback.test,distance.from.source==distance.index.test[1])
      
       
      # Assign the actual receive level (not zeroed) to new variable
      actual.receive.level <- temp.recorder.received$PowerDb
      if(length(actual.receive.level)==0){
        print( small.sample.playback.test$Sound.Type[1] )
      }
      
      # Assign zeroed receive level to new variable 
      zero.receive.level <- temp.recorder.received$PowerDb.zero
      
      # Assign 'source' level to new variable
      source.level <- temp.recorder.source$PowerDb.zero
      
      # Assign distance to new variable
      distance <- as.numeric(temp.recorder.received$distance.from.source)
      
      isolate.distance <- which.min(abs(PredictedSpreadingRungan$Dist2 - distance))
      
      PredictedSpreadingRunganTemp <- PredictedSpreadingRungan[isolate.distance,]
      
      ActualDbDifference <- temp.recorder.source$PowerDb  - temp.recorder.received$PowerDb  
      
      ExcessAttenuation <-  ActualDbDifference -PredictedSpreadingRunganTemp$dBLoss_Spherical
      
      # Assign noise level estimate to new variable
      noise.level <- temp.recorder.received$NoisevalueDb
      
      # Calculate the distance ratio for propagation loss equation
      dist.ratio <- log10(distance/dist.to.playback)
      
      # Calculate the 'magic x'
      magic.x <-  zero.receive.level /dist.ratio
      
      # dB per doubling distance
      dBdoubledist <- magic.x*log10(10/5)
      print(dBdoubledist)

      # Assign sound type to new variable
      Sound.type <- temp.recorder.received$Sound.Type
      
      # Assign time  to new variable
      time <- temp.recorder.received$time
      
      # Assign date to new variable
      date <- temp.recorder.received$date
      
      # Habitat type
      habitat <- temp.playback.data$Habitat
      playback.num <- temp.playback.data$PB_No
      Loc_Name <- temp.playback.data$Loc_Name
      ARU_ID <- temp.playback.data$ARU_ID
      # Combine all into a new temp dataframe
      temp.df <- cbind.data.frame(zero.receive.level,actual.receive.level,source.level,distance,Sound.type,time,date,magic.x,noise.level,ExcessAttenuation,habitat,playback.num,dBdoubledist, Loc_Name, ARU_ID)
      
      # Combine all observations into one dataframe
      observed.prop.lossRungan <- rbind.data.frame(observed.prop.lossRungan,temp.df)
      
   # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
       }
    
  }
#}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

unique(observed.prop.lossRungan$Loc_Name)
table(observed.prop.lossRungan$time)

# Remove PBs
observed.prop.lossRungan <- droplevels(subset(observed.prop.lossRungan,playback.num!=1))
observed.prop.lossRungan <- droplevels(subset(observed.prop.lossRungan,playback.num!=110))

observed.prop.lossRungan$Call.category <- str_split_fixed(observed.prop.lossRungan$Sound.type,pattern = '_',n=3)[,2]
observed.prop.lossRunganSubset <- subset(observed.prop.lossRungan,Sound.type=="1_Pmor_P2" |Sound.type== "1_Hfunstart")
ggscatter(data=observed.prop.lossRunganSubset, y='distance',
          x='magic.x',color = 'Call.category',facet.by = 'habitat')#+geom_jitter(width = 1.5, height = 1)


# Modeling ----------------------------------------------------------------


# Prep Rungan Data
observed.prop.lossRunganSubset <- read.csv('observed.prop.lossRunganAugust1.csv')



observed.prop.lossRunganSubset$Species <- 
  recode(observed.prop.lossRunganSubset$Call.category, Hfunstart = "NGreyGibbon",
         Hfuntrill = "NGreyGibbon",Pmor='OrangSabah',PwurP='OrangCKalimantan',
         PwurS='OrangCKalimantan',Halbstart='WhiteBeardGibbon',Halbend='WhiteBeardGibbon',Halbpeak='WhiteBeardGibbon')


# Add time category
RunganTimeCats <- read.csv('Rungan_playbackNo_timeCat.csv')

TimeCatsRunganList <- list()
for(a in 1:nrow(observed.prop.lossRunganSubset)){
  Temprow <-  observed.prop.lossRunganSubset[a,]
  TempTime <- subset(RunganTimeCats,recorder== Temprow$playback.num)
  TimeCatsRunganList[[a]] <- TempTime$time.cat
}

observed.prop.lossRunganSubset$TimeCat <- unlist(TimeCatsRunganList)
observed.prop.lossRunganSubset$site <- rep('Rungan',nrow(observed.prop.lossRunganSubset)) 

observed.prop.lossRunganSubset$distance <- log10(observed.prop.lossRunganSubset$distance)
# Do we still need this?? characterize propagation loss as a function of frequency & distance at each site ------------------
## Build 2 models: 1 for Maliau, 1 for Rungan
## Loss of dB from reference recorder = signal - noise of selection
## Frequency = (log) center frequency of selection
## Distance = (log) distance from reference recorder
observed.prop.lossRunganSubset$time <- as.factor(observed.prop.lossRunganSubset$time)
observed.prop.lossRunganSubset$Loc_Name <- as.factor(observed.prop.lossRunganSubset$Loc_Name)
observed.prop.lossRunganSubset$ARU_ID <- as.factor(observed.prop.lossRunganSubset$ARU_ID)
observed.prop.lossRunganSubset$Call.category <- (as.factor(observed.prop.lossRunganSubset$Call.category))
levels(observed.prop.lossRunganSubset$Call.category)

# Include recorder ID as random effect
Rungan.lmm.prop.loss.null <- lmer(magic.x ~  (1|Loc_Name), data=observed.prop.lossRunganSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Rungan.lmm.prop.loss.full <- lmer(magic.x ~ distance + habitat + Species + TimeCat + (1|Loc_Name), data=observed.prop.lossRunganSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Rungan.lmm.prop.loss.full.aru <- lmer(magic.x ~ distance + habitat*Species + TimeCat + (1|Loc_Name)+ (1|ARU_ID), data=observed.prop.lossRunganSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)

Rungan.lmm.prop.loss.notime <- lmer(magic.x ~ distance + habitat + Species  + (1|Loc_Name), data=observed.prop.lossRunganSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Rungan.lmm.prop.loss.nohabitat <- lmer(magic.x ~ distance  + Species + TimeCat +(1|Loc_Name), data=observed.prop.lossRunganSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Rungan.lmm.prop.loss.nodistance <- lmer(magic.x ~  habitat + Species + TimeCat + (1|Loc_Name), data=observed.prop.lossRunganSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)

bbmle::AICctab(Rungan.lmm.prop.loss.null,Rungan.lmm.prop.loss.full,Rungan.lmm.prop.loss.full.aru,
               Rungan.lmm.prop.loss.notime,Rungan.lmm.prop.loss.nohabitat,Rungan.lmm.prop.loss.nodistance, weights=T)


summary(Rungan.lmm.prop.loss.full)
hist(resid(Rungan.lmm.prop.loss.full.aru))
sjPlot::plot_model(Rungan.lmm.prop.loss.full,intercept=F,sort.est = TRUE)+ggtitle('Rungan Propogation Loss')+theme_bw()+geom_hline(yintercept = 0)


ggpubr::ggboxplot(data=observed.prop.lossRunganSubset,
                  x='habitat',y='magic.x',fill='TimeCat')+ylab('Propagation Loss')+
  xlab('Habitat')

ggpubr::ggboxplot(data=observed.prop.lossRunganSubset,
                  x='habitat',y='magic.x',fill='habitat')+ylab('Propagation Loss')+
  xlab('Habitat')

ggpubr::ggboxplot(data=observed.prop.lossRunganSubset,
                  x='distance',y='magic.x',fill='Call.category')+ylab('Propagation Loss')+
  xlab('Distance (log)')

ggpubr::ggboxplot(data=observed.prop.lossRunganSubset,
                  x='distance',y='actual.receive.level',fill='Call.category')+ylab('Propagation Loss')+
  xlab('Distance (log)')
