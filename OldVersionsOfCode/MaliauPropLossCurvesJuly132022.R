library(geosphere)
library(XML)
library(ggplot2)
library(stringr)

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

# Check output
dput((dist.mat+17)[,1])

dist.to.playback.maliau <- 26.4

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
    

    distance.from.source <- dist.mat[c(small.sample.playback.test$recorder),c(small.sample.playback.test$recorder[1])]

     
    
    # Loop to calculate propagation loss; note the index starts at 2 since we use the closest one as the reference
    for(c in 2:length(recorder.index.test)){tryCatch({ 
      print(c)
      # Isolate the recorder that we will use to estimate receive levels
      temp.recorder.received <- subset(small.sample.playback.test,recorder==recorder.index.test[c])
      
      # Isolate the recorder we consider as the 'source' for our relative calculations
      temp.recorder.source <- subset(small.sample.playback.test,recorder==recorder.index.test[1])
      
      # Based on our distance matrix above calculate the distance between the two recorders
      distance.from.source <- dist.mat[c(temp.recorder.received$recorder),c(small.sample.playback.test$recorder[1])]
      
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
      }
    
  }
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}


observed.prop.lossMaliau$ExcessAttenuation <- round(as.numeric(observed.prop.lossMaliau$ExcessAttenuation),2)
observed.prop.lossMaliau$time <- as.factor(observed.prop.lossMaliau$time)
observed.prop.lossMaliau$Call.category <- str_split_fixed(observed.prop.lossMaliau$Sound.type,pattern = '_',n=3)[,2]
observed.prop.lossMaliau$distance <- round(observed.prop.lossMaliau$distance,0)

# Prop loss curves -----------------------------------------------------------------

# Gibbon prop loss --------------------------------------------------------

observed.prop.lossMaliauGibbons <- subset(observed.prop.lossMaliau,Call.category=="Hfunstart" |Call.category=="Hfuntrill" |
                                            Call.category=="Halbstart" |Call.category=="Halbpeak" )

observed.prop.lossMaliauGibbons$Call.category <- 
  recode(observed.prop.lossMaliauGibbons$Call.category, Hfunstart = "NGreyGibbon",
         Hfuntrill = "NGreyGibbon",Halbstart='WhiteBeardGibbon',Halbend='WhiteBeardGibbon',Halbpeak='WhiteBeardGibbon')

observed.prop.lossMaliauGibbons <- subset(observed.prop.lossMaliauGibbons,distance > 100)

uniquegibbons <- unique(observed.prop.lossMaliauGibbons$Call.category)
gibbondB <- 113


for(d in 1:length(uniquegibbons)){
  
  observed.prop.lossMaliaugibbonstemp <- subset(observed.prop.lossMaliauGibbons,Call.category==uniquegibbons[d])
  
  gibbonpropdBMaliau <- data.frame()
  gibbonNoisedf <- list()
  for(e in 1:nrow(observed.prop.lossMaliaugibbonstemp)){
    
    temp.magicx.row <- observed.prop.lossMaliaugibbonstemp[e,]
    
    gibbonNoisedf[[e]] <- temp.magicx.row$noise.level
    
    playback.line.1 <- temp.magicx.row$magic.x
    
    # Set the equations for observed, spherical and cylindrical spreading
    eq1 <- function(x){ playback.line.1*log10(x)}
    
    # Create a series of points based on the above equations
    Estimated1 <- cbind.data.frame(seq(1:500),eq1(1:500),rep('Estimated',500))
    colnames(Estimated1) <- c("X","Value","Label")
    
    Estimated1$X <- Estimated1$X -1
    
    scaledB <- gibbondB - Estimated1$Value[2] 
    Estimated1$Value <- scaledB+Estimated1$Value
    newtemprow <- rbind.data.frame(Estimated1,temp.magicx.row$Call.category)
    
    gibbonpropdBMaliau <- rbind.data.frame(gibbonpropdBMaliau,newtemprow)
  }
  
  gibbonpropdBMaliau$Value <- as.numeric(gibbonpropdBMaliau$Value)
  gibbonpropdBMaliau$X <- as.numeric(gibbonpropdBMaliau$X)
  
  gibbonpropdBMaliauCI <- gibbonpropdBMaliau %>%
    group_by(X) %>%
    summarise(meangibbon = mean(Value, na.rm = TRUE),
              sdgibbon = sd(Value, na.rm = TRUE),
              ngibbon = n()) %>%
    mutate(segibbon = sdgibbon / sqrt(ngibbon),
           lower.cigibbon = meangibbon - qt(1 - (0.05 / 2), ngibbon - 1) * segibbon,
           upper.cigibbon = meangibbon + qt(1 - (0.05 / 2), ngibbon - 1) * segibbon)
  
  
  noise.val <- median(unlist(gibbonNoisedf))
  
  gibbonplot <- ggplot(gibbonpropdBMaliauCI, aes(x = X, y = meangibbon, group = 1)) + 
    geom_line(col='red') + 
    geom_ribbon(aes(ymin = lower.cigibbon, ymax = upper.cigibbon), alpha = 0.25)+
    geom_hline(yintercept=noise.val,linetype="dashed", color = "black")+
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
    xlab("Distance from source (m)") + ylab("Amplitude (dB)")+theme_bw()+ggtitle(paste('Maliau',uniquegibbons[d], gibbondB, 'dB @ 1 m Source level'))+
    ylim(25,125)
  print(gibbonplot)
}


# Plot for ABS conference -------------------------------------------------

# Set the equations for observed, spherical and cylindrical spreading
eq1 <- function(x){ -25*log10(x)}
eq2 <- function(x){ -20*log10(x)}
eq3 <- function(x){ -10*log10(x)}

# Create a series of points based on the above equations
Estimated1 <- cbind.data.frame(seq(1:500),eq1(1:500),rep('Excess',500))
colnames(Estimated1) <- c("X","Value","Label")
Spherical <- cbind.data.frame(seq(1:500),eq2(1:500),rep('Spherical',500))
colnames(Spherical) <- c("X","Value","Label")
Cylindrical <-  cbind.data.frame(seq(1:500),eq3(1:500),rep('Cylindrical',500))
colnames(Cylindrical) <- c("X","Value","Label")

# Combine all three into a single dataframe
attenuation.df <- rbind.data.frame(Estimated1,Spherical)

# Convert to reasonable starting value
attenuation.df$Value <-  attenuation.df$Value+110

# Plot the results
ggplot(data = attenuation.df,aes(x=X, y=Value,group=Label, colour=Label,linetype=Label))+
  geom_line() +theme_bw() + scale_color_manual(values = c("red","black","darkgray"))+
  theme(legend.title = element_blank())+ 
  scale_linetype_manual(values=c( "solid","twodash", "dotted"))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
  xlab("Distance from source (m)") + ylab("Amplitude (dB)")+
  geom_hline(yintercept=45,color='blue')+
  # ylim(-120,0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


