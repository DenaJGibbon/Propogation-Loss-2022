library(lme4)
library(lmerTest)
library(ggplot2)
library(MASS)
library(multcomp)
library(FSA)
library(DHARMa)
library(ggpubr)

# Prep Rungan Data
RunganModelingData <- read.csv('observed.prop.lossRunganJuly18.csv')


# Focus only on primates
observed.prop.lossRunganSubset <- droplevels(subset(RunganModelingData,
                                         Call.category=="Hfunstart" |
                                           Call.category=="Hfuntrill" |
                                           Call.category=="Halbstart" |
                                           Call.category=="Halbpeak" |
                                           Call.category=="Halbend" |
                                           Call.category=="Pmor" |
                                           Call.category=="PwurP" |
                                           Call.category=="PwurS" ))


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


# Prep Maliau Data
MaliauModelingData <- read.csv('observed.prop.lossMaliau.csv')

# Focus only on primates
observed.prop.lossMaliauSubset <- droplevels(subset(MaliauModelingData,
                                                    Call.category=="Hfunstart" |
                                                      Call.category=="Hfuntrill" |
                                                      Call.category=="Halbstart" |
                                                      Call.category=="Halbpeak" |
                                                      Call.category=="Halbend" |
                                                      Call.category=="Pmor"))


# Add time category from 'TimeCategories-Maliau-Rungan.csv'
observed.prop.lossMaliauSubset$time <- as.factor(observed.prop.lossMaliauSubset$time)

MaliauTimes <- c("600", 
                  "640", "800", "840", "1000", "1040", "1440", "1600", "1640")

MaliauCategories <- c("Dawn", 
                 "Dawn", "Morning", "Morning", "Morning", "Morning", "Afternoon", "Afternoon", "Afternoon")

MaliauTimeCats <- cbind.data.frame(MaliauTimes,MaliauCategories)
MaliauTimeCats$MaliauTimes <- as.factor(MaliauTimeCats$MaliauTimes )
MaliauTimeCats$MaliauCategories <- as.factor(MaliauTimeCats$MaliauCategories)


TimeCatsMaliauList <- list()
for(a in 1:nrow(observed.prop.lossMaliauSubset)){
  Temprow <-  observed.prop.lossMaliauSubset[a,]

  TempTime <- subset(MaliauTimeCats,MaliauTimes==
                       as.factor(Temprow$time) )
  
  TimeCatsMaliauList[[a]] <- TempTime$MaliauCategories
}

observed.prop.lossMaliauSubset$TimeCat <- unlist(TimeCatsMaliauList)
observed.prop.lossMaliauSubset$playback.num <- paste(observed.prop.lossMaliauSubset$date,observed.prop.lossMaliauSubset$time,sep='_')
observed.prop.lossMaliauSubset$habitat <- rep('D',nrow(observed.prop.lossMaliauSubset))
observed.prop.lossMaliauSubset$site <- rep('Maliau',nrow(observed.prop.lossMaliauSubset)) 
observed.prop.lossMaliauSubset$distance <- as.factor(observed.prop.lossMaliauSubset$distance)
levels(observed.prop.lossMaliauSubset$distance) <- c("50", "100","150","200",'250','300','350')
observed.prop.lossMaliauSubset$distance <- as.numeric(as.character(observed.prop.lossMaliauSubset$distance))

observed.prop.lossMaliauSubset$distance <- log10(observed.prop.lossMaliauSubset$distance)


### Z scaling frequency & distance data??
#data$center.frequency <- scale(data$center.frequency, center = TRUE, scale = TRUE)
#PropogationLossModelingDF$distance <- scale(PropogationLossModelingDF$distance, center = TRUE, scale = TRUE)



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


Rungan.lmm.prop.loss.null <- lmer(magic.x ~  (1|Loc_Name), data=observed.prop.lossRunganSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Rungan.lmm.prop.loss.full <- lmer(magic.x ~ distance + habitat + Call.category + TimeCat + (1|Loc_Name), data=observed.prop.lossRunganSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Rungan.lmm.prop.loss.notime <- lmer(magic.x ~ distance + habitat + Call.category  + (1|Loc_Name), data=observed.prop.lossRunganSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Rungan.lmm.prop.loss.nohabitat <- lmer(magic.x ~ distance  + Call.category + TimeCat +(1|Loc_Name), data=observed.prop.lossRunganSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Rungan.lmm.prop.loss.nodistance <- lmer(magic.x ~  habitat + Call.category + TimeCat + (1|Loc_Name), data=observed.prop.lossRunganSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)

bbmle::AICctab(Rungan.lmm.prop.loss.null,Rungan.lmm.prop.loss.full,
      Rungan.lmm.prop.loss.notime,Rungan.lmm.prop.loss.nohabitat,Rungan.lmm.prop.loss.nodistance, weights=T)


summary(Rungan.lmm.prop.loss.full)
hist(resid(Rungan.lmm.prop.loss.full))
sjPlot::plot_model(Rungan.lmm.prop.loss.full,intercept=F,sort.est = TRUE)+ggtitle('Rungan Propogation Loss')

ggpubr::ggboxplot(data=observed.prop.lossRunganSubset,
                  x='habitat',y='magic.x',fill='habitat')+ylab('Propagation Loss')+
                    xlab('Habitat')

ggpubr::ggboxplot(data=observed.prop.lossRunganSubset,
                  x='Call.category',y='magic.x',fill='Call.category')+ylab('Propagation Loss')+
  xlab('Call category')



observed.prop.lossMaliauSubset <- subset(observed.prop.lossMaliauSubset,magic.x >= -50 & magic.x <= -10)

#observed.prop.lossMaliauSubset <- subset(observed.prop.lossMaliauSubset, time==600 |time==800 |time==800 |time==1000| time==1600 )

Maliau.lmm.prop.loss.null <- lmer(magic.x ~  (1|date), data=observed.prop.lossMaliauSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Maliau.lmm.prop.loss.full <- lmer(magic.x ~ distance  + Call.category + TimeCat + (1|date), data=observed.prop.lossMaliauSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Maliau.lmm.prop.loss.notime <- lmer(magic.x ~ distance  + Call.category  + (1|date), data=observed.prop.lossMaliauSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Maliau.lmm.prop.loss.nodistance <- lmer(magic.x ~  Call.category + TimeCat + (1|date), data=observed.prop.lossMaliauSubset) # + (Call.Type|recorder.ID + Call.Type|recorder.location)

bbmle::AICctab(Maliau.lmm.prop.loss.null,Maliau.lmm.prop.loss.full,
               Maliau.lmm.prop.loss.notime,Maliau.lmm.prop.loss.nodistance, weights=T)

summary(Maliau.lmm.prop.loss.full)
sjPlot::plot_model(Maliau.lmm.prop.loss.full,intercept=F,sort.est = TRUE)+ggtitle('Maliau propogation loss')+ theme_bw()+geom_hline(yintercept = 0)

ggpubr::ggboxplot(data=observed.prop.lossMaliauSubset,
                  x='TimeCat',y='magic.x',fill='TimeCat')

ggpubr::ggboxplot(data=observed.prop.lossMaliauSubset,
                  x='Call.category',y='magic.x',
                  fill='Call.category') +ylab('Propogation loss')


# 1. compare prop loss across habitats, sound types, times of day ------------------
# Combined sites ----------------------------------------------------------

# Prep Rungan Data
RunganModelingData <- read.csv('observed.prop.lossRunganJuly152022.csv')

# Focus only on primates
observed.prop.lossRunganSubset <- droplevels(subset(RunganModelingData,
                                                    Call.category=="Hfunstart" |
                                                      Call.category=="Hfuntrill" |
                                                      Call.category=="Halbstart" |
                                                      Call.category=="Halbpeak" |
                                                      Call.category=="Halbend" |
                                                      Call.category=="Pmor" ))


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

table(observed.prop.lossRunganSubset$distance)
table(observed.prop.lossRunganSubset$time)

# Prep Maliau Data
MaliauModelingData <- read.csv('observed.prop.lossMaliau.csv')

# Focus only on primates
observed.prop.lossMaliauSubset <- droplevels(subset(MaliauModelingData,
                                                    Call.category=="Hfunstart" |
                                                      Call.category=="Hfuntrill" |
                                                      Call.category=="Halbstart" |
                                                      Call.category=="Halbpeak" |
                                                      Call.category=="Halbend" |
                                                      Call.category=="Pmor" ))


# Add time category from 'TimeCategories-Maliau-Rungan.csv'
observed.prop.lossMaliauSubset$time <- as.factor(observed.prop.lossMaliauSubset$time)

MaliauTimes <- c("600", 
                 "640", "800", "840", "1000", "1040", "1440", "1600", "1640")

MaliauCategories <- c("Dawn", 
                      "Dawn", "Morning", "Morning", "Morning", "Morning", "Afternoon", "Afternoon", "Afternoon")

MaliauTimeCats <- cbind.data.frame(MaliauTimes,MaliauCategories)
MaliauTimeCats$MaliauTimes <- as.factor(MaliauTimeCats$MaliauTimes )
MaliauTimeCats$MaliauCategories <- as.factor(MaliauTimeCats$MaliauCategories)


TimeCatsMaliauList <- list()
for(a in 1:nrow(observed.prop.lossMaliauSubset)){
  Temprow <-  observed.prop.lossMaliauSubset[a,]
  
  TempTime <- subset(MaliauTimeCats,MaliauTimes==
                       as.factor(Temprow$time) )
  
  TimeCatsMaliauList[[a]] <- TempTime$MaliauCategories
}

observed.prop.lossMaliauSubset$TimeCat <- unlist(TimeCatsMaliauList)
observed.prop.lossMaliauSubset$playback.num <- paste(observed.prop.lossMaliauSubset$date,observed.prop.lossMaliauSubset$time,sep='_')
observed.prop.lossMaliauSubset$habitat <- rep('D',nrow(observed.prop.lossMaliauSubset))
observed.prop.lossMaliauSubset$site <- rep('Maliau',nrow(observed.prop.lossMaliauSubset)) 
observed.prop.lossMaliauSubset$distance <- as.factor(observed.prop.lossMaliauSubset$distance)
levels(observed.prop.lossMaliauSubset$distance) <- c("50", "100","150","200",'250','300','350')
observed.prop.lossMaliauSubset$distance <- as.numeric(as.character(observed.prop.lossMaliauSubset$distance))
observed.prop.lossMaliauSubset <- subset(observed.prop.lossMaliauSubset,magic.x >= -45 & magic.x <= -15)

# Combine Rungan and Maliau data for modelling


PropogationLossModelingDF <- rbind.data.frame(observed.prop.lossRunganSubset[,-c(14,15)],observed.prop.lossMaliauSubset)

## LMM or GAM? AIC model selection (include all 2-way interactions?)
## Response: Prop.Loss = dB loss per doubling distance
## Site = Maliau or Mungku Baru
## Habitat = Dipterocarp, Kerangas, Mixed Swamp, Low Pole
## Call type = Hfun start, Hfun peak, Halb start, Halb peak, Pmor roar, Pwur roar, Pwur sigh
## Time = predawn, dawn, morning, midday, afternoon, predusk, dusk, night
## Random effects = recorder ID, recorder location

colnames(PropogationLossModelingDF)
unique(PropogationLossModelingDF$distance)

# Subset to compare only comparable distances
unique(PropogationLossModelingDF$distance)
str(PropogationLossModelingDF$distance)
PropogationLossModelingDFCompDist <- subset(PropogationLossModelingDF,distance==100 | distance==250)

PropogationLossModelingDFCompDist$distance <- log10(PropogationLossModelingDFCompDist$distance)
PropogationLossModelingDFCompDist$site <- as.factor(PropogationLossModelingDFCompDist$site)
PropogationLossModelingDFCompDist$TimeCat <- as.factor(PropogationLossModelingDFCompDist$TimeCat)
levels(PropogationLossModelingDFCompDist$TimeCat)

## Example model = LMM
str(PropogationLossModelingDFCompDist)

PropogationLossModelingDFCompDist$habitat <- 
  as.factor(PropogationLossModelingDFCompDist$habitat)

unique(PropogationLossModelingDFCompDist$Call.category)


Combined.lmm.prop.loss.null <- lmer(magic.x ~  (1|date), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.full <- lmer(magic.x ~ distance + habitat + Call.category + TimeCat + (1|date), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.notime <- lmer(magic.x ~ distance + habitat + Call.category  + (1|date), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.nohabitat <- lmer(magic.x ~ distance  + Call.category + TimeCat +(1|date), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.nodistance <- lmer(magic.x ~  habitat + Call.category + TimeCat + (1|date), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)

bbmle::AICctab(Combined.lmm.prop.loss.null,Combined.lmm.prop.loss.full,
               Combined.lmm.prop.loss.notime,Combined.lmm.prop.loss.nohabitat,Combined.lmm.prop.loss.nodistance, weights=T)


summary(lmm.prop.loss.full)
sjPlot::plot_model(lmm.prop.loss.full,intercept=F,sort.est = TRUE)+ggtitle('All sites model coefficients')+theme_bw()+geom_hline(yintercept=0)

ggboxplot(data=PropogationLossModelingDFCompDist,x='Call.category',y='magic.x',fill ='habitat' )

plot(lmm.prop.loss.full)

require(DHARMa)
fittedModel = lmm.prop.loss.full
simulationOutput=simulateResiduals(fittedModel=fittedModel)
plot(simulationOutput)

### Use this for posthoc tests for variables of interest, e.g., call type
summary(glht(lmm.prop.loss.full, linfct = mcp(Call.Type = "Tukey")), test = adjusted("holm"))

### Plotting data
plot.prop.loss <- ggboxplot(data=typical.pulses,x = 'Call.Type', y = 'Prop.Loss', fill = 'Habitat', font.label = list(size = 16, face = "plain"), xlab = "Call type", ylab = "dB loss per doubling distance")


# 2  compare background noise across habitats, sound types, time ------------------
## LMM or GAM? AIC model selection (include all 2-way interactions?)
## Response: Noise = SPL of background noise in third-octave bands
## Frequency bin = Third-octave bands? or avg. center frequency of each call type?
## Site = Maliau or Mungku Baru
## Habitat = Dipterocarp, Kerangas, Mixed Swamp, Low Pole
## Time = predawn, dawn, morning, midday, afternoon, predusk, dusk, night
## Random effects = recorder ID, recorder location



# 3. Boxplots of magic x -----------------------------------------
head(PropogationLossModelingDF)

ggpubr::ggboxplot(data=PropogationLossModelingDF,
                  x='Call.category',y='magic.x',facet.by ='habitat' )+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


table(PropogationLossModelingDF$site,PropogationLossModelingDF$TimeCat)


# Boxplots of detection range ---------------------------------------------
gibbondB <- 100

detectionrange.df <- data.frame()

for(a in 1:nrow(PropogationLossModelingDFCompDist)){
  TempRow <- PropogationLossModelingDFCompDist[a,]
  TempMagic.x <- TempRow$magic.x
  Temp.background <- TempRow$noise.level
  print(TempRow[a,]$noise.level)
  if(is.na(Temp.background)==F){
  # Set the equations for observed, spherical and cylindrical spreading
  eq1 <- function(x){ TempMagic.x*log10(x)}
  
  Estimated1 <- cbind.data.frame(seq(1:1000),eq1(1:1000),rep('Estimated',1000))
  colnames(Estimated1) <- c("X","Value","Label")
  
  Estimated1$X <- Estimated1$X -1
  
  scaledB <- gibbondB - Estimated1$Value[2] 
  Estimated1$Value <- scaledB+Estimated1$Value
  
  TempRow$detect.distance <- 
    which(abs(Estimated1$Value - Temp.background) == min(abs(Estimated1$Value - Temp.background)))
  
  
  detectionrange.df <- rbind.data.frame(detectionrange.df,TempRow)
  }
}


detectionrange.df.rungan <- subset(detectionrange.df,site=='Rungan')
ggboxplot(data=detectionrange.df.rungan,x='Call.category',y='detect.distance',
          fill = 'Call.category')+ylab('Detection Range (m)')+ggtitle('Rungan')

detectionrange.df.maliau <- subset(detectionrange.df,site=='Maliau')
ggboxplot(data=detectionrange.df.maliau,x='Call.category',y='detect.distance',
          fill = 'Call.category')+ylab('Detection Range (m)')+ggtitle('Maliau')

