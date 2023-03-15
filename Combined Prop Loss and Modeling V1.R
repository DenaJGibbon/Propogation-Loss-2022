library(lme4)
library(lmerTest)
library(ggplot2)
library(MASS)
library(multcomp)
library(FSA)
library(DHARMa)
library(ggpubr)
library(dplyr)

# Combine Rungan and Maliau data for modeling
PropogationLossModelingDF <-
  rbind.data.frame(
    observed.prop.lossRunganOutRM[, c(
      "zero.receive.level",
      "actual.receive.level",
      "magic.x",
      "source.level",
      "log.distance",
      'distance',
      "Sound.type",
      "date",
      "dBdoubledist",
      "ExcessAttenuation",
      "dBdoubledist",
      "Call.category",
      "Species",
      "TimeCat",
      "playback.num",
      "habitat",
      "site",
      "noise.level"
    )],
    observed.prop.lossMaliauOutRMdf[, c(
      "zero.receive.level",
      "actual.receive.level",
      "magic.x",
      "source.level",
      "log.distance",
      'distance',
      "Sound.type",
      "date",
      "dBdoubledist",
      "ExcessAttenuation",
      "dBdoubledist", "Call.category", "Species", "TimeCat", "playback.num", 
          "habitat", "site","noise.level")])



## LMM or GAM? AIC model selection (include all 2-way interactions?)
## Response: Prop.Loss = dB loss per doubling log.distance
## Site = Maliau or Mungku Baru
## Habitat = Dipterocarp, Kerangas, Mixed Swamp, Low Pole
## Call type = Hfun start, Hfun peak, Halb start, Halb peak, Pmor roar, Pwur roar, Pwur sigh
## Time = predawn, dawn, morning, midday, afternoon, predusk, dusk, night
## Random effects = recorder ID, recorder location

colnames(PropogationLossModelingDF)
unique(PropogationLossModelingDF$log.distance)
unique(PropogationLossModelingDF$distance)
# Subset to compare only comparable log.distances
unique(PropogationLossModelingDF$log.distance)
str(PropogationLossModelingDF$log.distance)
PropogationLossModelingDF$Species <- as.factor(PropogationLossModelingDF$Species)

PropogationLossModelingDF$distance <- as.factor(PropogationLossModelingDF$distance)

PropogationLossModelingDFCompDist <- subset(PropogationLossModelingDF,distance=='89' | distance=='100.498756211209'| distance=='250.199920063936'| distance=='282')
PropogationLossModelingDFCompDistOrang <- subset(PropogationLossModelingDFCompDist,Species!="NGreyGibbon" & Species!="WhiteBeardGibbon")


PropogationLossModelingDFCompDist <- subset(PropogationLossModelingDFCompDist,Species=="NGreyGibbon" | Species=="WhiteBeardGibbon")

PropogationLossModelingDFCompDist$site <- as.factor(PropogationLossModelingDFCompDist$site)
PropogationLossModelingDFCompDist$TimeCat <- as.factor(PropogationLossModelingDFCompDist$TimeCat)
levels(PropogationLossModelingDFCompDist$TimeCat)
hist((PropogationLossModelingDFCompDist$dBdoubledist))

## Example model = LMM
str(PropogationLossModelingDFCompDist)

PropogationLossModelingDFCompDist$habitat <- 
  as.factor(PropogationLossModelingDFCompDist$habitat)

unique(PropogationLossModelingDFCompDist$distance)

# Remove outliers for now
# PropogationLossModelingDFCompDist <- 
#   subset(PropogationLossModelingDFCompDist,dBdoubledist >= -40 & dBdoubledist <= -20)

#PropogationLossModelingDFCompDist$habitat <- factor(PropogationLossModelingDFCompDist$habitat, levels = c("K", "LP", "MS", "D"))

Combined.lmm.prop.loss.null <- lmer(dBdoubledist ~  (1|site/playback.num), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.full <- lmer(dBdoubledist ~ log.distance + habitat + Species + TimeCat + (1|site/playback.num), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.notime <- lmer(dBdoubledist ~ log.distance + habitat + Species  + (1|site/playback.num), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.nohabitat <- lmer(dBdoubledist ~ log.distance  + Species + TimeCat +(1|site/playback.num), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.nolog.distance <- lmer(dBdoubledist ~  habitat + Species + TimeCat + (1|site/playback.num), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
#Combined.lmm.prop.loss.interaction <- lmer(dBdoubledist ~  habitat*Species+ log.distance+ (1|site/playback.num), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.habitatonly <- lmer(dBdoubledist ~  habitat + (1|site/playback.num), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.habitatdist <- lmer(dBdoubledist ~  habitat + log.distance + (1|site/playback.num), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
#Combined.lmm.prop.loss.habitatdist <- lmer(dBdoubledist ~  habitat + log.distance + Call.category+(1|site/playback.num), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)

bbmle::AICctab(Combined.lmm.prop.loss.null,Combined.lmm.prop.loss.full,Combined.lmm.prop.loss.habitatdist,
               Combined.lmm.prop.loss.notime,Combined.lmm.prop.loss.nohabitat,Combined.lmm.prop.loss.nolog.distance, Combined.lmm.prop.loss.habitatonly,weights=T)


hist(resid(Combined.lmm.prop.loss.notime  ))

simulationOutput <- simulateResiduals(fittedModel = Combined.lmm.prop.loss.notime , plot = F)
plotQQunif(simulationOutput)

MaliauMungkBaruComb <- sjPlot::plot_model(Combined.lmm.prop.loss.notime  ,sort.est = F)+ggtitle('Combined Maliau and Mungku Baru')+
  theme_bw()+geom_hline(yintercept=0,color='black')
MaliauMungkBaruComb

ggboxplot(data=PropogationLossModelingDFCompDist,x='Species',y='dBdoubledist',fill='habitat')
  
Plotlist <-sjPlot::plot_model(Combined.lmm.prop.loss.notime, type='pred')

# Marginal effects tells us how a dependent variable (outcome) changes when a specific independent variable (explanatory variable) changes.
CombinedGibbonMarginalEffectsDist <- Plotlist[[1]]+ylab('Propagation loss \n (dB decrease per doubling distance)')+
  ggtitle('')+xlab('Log Distance')+theme_bw()

CombinedGibbonMarginalEffectsHabitat <- Plotlist[[2]]+ylab('Propagation loss \n (dB decrease per doubling distance)')+
  ggtitle('')+xlab('Habitat')+theme_bw()

CombinedGibbonMarginalEffectsSpecies <- Plotlist[[3]]+ylab('Propagation loss \n (dB decrease per doubling distance)')+
  ggtitle('')+theme_bw()

cowplot::plot_grid(CombinedGibbonMarginalEffectsDist,
                   CombinedGibbonMarginalEffectsHabitat,
                   CombinedGibbonMarginalEffectsSpecies,
                   labels=c('A','B','C'),label_x = 0.9)

ggboxplot(data=PropogationLossModelingDFCompDist,x='Species',y='dBdoubledist',fill ='habitat',
          outlier.shape = NA)+ylab('Propagation loss')



# 2  compare background DetectRangeGibbon across habitats, sound types, time ------------------
## LMM or GAM? AIC model selection (include all 2-way interactions?)
## Response: DetectRangeGibbon = SPL of background DetectRangeGibbon in third-octave bands
## Frequency bin = Third-octave bands? or avg. center frequency of each call type?
## Site = Maliau or Mungku Baru
## Habitat = Dipterocarp, Kerangas, Mixed Swamp, Low Pole
## Time = predawn, dawn, morning, midday, afternoon, predusk, dusk, night
## Random effects = recorder ID, recorder location



# 3. Boxplots of magic x -----------------------------------------
head(PropogationLossModelingDF)

ggpubr::ggboxplot(data=PropogationLossModelingDFCompDist,
                  x='Species',y='dBdoubledist',#facet.by ='habitat',
                  fill='habitat' )




# Boxplots of detection range for gibbons ---------------------------------------------

sourcelevels <- c(113.9)

detectionrange.df <- data.frame()

for(a in 1:nrow(PropogationLossModelingDFCompDist)){
  for(b in 1:length(sourcelevels)){
  TempRow <- PropogationLossModelingDFCompDist[a,]
  Subsetformedianbackground <- subset(PropogationLossModelingDF, site==TempRow$site & Call.category==TempRow$Call.category)
  
  TempdBmagicx <- TempRow$magic.x
  Temp.background <-  median(Subsetformedianbackground$noise.level)
  if(is.na(Temp.background)==F){
  # Set the equations for observed, spherical and cylindrical spreading
  eq1 <- function(x){ TempdBmagicx*log10(x)}
  
  Estimated1 <- cbind.data.frame(seq(1:4000),eq1(1:4000),rep('Estimated',4000))
  colnames(Estimated1) <- c("X","Value","Label")
  
  Estimated1$X <- Estimated1$X -1
  
  gibbondB <- sourcelevels[b]
  scaledB <- gibbondB - Estimated1$Value[2] 
  Estimated1$Value <- scaledB+Estimated1$Value
  
  TempRow$detect.log.distance <- 
    which.min(abs(Estimated1$Value - TempRow$noise.level))
  
  TempRow$gibbondB <- gibbondB
  
  detectionrange.df <- rbind.data.frame(detectionrange.df,TempRow)
  }
  }
}

detectionrange.df.gibbon <- droplevels(subset(detectionrange.df,Species=='NGreyGibbon'|Species=="WhiteBeardGibbon"))
detectionrange.df.gibbon <- subset(detectionrange.df.gibbon,detect.log.distance < 1500)

ggboxplot(data=detectionrange.df.gibbon,x='habitat',y='detect.log.distance',
          fill = 'Species',facet.by = 'gibbondB',outlier.shape =NA)+
  ylab('Detection Range (m)')+ggtitle('Gibbon detection range')+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1))+xlab('')#+ylim(0,1000)

detectionrange.df.gibbon$detect.log.distance <-log(detectionrange.df.gibbon$detect.log.distance)

# Gibbon detection range model selection
DetectRangeGibbonFullModel <- lmer(detect.log.distance ~ habitat+ TimeCat+  Species+(1|site/date),
                       data=detectionrange.df.gibbon)

DetectRangeGibbonInterModel <- lmer(detect.log.distance ~ habitat*Species+(1|site/date),
                                    data=detectionrange.df.gibbon)

DetectRangeGibbonNullModel <- lmer(detect.log.distance ~  (1|site/date),
                       data=detectionrange.df.gibbon)

AICctab(DetectRangeGibbonFullModel,DetectRangeGibbonNullModel,DetectRangeGibbonInterModel,weights=T)

sjPlot::plot_model(DetectRangeGibbonFullModel, type='est')+theme_bw()+geom_hline(yintercept = 0)+ ggtitle('Background DetectRangeGibbon coefficient plot')

hist(resid(DetectRangeGibbonFullModel))

simulationOutputDetectRangeGibbon <- simulateResiduals(fittedModel = DetectRangeGibbonFullModel, plot = F)
plotQQunif(simulationOutputDetectRangeGibbon)

# Orangutan

sourcelevels <- c(105)

detectionrange.df.Orangutan <- data.frame()

for(a in 1:nrow(PropogationLossModelingDFCompDistOrang)){
  for(b in 1:length(sourcelevels)){
    TempRow <- PropogationLossModelingDFCompDistOrang[a,]
    Subsetformedianbackground <- subset( PropogationLossModelingDFCompDistOrang, site==TempRow$site & Call.category==TempRow$Call.category)
    
    TempdBmagicx <- TempRow$magic.x
    Temp.background <-  median(Subsetformedianbackground$noise.level)
    if(is.na(Temp.background)==F){
      # Set the equations for observed, spherical and cylindrical spreading
      eq1 <- function(x){ TempdBmagicx*log10(x)}
      
      Estimated1 <- cbind.data.frame(seq(1:4000),eq1(1:4000),rep('Estimated',4000))
      colnames(Estimated1) <- c("X","Value","Label")
      
      Estimated1$X <- Estimated1$X -1
      
      orangutandB <- sourcelevels[b]
      scaledB <- orangutandB - Estimated1$Value[2] 
      Estimated1$Value <- scaledB+Estimated1$Value
      
      TempRow$detect.log.distance <- 
        which.min(abs(Estimated1$Value - TempRow$noise.level))
      
      TempRow$orangutandB <- orangutandB
      
      detectionrange.df.Orangutan <- rbind.data.frame(detectionrange.df.Orangutan,TempRow)
    }
  }
}

hist(detectionrange.df.Orangutan$detect.log.distance)
detectionrange.df.Orangutan <- subset(detectionrange.df.Orangutan,detect.log.distance < 3000)

ggboxplot(data=detectionrange.df.Orangutan,x='habitat',y='detect.log.distance',
          fill = "Call.category",facet.by = 'OrangutandB',outlier.shape =NA)+
  ylab('Detection Range (m)')+ggtitle('Orangutan detection range')+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1))+xlab('')

ggboxplot(data=detectionrange.df.Orangutan,x='habitat',y='dBdoubledist',
          fill = 'Call.category',outlier.shape =NA)+
  ylab('dB decrease per doubling log.distance')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+xlab('')

ggboxplot(data=detectionrange.df.Orangutan,x='habitat',y='detect.log.distance',
          fill = 'Species',outlier.shape =NA)+
  ylab('Detection Range (m)')+ggtitle('Orangutan detection range')+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1))+xlab('')#+ylim(0,1000)

detectionrange.df.Orangutan$detect.log.distance <-log(detectionrange.df.Orangutan$detect.log.distance)

# Orangutan detection range model selection
DetectRangeOrangutanFullModel <- lmer(detect.log.distance ~ habitat+ TimeCat+  Species+(1|date),
                                   data=detectionrange.df.Orangutan)

DetectRangeOrangutanInterModel <- lmer(detect.log.distance ~ habitat*Species+(1|date),
                                    data=detectionrange.df.Orangutan)

DetectRangeOrangutanNullModel <- lmer(detect.log.distance ~  (1|date),
                                   data=detectionrange.df.Orangutan)

AICctab(DetectRangeOrangutanFullModel,DetectRangeOrangutanNullModel,DetectRangeOrangutanInterModel,weights=T)

sjPlot::plot_model(DetectRangeOrangutanInterModel, type='est')+theme_bw()+geom_hline(yintercept = 0)+ ggtitle('Detection Range Orangutan coefficient plot')

hist(resid(DetectRangeOrangutanFullModel))

simulationOutputDetectRangeOrangutan <- simulateResiduals(fittedModel = DetectRangeOrangutanFullModel, plot = F)
plotQQunif(simulationOutputDetectRangeOrangutan)

