library(lme4)
library(lmerTest)
library(ggplot2)
library(MASS)
library(multcomp)
library(FSA)
library(DHARMa)
library(ggpubr)
library(dplyr)

# Combine Rungan and Maliau data for modelling


PropogationLossModelingDF <- rbind.data.frame(observed.prop.lossRunganSubset[,c("zero.receive.level", "actual.receive.level", "source.level", 
                                                                                "distance", "Sound.type",  "date", "magic.x", "ExcessAttenuation", 
                                                                                "dBdoubledist", "Call.category", "Species", "TimeCat", "playback.num", 
                                                                                "habitat", "site")],observed.prop.lossMaliauSubset[,c("zero.receive.level", "actual.receive.level", "source.level", 
                                                                                                                                      "distance", "Sound.type",  "date", "magic.x", "ExcessAttenuation", 
                                                                                                                                      "dBdoubledist", "Call.category", "Species", "TimeCat", "playback.num", 
                                                                                                                                      "habitat", "site")])

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
PropogationLossModelingDFCompDist <- subset(PropogationLossModelingDF,distance==2 | distance==2.39794)
PropogationLossModelingDFCompDist <- subset(PropogationLossModelingDF,Species=="NGreyGibbon" | Species=="WhiteBeardGibbon")

PropogationLossModelingDFCompDist$distance <- log10(PropogationLossModelingDFCompDist$distance)
PropogationLossModelingDFCompDist$site <- as.factor(PropogationLossModelingDFCompDist$site)
PropogationLossModelingDFCompDist$TimeCat <- as.factor(PropogationLossModelingDFCompDist$TimeCat)
levels(PropogationLossModelingDFCompDist$TimeCat)

## Example model = LMM
str(PropogationLossModelingDFCompDist)

PropogationLossModelingDFCompDist$habitat <- 
  as.factor(PropogationLossModelingDFCompDist$habitat)

unique(PropogationLossModelingDFCompDist$Call.category)

# Remove outliers for now
# PropogationLossModelingDFCompDist <- 
#   subset(PropogationLossModelingDFCompDist,magic.x >= -40 & magic.x <= -20)

Combined.lmm.prop.loss.null <- lmer(magic.x ~  (1|date), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.full <- lmer(magic.x ~ distance + habitat + Species + TimeCat + (1|date), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.notime <- lmer(magic.x ~ distance + habitat + Species  + (1|date), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.nohabitat <- lmer(magic.x ~ distance  + Species + TimeCat +(1|date), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.nodistance <- lmer(magic.x ~  habitat + Species + TimeCat + (1|date), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.interaction <- lmer(magic.x ~  habitat*Species + TimeCat + (1|date), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.habitatonly <- lmer(magic.x ~  habitat + (1|date), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)
Combined.lmm.prop.loss.habitatdist <- lmer(magic.x ~  habitat + distance + (1|date), data=PropogationLossModelingDFCompDist) # + (Call.Type|recorder.ID + Call.Type|recorder.location)

bbmle::AICctab(Combined.lmm.prop.loss.null,Combined.lmm.prop.loss.full,Combined.lmm.prop.loss.interaction,Combined.lmm.prop.loss.habitatdist,
               Combined.lmm.prop.loss.notime,Combined.lmm.prop.loss.nohabitat,Combined.lmm.prop.loss.nodistance, Combined.lmm.prop.loss.habitatonly,weights=T)


hist(resid(Combined.lmm.prop.loss.full ))

sjPlot::plot_model(Combined.lmm.prop.loss.interaction ,intercept=F,sort.est = TRUE)

ggboxplot(data=PropogationLossModelingDFCompDist,x='Call.category',y='magic.x',fill ='habitat',
          outlier.shape = NA)+ylab('Propagation loss')

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

sourcelevels <- c(95, 100, 105)

detectionrange.df <- data.frame()

for(a in 1:nrow(PropogationLossModelingDF)){
  for(b in 1:length(sourcelevels)){
  TempRow <- PropogationLossModelingDF[a,]
  Subsetformedianbackground <- subset(PropogationLossModelingDF, site==TempRow$site & Call.category==TempRow$Call.category)
  
  TempMagic.x <- TempRow$magic.x
  Temp.background <-  median(Subsetformedianbackground$noise.level)
  if(is.na(Temp.background)==F){
  # Set the equations for observed, spherical and cylindrical spreading
  eq1 <- function(x){ TempMagic.x*log10(x)}
  
  Estimated1 <- cbind.data.frame(seq(1:2500),eq1(1:2500),rep('Estimated',2500))
  colnames(Estimated1) <- c("X","Value","Label")
  
  Estimated1$X <- Estimated1$X -1
  
  gibbondB <- sourcelevels[b]
  scaledB <- gibbondB - Estimated1$Value[2] 
  Estimated1$Value <- scaledB+Estimated1$Value
  
  TempRow$detect.distance <- 
    which(abs(Estimated1$Value - Temp.background) == min(abs(Estimated1$Value - Temp.background)))
  
  TempRow$gibbondB <- gibbondB
  
  detectionrange.df <- rbind.data.frame(detectionrange.df,TempRow)
  }
  }
}


detectionrange.df.rungan <- subset(detectionrange.df,site=='Rungan')
ggboxplot(data=detectionrange.df.rungan,x='Call.category',y='detect.distance',
          fill = 'Call.category',facet.by = 'gibbondB',outlier.shape =NA)+ylab('Detection Range (m)')+ggtitle('Rungan')+ theme(axis.text.x = element_text(angle = 45, hjust=1))+xlab('')

detectionrange.df.maliau <- subset(detectionrange.df,site=='Maliau')
ggboxplot(data=detectionrange.df.maliau,x='Call.category',y='detect.distance',
          fill = 'Call.category',facet.by = 'gibbondB',outlier.shape =NA)+ylab('Detection Range (m)')+ggtitle('Maliau')+ theme(axis.text.x = element_text(angle = 45, hjust=1))+xlab('')

