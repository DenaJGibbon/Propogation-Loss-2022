library(bbmle)

# Combine Rungan and Maliau data for modeling
NoiseModelingDF <-
  rbind.data.frame(
    observed.prop.lossRunganOutRM[, c(
      "Species",
      "TimeCat",
      "playback.num",
      "habitat",
      "site",
      "noise.level"
    )],
    observed.prop.lossMaliauOutRM[, c("Species", "TimeCat", "playback.num", 
      "habitat", "site","noise.level")])


NoiseModelingDF <- PropogationLossModelingDF

# Create boxplots by Species ----------------------------------------------
# # Remove outlies
# outliers <- boxplot.stats(NoiseModelingDF$noise.level)$out
# IndexRM <- which(NoiseModelingDF$noise.level %in% outliers)
# NoiseModelingDF <- NoiseModelingDF[-IndexRM,]

#NoiseModelingDF$noise.level <- log(NoiseModelingDF$noise.level)

NoiseModelingDF$Species <- as.factor(NoiseModelingDF$Species)

NoiseModelingDF$Species  <- factor(NoiseModelingDF$Species , levels = c("NGreyGibbon","WhiteBeardGibbon",
                                                                        "OrangSabah" , "OrangCKalimantan" ))


ggpubr::ggboxplot(x='Species',y='noise.level',fill='habitat',
                  data=NoiseModelingDF,outlier.shape = NA)+
                  ylab('Ambient noise \n (dB re 20 µPa) ')+
                  theme(legend.title=element_blank())#+ stat_compare_means(paired = TRUE) 


NoiseFullModel <- lmer(noise.level ~ habitat+ TimeCat+  Species+(1|site/date),
                         data=NoiseModelingDF)

NoiseInterModel <- lmer(noise.level ~ habitat*Species+(1|site/date),
                       data=NoiseModelingDF)

NoiseNullModel <- lmer(noise.level ~  (1|site/date),
                       data=NoiseModelingDF)

AICctab(NoiseFullModel,NoiseNullModel,NoiseInterModel,weights=T)

sjPlot::plot_model(NoiseInterModel, type='est')+theme_bw()+geom_hline(yintercept = 0)+ ggtitle('Background noise coefficient plot')

hist(resid(NoiseInterModel))

simulationOutputNoise <- simulateResiduals(fittedModel = NoiseFullModel, plot = F)
plotQQunif(simulationOutputNoise)
