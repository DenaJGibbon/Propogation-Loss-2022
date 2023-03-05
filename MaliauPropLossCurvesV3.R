library(plyr)
library(ggplot2)
# Prop loss curves -----------------------------------------------------------------
observed.prop.lossMaliauGibbons <- droplevels(subset(observed.prop.lossMaliauSubset,
                                          Species != "OrangSabah"))#read.csv('observed.prop.lossMaliauAugust4.csv')

# Gibbon prop loss --------------------------------------------------------
# Set the equations for observed, spherical and cylindrical spreading
eq1 <- function(x){ median(observed.prop.lossMaliauGibbons$magic.x) *log10(x)}
eq2 <- function(x){ -20*log10(x)}
eq3 <- function(x){ -10*log10(x)}

uniquegibbons <- unique(observed.prop.lossMaliauGibbons$Species)
gibbondB <- 107.1
gibbondB <- 113.9

gibbonplotlist <- list()

for(d in 1:length(uniquegibbons)){
  
  observed.prop.lossMaliaugibbonstemp <- subset(observed.prop.lossMaliauGibbons,Species==uniquegibbons[d])
  
  gibbonpropdBMaliau <- data.frame()
  gibbonNoisedf <- list()
  # gibbonNoisedfMorning <- list()
  # gibbonNoisedfAfternoon <- list()
  
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
    dplyr::summarise(meangibbon = mean(Value, na.rm = TRUE),
              sdgibbon = sd(Value, na.rm = TRUE),
              ngibbon = n()) %>%
    mutate(segibbon = sdgibbon / sqrt(ngibbon),
           lower.cigibbon = meangibbon - qt(1 - (0.05 / 2), ngibbon - 1) * segibbon,
           upper.cigibbon = meangibbon + qt(1 - (0.05 / 2), ngibbon - 1) * segibbon)
  
  
  noise.val <- median(unlist(gibbonNoisedf))
  senoise = sd(unlist(gibbonNoisedf)) / sqrt(length(gibbonNoisedf))
  gibbonpropdBMaliauCI$lower.cinoise = noise.val - qt(1 - (0.05 / 2), length(gibbonNoisedf) - 1) * senoise
  gibbonpropdBMaliauCI$upper.cinoise = noise.val + qt(1 - (0.05 / 2), length(gibbonNoisedf) - 1) * senoise
  
  gibbonplot <- ggplot(gibbonpropdBMaliauCI, aes(x = X, y = meangibbon, group = 1)) + 
    geom_line(col='red') + 
    geom_ribbon(aes(ymin = lower.cigibbon, ymax = upper.cigibbon), alpha = 0.25)+
    #geom_hline(yintercept=noise.val,linetype="dashed", color = "black")+
    geom_ribbon(aes(ymin = lower.cinoise, ymax = upper.cinoise), alpha = 0.5, color='blue')+
    
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
    xlab("Distance from source (m)") + ylab("Amplitude (dB)")+theme_bw()+ggtitle(paste('D',uniquegibbons[d]))+
    ylim(25,gibbondB)
  
   gibbonplotlist[[d]] <- gibbonplot
}


cowplot::plot_grid(gibbonplotlist[[1]],gibbonplotlist[[2]])


