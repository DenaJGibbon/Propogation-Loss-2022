library(plyr)
# Prop loss curves -----------------------------------------------------------------
observed.prop.lossMaliauGibbons <- observed.prop.lossMaliau #read.csv('observed.prop.lossMaliauAugust4.csv')

# Gibbon prop loss --------------------------------------------------------

observed.prop.lossMaliauGibbons$Species <- 
  recode(observed.prop.lossMaliauGibbons$Call.category, Hfunstart = "NGreyGibbon",
         Hfuntrill = "NGreyGibbon",Halbstart='WhiteBeardGibbon',Halbend='WhiteBeardGibbon',Halbpeak='WhiteBeardGibbon')


uniquegibbons <- unique(observed.prop.lossMaliauGibbons$Species)
gibbondB <- 100


for(d in 1:length(uniquegibbons)){
  
  observed.prop.lossMaliaugibbonstemp <- subset(observed.prop.lossMaliauGibbons,Species==uniquegibbons[d])
  
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
    dplyr::summarise(meangibbon = mean(Value, na.rm = TRUE),
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


