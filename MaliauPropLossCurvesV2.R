library(plyr)
library(ggplot2)
# Prop loss curves -----------------------------------------------------------------
observed.prop.lossMaliauGibbons <- observed.prop.lossMaliau #read.csv('observed.prop.lossMaliauAugust4.csv')

# Gibbon prop loss --------------------------------------------------------

observed.prop.lossMaliauGibbons$Species <- 
  recode(observed.prop.lossMaliauGibbons$Call.category, Hfunstart = "NGreyGibbon",
         Hfuntrill = "NGreyGibbon",Halbstart='WhiteBeardGibbon',Halbend='WhiteBeardGibbon',Halbpeak='WhiteBeardGibbon')

observed.prop.lossMaliauGibbons$time <- as.numeric(as.character(observed.prop.lossMaliauGibbons$time))

observed.prop.lossMaliauGibbons <-  observed.prop.lossMaliauGibbons %>%
  mutate(TimeCat = case_when(
    time <= 700  ~ 'Dawn',
    time >= 800 & time <= 1100 ~ 'Morning',
    TRUE ~ 'Afternoon'
  ))

uniquegibbons <- unique(observed.prop.lossMaliauGibbons$Species)
gibbondB <- 106


for(d in 1:length(uniquegibbons)){
  
  observed.prop.lossMaliaugibbonstemp <- subset(observed.prop.lossMaliauGibbons,Species==uniquegibbons[d])
  
  gibbonpropdBMaliau <- data.frame()
  gibbonNoisedf <- list()
  # gibbonNoisedfMorning <- list()
  # gibbonNoisedfAfternoon <- list()
  
  for(e in 1:nrow(observed.prop.lossMaliaugibbonstemp)){
    
    temp.magicx.row <- observed.prop.lossMaliaugibbonstemp[e,]
 
    gibbonNoisedf[[e]] <- temp.magicx.row$noise.level
    # if(temp.magicx.row$TimeCat=='Dawn'){
    # gibbonNoisedfDawn[[e]] <- temp.magicx.row$noise.level
    # }
    # 
    # if(temp.magicx.row$TimeCat=='Morning'){
    #   gibbonNoisedfMorning[[e]] <- temp.magicx.row$noise.level
    # }
    # 
    # if(temp.magicx.row$TimeCat=='Afternoon'){
    #   gibbonNoisedfAfternoon[[e]] <- temp.magicx.row$noise.level
    # }
    # 
    playback.line.1 <- temp.magicx.row$magic.x
    
    # Set the equations for observed, spherical and cylindrical spreading
    eq1 <- function(x){ playback.line.1*log10(x)}
    
    # Create a series of points based on the above equations
    Estimated1 <- cbind.data.frame(seq(3:500),eq1(3:500),rep('Estimated',500-2))
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
  # noise.val.morning <- median(unlist(gibbonNoisedfMorning))
  # noise.val.afternoon <- median(unlist(gibbonNoisedfAfternoon))
  # 
  gibbonplot <- ggplot(gibbonpropdBMaliauCI, aes(x = X, y = meangibbon, group = 1)) + 
    geom_line(col='red') + 
    geom_ribbon(aes(ymin = lower.cigibbon, ymax = upper.cigibbon), alpha = 0.25)+
    geom_hline(yintercept=noise.val.dawn,linetype="dashed", color = "black")+
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
    xlab("Distance from source (m)") + ylab("Amplitude (dB)")+theme_bw()+ggtitle(paste('Maliau',uniquegibbons[d], gibbondB, 'dB @ 3 m Source level'))+
    ylim(25,gibbondB)
  
  Spherical <- cbind.data.frame(seq(3:500),eq2(3:500),rep('Spherical',500-2))
  colnames(Spherical) <- c("X","Value","Label")
  
  Spherical$Value <- Spherical$Value + 120.5759
  
  gibbonplot +geom_line(data=Spherical,aes(x = X, y = Value),col='black',linetype='dashed')
  print(gibbonplot)
}


# Plot for ABS conference -------------------------------------------------

# Set the equations for observed, spherical and cylindrical spreading
eq1 <- function(x){ median(observed.prop.lossMaliauGibbons$magic.x) *log10(x)}
eq2 <- function(x){ -20*log10(x)}
eq3 <- function(x){ -10*log10(x)}

# Create a series of points based on the above equations
Estimated1 <- cbind.data.frame(seq(3:500),eq1(3:500),rep('Excess',500-2))
colnames(Estimated1) <- c("X","Value","Label")
Spherical <- cbind.data.frame(seq(3:500),eq2(3:500),rep('Spherical',500-2))
colnames(Spherical) <- c("X","Value","Label")
Cylindrical <-  cbind.data.frame(seq(3:500),eq3(3:500),rep('Cylindrical',500-2))
colnames(Cylindrical) <- c("X","Value","Label")

# Combine all three into a single dataframe
attenuation.df <- rbind.data.frame(Estimated1,Spherical)

# Convert to reasonable starting value
attenuation.df$Value <-  attenuation.df$Value + 120.5759

# Plot the results
ggplot(data = attenuation.df,aes(x=X, y=Value,group=Label, colour=Label,linetype=Label))+
  geom_line() +theme_bw() + scale_color_manual(values = c("red","black","darkgray"))+
  theme(legend.title = element_blank())+ 
  scale_linetype_manual(values=c( "solid","twodash", "dotted"))+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
  xlab("Distance from source (m)") + ylab("Amplitude (dB)")+
  geom_hline(yintercept=noise.val,color='blue')+
  # ylim(-120,0)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



