
# Prop loss curves -----------------------------------------------------------------

# Prep Rungan Data
observed.prop.lossRunganSubset <- read.csv('observed.prop.lossRunganAugust1.csv')


# Gibbon prop loss --------------------------------------------------------

observed.prop.lossRunganGibbons <- subset(observed.prop.lossRunganSubset,Call.category=="Hfunstart" |Call.category=="Hfuntrill" |
                                           Call.category=="Halbstart" |Call.category=="Halbpeak" )


observed.prop.lossRunganGibbons$Species <- 
  recode(observed.prop.lossRunganGibbons$Call.category, Hfunstart = "NGreyGibbon",
         Hfuntrill = "NGreyGibbon",Halbstart='WhiteBeardGibbon',Halbend='WhiteBeardGibbon',Halbpeak='WhiteBeardGibbon')

observed.prop.lossRunganGibbons <- subset(observed.prop.lossRunganGibbons,distance<500)
uniquegibbons <- unique(observed.prop.lossRunganGibbons$Species)
gibbondB <- 100


for(d in 1:length(uniquegibbons)){
  
  observed.prop.lossRungangibbonstemp <- subset(observed.prop.lossRunganGibbons,Species==uniquegibbons[d])
  
  gibbonpropdBRungan <- data.frame()
  gibbonNoisedf <- list()
  for(e in 1:nrow(observed.prop.lossRungangibbonstemp)){
    
    temp.magicx.row <- observed.prop.lossRungangibbonstemp[e,]
    
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
    
    gibbonpropdBRungan <- rbind.data.frame(gibbonpropdBRungan,newtemprow)
  }
  
  gibbonpropdBRungan$Value <- as.numeric(gibbonpropdBRungan$Value)
  gibbonpropdBRungan$X <- as.numeric(gibbonpropdBRungan$X)
  
  gibbonpropdBRunganCI <- gibbonpropdBRungan %>%
    group_by(X) %>%
    summarise(meangibbon = mean(Value, na.rm = TRUE),
              sdgibbon = sd(Value, na.rm = TRUE),
              ngibbon = n()) %>%
    mutate(segibbon = sdgibbon / sqrt(ngibbon),
           lower.cigibbon = meangibbon - qt(1 - (0.05 / 2), ngibbon - 1) * segibbon,
           upper.cigibbon = meangibbon + qt(1 - (0.05 / 2), ngibbon - 1) * segibbon)
  
  
  noise.val <- median(unlist(gibbonNoisedf))
  
  gibbonplot <- ggplot(gibbonpropdBRunganCI, aes(x = X, y = meangibbon, group = 1)) + 
    geom_line(col='red') + 
    geom_ribbon(aes(ymin = lower.cigibbon, ymax = upper.cigibbon), alpha = 0.25)+
    geom_hline(yintercept=noise.val,linetype="dashed", color = "black")+
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
    xlab("Distance from source (m)") + ylab("Amplitude (dB)")+theme_bw()+ggtitle(paste('Rungan',uniquegibbons[d], gibbondB, 'dB @ 1 m Source level'))+ylim(25,125)
  print(gibbonplot)
}


# Orangutan prop loss --------------------------------------------------------

observed.prop.lossRunganorangutans <- subset(observed.prop.lossRungan,Call.category=="Pmor" |Call.category=="PwurP" |
                                            Call.category=="PwurS" )

observed.prop.lossRunganorangutans$Species <- 
  recode(observed.prop.lossRunganorangutans$Call.category, Pmor='OrangSabah',PwurP='OrangCKalimantan',
         PwurS='OrangCKalimantan')

uniqueorangutans <- unique(observed.prop.lossRunganorangutans$Species)
orangutandB <- 100


for(d in 1:length(uniqueorangutans)){
  
  observed.prop.lossRunganorangutanstemp <- subset(observed.prop.lossRunganorangutans,Species==uniqueorangutans[d])
  
  orangutanpropdBRungan <- data.frame()
  orangutanNoisedf <- list()
  for(e in 1:nrow(observed.prop.lossRunganorangutanstemp)){
    
    temp.magicx.row <- observed.prop.lossRunganorangutanstemp[e,]
    
    orangutanNoisedf[[e]] <- temp.magicx.row$noise.level
    
    playback.line.1 <- temp.magicx.row$magic.x
    
    # Set the equations for observed, spherical and cylindrical spreading
    eq1 <- function(x){ playback.line.1*log10(x)}
    
    # Create a series of points based on the above equations
    Estimated1 <- cbind.data.frame(seq(1:500),eq1(1:500),rep('Estimated',500))
    colnames(Estimated1) <- c("X","Value","Label")
    
    Estimated1$X <- Estimated1$X -1
    
    scaledB <- orangutandB - Estimated1$Value[2] 
    Estimated1$Value <- scaledB+Estimated1$Value
    newtemprow <- rbind.data.frame(Estimated1,temp.magicx.row$Call.category)
    
    orangutanpropdBRungan <- rbind.data.frame(orangutanpropdBRungan,newtemprow)
  }
  
  orangutanpropdBRungan$Value <- as.numeric(orangutanpropdBRungan$Value)
  orangutanpropdBRungan$X <- as.numeric(orangutanpropdBRungan$X)
  
  orangutanpropdBRunganCI <- orangutanpropdBRungan %>%
    group_by(X) %>%
    summarise(meanorangutan = mean(Value, na.rm = TRUE),
              sdorangutan = sd(Value, na.rm = TRUE),
              norangutan = n()) %>%
    mutate(seorangutan = sdorangutan / sqrt(norangutan),
           lower.ciorangutan = meanorangutan - qt(1 - (0.05 / 2), norangutan - 1) * seorangutan,
           upper.ciorangutan = meanorangutan + qt(1 - (0.05 / 2), norangutan - 1) * seorangutan)
  
  
  noise.val <- median(unlist(orangutanNoisedf))
  
  orangutanplot <- ggplot(orangutanpropdBRunganCI, aes(x = X, y = meanorangutan, group = 1)) + 
    geom_line(col='red') + 
    geom_ribbon(aes(ymin = lower.ciorangutan, ymax = upper.ciorangutan), alpha = 0.25)+
    geom_hline(yintercept=noise.val,linetype="dashed", color = "black")+
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
    xlab("Distance from source (m)") + ylab("Amplitude (dB)")+theme_bw()+ggtitle(paste('Rungan',uniqueorangutans[d]))
  print(orangutanplot)
}
