
# Prop loss curves -----------------------------------------------------------------

# Prep Rungan Data
#observed.prop.lossRunganSubset <- read.csv('observed.prop.lossRunganAugust1.csv')
observed.prop.lossRunganGibbons<- droplevels(subset(observed.prop.lossRunganOutRM,
                                                     Species != "OrangSabah" & Species != "OrangCKalimantan"))#read.csv('observed.prop.lossRunganAugust4.csv')


# gibbonRungan prop loss --------------------------------------------------------
# Set the equations for observed, spherical and cylindrical spreading
eq1 <- function(x){ median(observed.prop.lossRunganGibbons$magic.x) *log10(x)}
eq2 <- function(x){ -20*log10(x)}
eq3 <- function(x){ -10*log10(x)}

uniquegibbonRungans <- unique(observed.prop.lossRunganGibbons$Species)
gibbonRungandB <- 107.1
gibbonRungandB <- 113.9

gibbonRunganplotlist <- list()

for(d in 1:length(uniquegibbonRungans)){
  
  observed.prop.lossRunganGibbonstemp <- subset(observed.prop.lossRunganGibbons,Species==uniquegibbonRungans[d])
  
  gibbonRunganpropdBRungan <- data.frame()
  gibbonRunganNoisedf <- list()
  # gibbonRunganNoisedfMorning <- list()
  # gibbonRunganNoisedfAfternoon <- list()
  
  for(e in 1:nrow(observed.prop.lossRunganGibbonstemp)){
    
    temp.magicx.row <- observed.prop.lossRunganGibbonstemp[e,]
    
    gibbonRunganNoisedf[[e]] <- temp.magicx.row$noise.level
    
    playback.line.1 <- temp.magicx.row$magic.x
    
    # Set the equations for observed, spherical and cylindrical spreading
    eq1 <- function(x){ playback.line.1*log10(x)}
    
    # Create a series of points based on the above equations
    Estimated1 <- cbind.data.frame(seq(1:500),eq1(1:500),rep('Estimated',500))
    colnames(Estimated1) <- c("X","Value","Label")
    
    Estimated1$X <- Estimated1$X -1
    
    scaledB <- gibbonRungandB - Estimated1$Value[2] 
    Estimated1$Value <- scaledB+Estimated1$Value
    newtemprow <- rbind.data.frame(Estimated1,temp.magicx.row$Call.category)
    
    gibbonRunganpropdBRungan <- rbind.data.frame(gibbonRunganpropdBRungan,newtemprow)
  }
  
  gibbonRunganpropdBRungan$Value <- as.numeric(gibbonRunganpropdBRungan$Value)
  gibbonRunganpropdBRungan$X <- as.numeric(gibbonRunganpropdBRungan$X)
  
  gibbonRunganpropdBRungan <- 
    subset(gibbonRunganpropdBRungan,Value > 10 )
  
  gibbonRunganpropdBRunganCI <- gibbonRunganpropdBRungan %>%
    group_by(X) %>%
    dplyr::summarise(meangibbonRungan = mean(Value, na.rm = TRUE),
                     sdgibbonRungan = sd(Value, na.rm = TRUE),
                     ngibbonRungan = n()) %>%
    mutate(segibbonRungan = sdgibbonRungan / sqrt(ngibbonRungan),
           lower.cigibbonRungan = meangibbonRungan - qt(1 - (0.05 / 2), ngibbonRungan - 1) * segibbonRungan,
           upper.cigibbonRungan = meangibbonRungan + qt(1 - (0.05 / 2), ngibbonRungan - 1) * segibbonRungan)
  
  
  noise.val <- median(unlist(gibbonRunganNoisedf))
  senoise = sd(unlist(gibbonRunganNoisedf)) / sqrt(length(gibbonRunganNoisedf))
  lower.cinoise = noise.val - qt(1 - (0.05 / 2), length(gibbonRunganNoisedf) - 1) * senoise
  upper.cinoise = noise.val + qt(1 - (0.05 / 2), length(gibbonRunganNoisedf) - 1) * senoise
  
  gibbonRunganplot <- ggplot(gibbonRunganpropdBRunganCI, aes(x = X, y = meangibbonRungan, group = 1)) + 
    geom_line(col='red') + 
    geom_ribbon(aes(ymin = lower.cigibbonRungan, ymax = upper.cigibbonRungan), alpha = 0.25)+
    #geom_hline(yintercept=noise.val,linetype="dashed", color = "black")+
    geom_ribbon(aes(ymin = lower.cinoise, ymax = upper.cinoise), alpha = 0.25)+
    
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
    xlab("Distance from source (m)") + ylab("Amplitude (dB)")+theme_bw()+ggtitle(paste('Mungku Baru',uniquegibbonRungans[d], gibbonRungandB, 'dB @ 1 m Source level'))+
    ylim(25,gibbonRungandB)+xlim(1,500)
  
  gibbonRunganplotlist[[d]] <- gibbonRunganplot
}


cowplot::plot_grid(gibbonRunganplotlist[[1]],gibbonRunganplotlist[[2]])


# Orangutan prop loss --------------------------------------------------------



# Prep Rungan Data
#observed.prop.lossRunganSubset <- read.csv('observed.prop.lossRunganAugust1.csv')
observed.prop.lossRunganOrangutans<- droplevels(subset(observed.prop.lossRunganOutRM ,
                                                    Species == "OrangSabah" | Species == "OrangCKalimantan"))#read.csv('observed.prop.lossRunganAugust4.csv')


# OrangutanRungan prop loss --------------------------------------------------------
# Set the equations for observed, spherical and cylindrical spreading
eq1 <- function(x){ median(observed.prop.lossRunganOrangutans$magic.x) *log10(x)}
eq2 <- function(x){ -20*log10(x)}
eq3 <- function(x){ -10*log10(x)}

uniqueOrangutanRungans <- unique(observed.prop.lossRunganOrangutans$Species)
OrangutanRungandB <- 107.1
OrangutanRungandB <- 105

OrangutanRunganplotlist <- list()

for(d in 1:length(uniqueOrangutanRungans)){
  
  observed.prop.lossRunganOrangutanstemp <- subset(observed.prop.lossRunganOrangutans,Species==uniqueOrangutanRungans[d])
  
  OrangutanRunganpropdBRungan <- data.frame()
  OrangutanRunganNoisedf <- list()
  # OrangutanRunganNoisedfMorning <- list()
  # OrangutanRunganNoisedfAfternoon <- list()
  
  for(e in 1:nrow(observed.prop.lossRunganOrangutanstemp)){
    
    temp.magicx.row <- observed.prop.lossRunganOrangutanstemp[e,]
    
    OrangutanRunganNoisedf[[e]] <- temp.magicx.row$noise.level
    
    playback.line.1 <- temp.magicx.row$magic.x
    
    # Set the equations for observed, spherical and cylindrical spreading
    eq1 <- function(x){ playback.line.1*log10(x)}
    
    # Create a series of points based on the above equations
    Estimated1 <- cbind.data.frame(seq(1:500),eq1(1:500),rep('Estimated',500))
    colnames(Estimated1) <- c("X","Value","Label")
    
    Estimated1$X <- Estimated1$X -1
    
    scaledB <- OrangutanRungandB - Estimated1$Value[2] 
    Estimated1$Value <- scaledB+Estimated1$Value
    newtemprow <- rbind.data.frame(Estimated1,temp.magicx.row$Call.category)
    
    OrangutanRunganpropdBRungan <- rbind.data.frame(OrangutanRunganpropdBRungan,newtemprow)
  }
  
  OrangutanRunganpropdBRungan$Value <- as.numeric(OrangutanRunganpropdBRungan$Value)
  OrangutanRunganpropdBRungan$X <- as.numeric(OrangutanRunganpropdBRungan$X)
  
  OrangutanRunganpropdBRungan <- 
    subset(OrangutanRunganpropdBRungan,Value > 10 )
  
  
  OrangutanRunganpropdBRunganCI <- OrangutanRunganpropdBRungan %>%
    group_by(X) %>%
    dplyr::summarise(meanOrangutanRungan = mean(Value, na.rm = TRUE),
                     sdOrangutanRungan = sd(Value, na.rm = TRUE),
                     nOrangutanRungan = n()) %>%
    mutate(seOrangutanRungan = sdOrangutanRungan / sqrt(nOrangutanRungan),
           lower.ciOrangutanRungan = meanOrangutanRungan - qt(1 - (0.05 / 2), nOrangutanRungan - 1) * seOrangutanRungan,
           upper.ciOrangutanRungan = meanOrangutanRungan + qt(1 - (0.05 / 2), nOrangutanRungan - 1) * seOrangutanRungan)
  
  
  noise.val <- median(unlist(OrangutanRunganNoisedf))
  print(noise.val)
  senoise = sd(unlist(OrangutanRunganNoisedf)) / sqrt(length(OrangutanRunganNoisedf))
  lower.cinoise = noise.val - qt(1 - (0.05 / 2), length(OrangutanRunganNoisedf) - 1) * senoise
  upper.cinoise = noise.val + qt(1 - (0.05 / 2), length(OrangutanRunganNoisedf) - 1) * senoise
  print(lower.cinoise)
  print(upper.cinoise)
  OrangutanRunganplot <- ggplot(OrangutanRunganpropdBRunganCI, aes(x = X, y = meanOrangutanRungan, group = 1)) + 
    geom_line(col='red') + 
    geom_ribbon(aes(ymin = lower.ciOrangutanRungan, ymax = upper.ciOrangutanRungan), alpha = 0.25)+
    #geom_hline(yintercept=noise.val,linetype="dashed", color = "black")+
    geom_ribbon(aes(ymin = lower.cinoise, ymax = upper.cinoise), alpha = 0.25)+
    
    theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
    xlab("Distance from source (m)") + ylab("Amplitude (dB)")+theme_bw()+ggtitle(paste('Mungku Baru',uniqueOrangutanRungans[d], OrangutanRungandB, 'dB @ 1 m'))+
    ylim(25,OrangutanRungandB)
  
  OrangutanRunganplotlist[[d]] <- OrangutanRunganplot
}


cowplot::plot_grid(gibbonRunganplotlist[[1]],gibbonRunganplotlist[[2]],
                   OrangutanRunganplotlist[[1]],OrangutanRunganplotlist[[2]])
