
# Prop loss curves -----------------------------------------------------------------
# V2 add by habitat

# Prep Rungan Data
#observed.prop.lossRunganSubset <- read.csv('observed.prop.lossRunganAugust1.csv')
observed.prop.lossRunganGibbons<- droplevels(subset(observed.prop.lossRunganOutRM,
                                                     Species != "OrangSabah" & Species != "OrangCKalimantan"))#read.csv('observed.prop.lossRunganAugust4.csv')

# GibbonRungan prop loss --------------------------------------------------------
# Set the equations for observed, spherical and cylindrical spreading
eq1 <- function(x){ median(observed.prop.lossRunganGibbons$magic.x) *log10(x)}
eq2 <- function(x){ -20*log10(x)}
eq3 <- function(x){ -10*log10(x)}

uniqueGibbonRungans <- unique(observed.prop.lossRunganGibbons$Species)
GibbonRungandB <- 113.9



GibbonRunganplotlistComb <- list()

for(d in 1:length(uniqueGibbonRungans)){
  
  observed.prop.lossRunganGibbonstemp <- subset(observed.prop.lossRunganGibbons,Species==uniqueGibbonRungans[d])
  
  unique.hab <- unique(observed.prop.lossRunganGibbonstemp$habitat)
  GibbonRunganplotlist<- list()
  for(a in 1:length(unique.hab)){
    observed.prop.lossRunganGibbonshab <-   
      subset(observed.prop.lossRunganGibbonstemp,habitat==unique.hab[a])
    GibbonRunganpropdBRungan <- data.frame()

    GibbonRunganNoisedf <- list()
    # GibbonRunganNoisedfMorning <- list()
    # GibbonRunganNoisedfAfternoon <- list()
    
    for(e in 1:nrow(observed.prop.lossRunganGibbonshab)){
      
      temp.magicx.row <- observed.prop.lossRunganGibbonshab[e,]
      
      GibbonRunganNoisedf[[e]] <- temp.magicx.row$noise.level
      
      playback.line.1 <- temp.magicx.row$magic.x
      
      # Set the equations for observed, spherical and cylindrical spreading
      eq1 <- function(x){ playback.line.1*log10(x)}
      
      # Create a series of points based on the above equations
      Estimated1 <- cbind.data.frame(seq(1:500),eq1(1:500),rep('Estimated',500))
      colnames(Estimated1) <- c("X","Value","Label")
      
      Estimated1$X <- Estimated1$X -1
      
      scaledB <- GibbonRungandB - Estimated1$Value[2] 
      Estimated1$Value <- scaledB+Estimated1$Value
      newtemprow <- rbind.data.frame(Estimated1,temp.magicx.row$Call.category)
      
      GibbonRunganpropdBRungan <- rbind.data.frame(GibbonRunganpropdBRungan,newtemprow)
    }
    
    GibbonRunganpropdBRungan$Value <- as.numeric(GibbonRunganpropdBRungan$Value)
    GibbonRunganpropdBRungan$X <- as.numeric(GibbonRunganpropdBRungan$X)
    
    GibbonRunganpropdBRungan <- 
      subset(GibbonRunganpropdBRungan,Value > 10 )
    
    GibbonRunganpropdBRunganCI <- GibbonRunganpropdBRungan %>%
      group_by(X) %>%
      dplyr::summarise(meanGibbonRungan = mean(Value, na.rm = TRUE),
                       sdGibbonRungan = sd(Value, na.rm = TRUE),
                       nGibbonRungan = n()) %>%
      mutate(seGibbonRungan = sdGibbonRungan / sqrt(nGibbonRungan),
             lower.ciGibbonRungan = meanGibbonRungan - qt(1 - (0.05 / 2), nGibbonRungan - 1) * seGibbonRungan,
             upper.ciGibbonRungan = meanGibbonRungan + qt(1 - (0.05 / 2), nGibbonRungan - 1) * seGibbonRungan)
    
    
    noise.val <- median(unlist(GibbonRunganNoisedf))
    senoise = sd(unlist(GibbonRunganNoisedf)) / sqrt(length(GibbonRunganNoisedf))
    GibbonRunganpropdBRunganCI$lower.cinoise = noise.val - qt(1 - (0.05 / 2), length(GibbonRunganNoisedf) - 1) * senoise
    GibbonRunganpropdBRunganCI$upper.cinoise = noise.val + qt(1 - (0.05 / 2), length(GibbonRunganNoisedf) - 1) * senoise
    
    GibbonRunganplot <- ggplot(GibbonRunganpropdBRunganCI, aes(x = X, y = meanGibbonRungan, group = 1)) + 
      geom_line(col='red') + 
      geom_ribbon(aes(ymin = lower.ciGibbonRungan, ymax = upper.ciGibbonRungan), alpha = 0.25)+
      #geom_hline(yintercept=noise.val,linetype="dashed", color = "black")+
      geom_ribbon(aes(ymin = lower.cinoise, ymax = upper.cinoise), alpha = 0.5, color='blue')+
      
      theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
      xlab("Distance from source (m)") + ylab("Amplitude (dB)")+theme_bw()+ggtitle(paste(unique.hab[a],uniqueGibbonRungans[d]))+
      ylim(25,GibbonRungandB)+xlim(0,500)
    
    GibbonRunganplotlist[[a]] <- GibbonRunganplot
    print(GibbonRunganplot)
  }
  GibbonRunganplotlistComb[[d]] <-GibbonRunganplotlist
}


cowplot::plot_grid(gibbonplotlist[[1]],
                   GibbonRunganplotlistComb[[1]][[1]],
                   GibbonRunganplotlistComb[[1]][[2]],
                   GibbonRunganplotlistComb[[1]][[3]],
                   gibbonplotlist[[2]],
                   GibbonRunganplotlistComb[[2]][[1]],
                   GibbonRunganplotlistComb[[2]][[2]],
                   GibbonRunganplotlistComb[[2]][[3]],
                   nrow=2,labels=c('A','B','C','D','E','F','G','H'),
                   label_x = 0.9)

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
OrangutanRungandB <- 105



OrangutanRunganplotlistComb <- list()

for(d in 1:length(uniqueOrangutanRungans)){
  
  observed.prop.lossRunganOrangutanstemp <- subset(observed.prop.lossRunganOrangutans,Species==uniqueOrangutanRungans[d])
  
  unique.hab <- unique(observed.prop.lossRunganOrangutanstemp$habitat)
  GibbonRunganplotlist<- list()
   for(a in 1:length(unique.hab)){
    observed.prop.lossRunganOrangutanshab <-   
      subset(observed.prop.lossRunganOrangutanstemp,habitat==unique.hab[a])
    OrangutanRunganpropdBRungan <- data.frame()
    OrangutanRunganNoisedf <- list()
    # OrangutanRunganNoisedfMorning <- list()
    # OrangutanRunganNoisedfAfternoon <- list()
    
    for(e in 1:nrow(observed.prop.lossRunganOrangutanshab)){
      
      temp.magicx.row <- observed.prop.lossRunganOrangutanshab[e,]
      
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
    senoise = sd(unlist(OrangutanRunganNoisedf)) / sqrt(length(OrangutanRunganNoisedf))
    OrangutanRunganpropdBRunganCI$lower.cinoise = noise.val - qt(1 - (0.05 / 2), length(OrangutanRunganNoisedf) - 1) * senoise
    OrangutanRunganpropdBRunganCI$upper.cinoise = noise.val + qt(1 - (0.05 / 2), length(OrangutanRunganNoisedf) - 1) * senoise
    
    OrangutanRunganplot <- ggplot(OrangutanRunganpropdBRunganCI, aes(x = X, y = meanOrangutanRungan, group = 1)) + 
      geom_line(col='red') + 
      geom_ribbon(aes(ymin = lower.ciOrangutanRungan, ymax = upper.ciOrangutanRungan), alpha = 0.25)+
      #geom_hline(yintercept=noise.val,linetype="dashed", color = "black")+
      geom_ribbon(aes(ymin = lower.cinoise, ymax = upper.cinoise), alpha = 0.5, color='blue')+
      
      theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))+
      xlab("Distance from source (m)") + ylab("Amplitude (dB)")+theme_bw()+ggtitle(paste(unique.hab[a],uniqueOrangutanRungans[d]))+
      ylim(25,OrangutanRungandB)+xlim(0,500)
    
    OrangutanRunganplotlist[[a]] <- OrangutanRunganplot
    print(OrangutanRunganplot)
  }
  OrangutanRunganplotlistComb[[d]] <-OrangutanRunganplotlist
}


cowplot::plot_grid(OrangutanRunganplotlistComb[[1]][[1]],OrangutanRunganplotlistComb[[1]][[2]],
                   OrangutanRunganplotlistComb[[1]][[3]],
                   OrangutanRunganplotlistComb[[2]][[1]],OrangutanRunganplotlistComb[[2]][[2]],
                   OrangutanRunganplotlistComb[[2]][[3]])
