MaliauDF <- read.csv("BackgroundNoiseRemovedMaliauAugust9AdaptiveNoise.csv")

TempSpectrograms <- list.files('/Users/denaclink/Desktop/RStudio Projects/Propogation-Loss-2022/NoiseSpectrograms5sec25ms',
                               pattern = '.pdf',
           full.names = T)

TempSpectrogramsShort <- list.files('/Users/denaclink/Desktop/RStudio Projects/Propogation-Loss-2022/NoiseSpectrograms5sec25ms',
                                    pattern = '.pdf',full.names = F)

Shortnames <- str_split_fixed(TempSpectrogramsShort,' ',n=3)[,2]
sort <- str_split_fixed(TempSpectrogramsShort,' ',n=4)[,3]

TempSpectrograms <- cbind.data.frame(TempSpectrograms,Shortnames,sort)
TempSpectrograms$sort <- as.numeric(TempSpectrograms$sort)

SelectionIDsMaliau$SelectionID <- seq(1,15,1)
  
UpdatedTempSpectrograms <- data.frame()

IndexShortnames <- unique(Shortnames)
for(z in 1:length(IndexShortnames) ){
  TempSubset <- subset(TempSpectrograms,Shortnames==IndexShortnames[z])
  TempSubset <- TempSubset[order(TempSubset$sort),]
  UpdatedTempSpectrograms <- rbind.data.frame(UpdatedTempSpectrograms,TempSubset)
}


# Stop at number 1588
UpdateMaliauDF <- read.csv('UpdateMaliauDFValidatedAugust22Noise5sec25ms.csv')

#UpdateMaliauDF <- data.frame()
for(a in 1342:nrow(UpdatedTempSpectrograms)) {
  TempSpec <- UpdatedTempSpectrograms[a,]
  TempIndex <- subset(SelectionIDsMaliau,SelectionID==TempSpec$sort)
  TempRow <-  subset(MaliauDF,file.name==TempSpec$Shortnames &
                       Sound.Type==TempIndex$Sound.Type)

 figure1.png <- magick::image_trim(magick::image_read(TempSpec$TempSpectrograms))
 print(figure1.png)
 print(TempRow[,c('Sound.Type','file.name')])
 target.signal <- readline(prompt = "For analysis? ")
 if(target.signal== 'break'){
   print('End loop')
   print(paste(a, 'out of', nrow(MaliauDF)))
   break
   
 }
 TempRow$target.signal <- target.signal
 UpdateMaliauDF <- rbind.data.frame(UpdateMaliauDF, TempRow)
 print(paste(a, 'out of', nrow(MaliauDF)))
 write.csv(UpdateMaliauDF,'UpdateMaliauDFValidatedAugust22Noise5sec25ms.csv',row.names = F)
}


