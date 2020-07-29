##create an R script to parse out WGCNA output and generate files for GOMWU and the likes##
##first we load our package##
library(tidyr)
##input WGCNA output##
data = read.csv("WGCNA_Diversity_Corrs.csv")
names(data)
##pull out just the rows we need##
data = data[c(1:2,5,7,9,11,13,17,19,21,23,25)]
##import all your genes##
genenames = read.csv("normalizedreads_TagSeq.csv", check.names = FALSE)
genenames=genenames[c(1:2)]

##ok start with the magenta module##
  ##pull out just the rows we want##
  magenta = data[data$moduleColor %in% c("magenta"), ] 
  magenta = magenta[c(1,7)]
  ##write this out to a csv##
  write.csv(magenta, "magenta_transcripts.csv", row.names =FALSE)
  ##create the GOMWU input##

  ##and merge##
  Merged = merge(magenta, genenames, by.x = "locus", by.y = "Sample", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "magenta_kme_GOMWU.csv", row.names = FALSE)
  
##do the same for the yellow module##
  ##pull out just the rows we want##
  yellow = data[data$moduleColor %in% c("yellow"), ] 
  names(yellow)
  yellow = yellow[c(1,6)]
  ##write this out to a csv##
  write.csv(yellow, "yellow_transcripts.csv", row.names =FALSE)
  ##create the GOMWU input##
  Merged = merge(yellow, genenames, by.x = "locus", by.y = "Sample", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "yellow_kme_GOMWU.csv", row.names = FALSE)

##the rest are for family analyses##
  ##pink module##
  pink = data[data$moduleColor %in% c("pink"), ] 
  names(pink)
  pink = pink[c(1,8)]
  ##write this out to a csv##
  write.csv(pink, "pink_transcripts.csv", row.names =FALSE)
  ##create the GOMWU input##
  Merged = merge(pink, genenames, by.x = "locus", by.y = "Sample", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "pink_kme_GOMWU.csv", row.names = FALSE)
  
  
  ##blue module##
  blue = data[data$moduleColor %in% c("blue"), ] 
  names(blue)
  blue = blue[c(1,9)]
  ##write this out to a csv##
  write.csv(blue, "blue_transcripts.csv", row.names =FALSE)
  ##create the GOMWU input##
  Merged = merge(blue, genenames, by.x = "locus", by.y = "Sample", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "blue_kme_GOMWU.csv", row.names = FALSE)
  
  ##black module##
  black = data[data$moduleColor %in% c("black"), ] 
  names(black)
  black = black[c(1,7)]
  ##write this out to a csv##
  write.csv(black, "black_transcripts.csv", row.names =FALSE)
  ##create the GOMWU input##
  Merged = merge(black, genenames, by.x = "locus", by.y = "Sample", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "black_kme_GOMWU.csv", row.names = FALSE)
  
  
  ##brown module##
  brown = data[data$moduleColor %in% c("brown"), ] 
  names(brown)
  brown = brown[c(1,4)]
  ##write this out to a csv##
  write.csv(brown, "brown_transcripts.csv", row.names =FALSE)
  ##create the GOMWU input##
  Merged = merge(brown, genenames, by.x = "locus", by.y = "Sample", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "brown_kme_GOMWU.csv", row.names = FALSE)

  
  
  ##green module##
  green = data[data$moduleColor %in% c("green"), ] 
  names(green)
  green = green[c(1,5)]
  ##write this out to a csv##
  write.csv(green, "green_transcripts.csv", row.names =FALSE)
  ##create the GOMWU input##
  Merged = merge(green, genenames, by.x = "locus", by.y = "Sample", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "green_kme_GOMWU.csv", row.names = FALSE)  
  
  
  
  ##purple module##
  purple = data[data$moduleColor %in% c("purple"), ] 
  names(purple)
  purple = purple[c(1,10)]
  ##write this out to a csv##
  write.csv(purple, "purple_transcripts.csv", row.names =FALSE)
  ##create the GOMWU input##
  Merged = merge(purple, genenames, by.x = "locus", by.y = "Sample", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "purple_kme_GOMWU.csv", row.names = FALSE) 
  
  
  ##turquoise module##
  turquoise = data[data$moduleColor %in% c("turquoise"), ] 
  names(turquoise)
  turquoise = turquoise[c(1,6)]
  ##write this out to a csv##
  write.csv(turquoise, "turquoise_transcripts.csv", row.names =FALSE)
  ##create the GOMWU input##
  Merged = merge(turquoise, genenames, by.x = "locus", by.y = "Sample", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "turquoise_kme_GOMWU.csv", row.names = FALSE) 
  
  ##red module##
  red = data[data$moduleColor %in% c("red"), ] 
  names(red)
  red = red[c(1,8)]
  ##write this out to a csv##
  write.csv(red, "red_transcripts.csv", row.names =FALSE)
  ##create the GOMWU input##
  Merged = merge(red, genenames, by.x = "locus", by.y = "Sample", all.y = TRUE)
  ##replace NA with 0, to indicate to GOMWU genes that aren't in our list of interest##
  Merged[is.na(Merged)] <- 0
  dim(Merged)
  names(Merged)
  ##grab just the gene names and the stat##
  Merged=Merged[c(1:2)]
  ##and write it out##
  write.csv(Merged, "red_kme_GOMWU.csv", row.names = FALSE) 
  
  