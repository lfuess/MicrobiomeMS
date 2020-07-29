##Generating a matrix of genes that are significantly correlated to Diversity for Tau value calculation##

##load in packages##
library(dplyr)
library(data.table)

##First let's load in a list of the genes which are significantlly correlated to Alpha diversity##
Genes = read.csv("SigCorrResults_withPs.csv")
colnames(Genes)[1] <- "Sample"
  ##load in your readcounts##
  Reads = read.csv("normalizedreads_TagSeq.csv", check.names = FALSE)
  ##remove the part of the transcript names after the period##
  Reads$Sample = gsub("\\..*","", Reads$Sample)
  ##merge them together##
  FullMatrix=merge(Genes, Reads,
                by = "Sample", all.x = FALSE,
                sort = TRUE,  no.dups = FALSE,)
##remove unnecessary columns##
  dim(FullMatrix)
  FullMatrix= FullMatrix[,c(1,3:389)]  
  FullMatrix2 <- FullMatrix[,-1]
  rownames(FullMatrix2) <- FullMatrix[,1]
##Merge in diversity data##
  FullMatrix2 = as.data.frame(t(FullMatrix2))
  FullMatrix2 = setDT(FullMatrix2, keep.rownames = TRUE)[]  
  ##load diversity data##
  diversity = read.csv("DiversityData2000.csv")
  diversity = diversity[c(1,3)]
  ##merge##
  Final=merge(FullMatrix2, diversity,
                   by.x = "rn", by.y = "SampleRNA", all.x = TRUE,
                   sort = TRUE,  no.dups = FALSE,)
  ##remove sample names##
  Final_all=Final[,c(2:1931)]
  ##write to csv##
  write.csv(Final_all,"Diversity_Tau_Matrix.csv", row.names = FALSE)
  
  
  
  
  
  ##We can also split this matrix up by infection status to see how cestode infection affects these relationships##
  
  ##so we need to merge in the info we want##
  MetaDataExp=read.csv("Masterdata_TagSeq_SDH_23July2018.csv", header=TRUE)
  names(MetaDataExp)  
  InfectionData=MetaDataExp[,c(2,31)]
  ##merge with previous matrix##
  Final_Infect=merge(Final, InfectionData,
              by.x = "rn", by.y = "sample_ID", all.x = TRUE,
              sort = TRUE,  no.dups = FALSE,)
  ##split it up##
  Final_True=Final_Infect[Final_Infect$worm_present == TRUE, ] 
  ##remove sample names##
  Final_True=Final_True[,c(2:1931)]
  ##write to csv##
  write.csv(Final_True,"Diversity_Infected_Matrix.csv", row.names = FALSE)
  ##Uninfected##
  Final_False=Final_Infect[Final_Infect$worm_present == FALSE, ] 
  ##remove sample names##
  Final_False=Final_False[,c(2:1931)]
  ##write to csv##
  write.csv(Final_False,"Diversity_Uninfected_Matrix.csv", row.names = FALSE)
  