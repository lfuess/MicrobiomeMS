##path analysis prep##
library(tidyr)
library(data.table)
library(dplyr)
  ##generate list of candidate genes
  genes = read.csv("Spartobacteria_Tau_Matrix.csv")
  genes = as.data.frame(t(genes))
  genes = setDT(genes, keep.rownames = TRUE)[]
  genes = genes[,c(1)]
  colnames(genes)[1] <- "Gene"
  genes = head(genes,-1)
  ##merge in expression data##
  reads = read.csv("normalizedreads_TagSeq.csv", check.names = FALSE)
  colnames(reads)[1] <- "Gene"
  reads$Gene = gsub("\\..*","",reads$Gene)
  genesexp = merge(genes,reads, by= "Gene", all.x = TRUE)  
  ##set row names##
  Merged2 <- data.frame(genesexp, row.names = 1, check.names = FALSE)
  ##transpose##
  expdat <- as.data.frame(t(Merged2))
  ##print rownames as column##
  expdat = setDT(expdat, keep.rownames = TRUE)[]
  ##rename that column##
  colnames(expdat)[1] <- "Sample"
  
##Now we go through the hassle of getting the meta-data##
  ##alpha diversity##
  T2 = read.csv("FamilyPropData.csv", header = TRUE)
  names(T2)
  FullMD = T2[c(1,265)]
  ##RNAseq metadata##
  MetaDataExp=read.csv("Masterdata_TagSeq_SDH_23July2018.csv", header=TRUE)
  names(MetaDataExp)  
  MetaDataExp=MetaDataExp[c(2,27,30:31,41)]  
  ##Merge everything together##
  FullExp=merge(expdat, MetaDataExp,
                by.x = "Sample", by.y = "sample_ID", all.x = FALSE,
                sort = TRUE,  no.dups = FALSE,)
  FinalData = merge(FullExp, FullMD,
                    by.x = "Sample", by.y = "SampleRNA", no.dups = FALSE)
  names(FinalData)
  FinalData = FinalData[,c(1, 490:494, 2:489)]
  names(FinalData)
  ##create a new column for log(mass)
  FinalData$LogMass = log(FinalData$weight)
  #create a new column transforming the proportion
  FinalData$Spart_T = asin(sqrt(FinalData$Spartobacteria_unclassified_Prop))
  names(FinalData)
  FinalData = FinalData[,c(1:4,495:496,7:494)]
  FinalData$CrossDir <- as.character(FinalData$CrossDir)
  FinalData$Sex <- as.character(FinalData$Sex)
  ##sub out values##
  FinalData$CrossDir[FinalData$CrossDir=="RBC"] <- 0.75
  FinalData$CrossDir[FinalData$CrossDir=="F2"] <- 0.50
  FinalData$CrossDir[FinalData$CrossDir=="GBC"] <- 0.25
  FinalData$Sex[FinalData$Sex=="M"] <- 1
  FinalData$Sex[FinalData$Sex=="F"] <- 0
  FinalData$worm_present[FinalData$worm_present=="FALSE"] <- 0
  FinalData$worm_present[FinalData$worm_present=="TRUE"] <- 1
  names(FinalData)
  colnames(FinalData)[1]="Sample"
  ##remove all rows with NAs-reduces our sample size down to 382##
  Final = FinalData[complete.cases(FinalData), ]
  ##write the Matrix for keeps##
  write.csv(Final, file = "PathMatrix_Spartobacteria_T.csv", row.names = FALSE)
  
  ##Ok move on to general data formatting in preparation for our loop##
  ##format PathMatrix##
  dat <- read.csv("PathMatrix_Spartobacteria_T.csv", header = T)
  dat2 <- dat[,-1]
  rownames(dat2) <- dat[,1]
  names(dat2)
  ##now we want to iterate over every column in the matrix, creating submatrixes and running stats##
  output = file("SEM_code_Spartobacteria.R")
  cat("", file=output, append=FALSE)
  cat("library(sem)", "\n", "library(semPlot)", "\n","\n","\n",
      'dat <- read.csv("PathMatrix_Spartobacteria_T.csv", header = T)', "\n",
      'dat2 <- dat[,-1]', "\n",
      'rownames(dat2) <- dat[,1]', "\n","\n","\n",
      file = output, append=TRUE)
  for (i in 6:493) {
    x=colnames(dat2)[i] 
    output = file("SEM_code_Spartobacteria.R", open = "a")
    cat((paste0("dat", i, "= dat2[,c(1:5,",i,")]","\n",
                "S <-cov(dat", i, ")", "\n","N <- dim(dat", i, ")[1]", "\n",
                "RAM <- specifyModel()", "\n",
                "CrossDir -> worm_present, CI, NA", "\n",
                "CrossDir -> LogMass,CM, NA", "\n",
                "CrossDir -> Spart_T,CP, NA", "\n",
                "CrossDir <-> CrossDir, C, NA", "\n",
                "CrossDir -> ", x, " ,CE, NA", "\n",
                "Sex -> ", x, " ,SE, NA", "\n",
                "Sex -> Spart_T,SP, NA", "\n",
                "Sex -> LogMass,SM, NA", "\n",
                "Sex <-> Sex, S, NA", "\n",
                "LogMass -> ", x, " ,ME, NA", "\n",
                "LogMass -> Spart_T,MP, NA", "\n",
                "LogMass <-> LogMass, M, NA", "\n",
                x, " <-> Spart_T, EP, NA", "\n",
                x, " <-> worm_present, EI, NA", "\n",
                x, " <-> ", x, " , E, NA", "\n",
                "worm_present <-> Spart_T, IP, NA", "\n",
                "worm_present <-> worm_present, I, NA", "\n",
                "Spart_T <-> Spart_T, P, NA", "\n","\n","\n",
                "sem.out <- sem(RAM, S, N)", "\n", "x=summary(sem.out)", "\n",
                "y=standardizedCoefficients(sem.out)", "\n",
                'capture.output(x, y, file = ("semoutput_' , x , '.txt"))', '\n',
                "\n","\n","\n", "##Onto Next##", "\n")), file=output)
    close(output)
  }

  