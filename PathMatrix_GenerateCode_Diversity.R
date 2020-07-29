##path analysis prep##
library(tidyr)
library(data.table)
library(dplyr)
  ##generate list of candidate genes
  genes = read.csv("Diversity_Tau_Matrix.csv")
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
  T2 = read.csv("alpha_2000.csv", header = TRUE)
  names(T2)
  FullMD = T2[c(3,4,6,12:13,18)]
  names(FullMD) = c("Sample","Date", "Number", "Sex", "worm_present", "Alpha_diversity_mean_2000")
  ##RNAseq metadata##
  MetaDataExp=read.csv("Masterdata_TagSeq_SDH_23July2018.csv", header=TRUE)
  names(MetaDataExp)  
  MetaDataExp=MetaDataExp[c(2,20,27,29:30,31,41)]  
  colnames(MetaDataExp)[2]="Date"
  ##Merge this with Data Matrix##
  FullExp=merge(expdat, MetaDataExp,
                by.x = "Sample", by.y = "sample_ID", all.x = FALSE,
                sort = TRUE,  no.dups = FALSE,)
  ##A few more housekeeping things to make the data mergable##
  ##fix the dates so they're consistent across data frames (FullMD has them listed as /201X)##
  FullExp <- as.data.frame(sapply(FullExp,gsub,pattern="/15",replacement="/2015"))
  FullExp <- as.data.frame(sapply(FullExp,gsub,pattern="9/14",replacement="9/2014"))
  FullExp <- as.data.frame(sapply(FullExp,gsub,pattern="0/14",replacement="0/2014"))
  FullExp <- as.data.frame(sapply(FullExp,gsub,pattern="7/14",replacement="7/2014"))
  FullExp <- as.data.frame(sapply(FullExp,gsub,pattern="3/14",replacement="3/2014"))
  FullExp <- as.data.frame(sapply(FullExp,gsub,pattern="4/14",replacement="4/2014"))
  FullExp <- as.data.frame(sapply(FullExp,gsub,pattern="2/14",replacement="2/2014"))
  FullMD <- as.data.frame(sapply(FullMD,gsub,pattern=" TRUE",replacement="TRUE"))
  ##perfect, now we just make sure the data classes match##
  sapply(FullExp,class)
  sapply(FullMD,class) 
  ##Ok, whew.. now for the tricky part: combining the FullMD and FullExp sets based on multiple matching columns##
  mergedData <- FullExp %>% left_join(FullMD, by=c("Date","Number","worm_present", "Sex"))
  FinalData = mergedData
  names(FinalData)
  FinalData = FinalData[,c(1:1930,1932,1934:1936,1938)]
  names(FinalData)
  ##create a new column for log(mass)
  sapply(FinalData, class)
  FinalData = transform(FinalData, weight = as.numeric(weight))
  FinalData$LogMass = log(FinalData$weight)
  FinalData = FinalData[,c(1:1933,1935:1936)]
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
  ##remove all rows with NAs-reduces our sample size down to 291##
  Final = FinalData[complete.cases(FinalData), ]
  ##write the Matrix for keeps##
  write.csv(Final, file = "PathMatrix.csv", row.names = FALSE)
  
  ##Ok move on to general data formatting in preparation for our loop##
  ##format PathMatrix##
  dat <- read.csv("PathMatrix.csv", header = T)
  dat2 <- dat[,-1]
  rownames(dat2) <- dat[,1]
  names(dat2)
  dat2 = dat2[,c(1930:1934, 1:1929)]
  names(dat2)
  ##now we want to iterate over every column in the matrix, creating submatrixes and running stats##
  output = file("SEM_Diversity_Code.R")
  cat("", file=output, append=FALSE)
  cat("library(sem)", "\n", "library(semPlot)", "\n","\n","\n",
      'dat <- read.csv("PathMatrix.csv", header = T)', "\n",
      'dat2 <- dat[,-1]', "\n",
      'rownames(dat2) <- dat[,1]', "\n",
      'dat2 = dat2[,c(1930:1934, 1:1929)]', "\n","\n","\n",
      file = output, append=TRUE)
  for (i in 6:1934) {
    x=colnames(dat2)[i] 
    output = file("SEM_Diversity_Code.R", open = "a")
    cat((paste0("dat", i, "= dat2[,c(1:5,",i,")]","\n",
                "S <-cov(dat", i, ")", "\n","N <- dim(dat", i, ")[1]", "\n",
                "RAM <- specifyModel()", "\n",
                "CrossDir -> worm_present, CI, NA", "\n",
                "CrossDir -> LogMass,CM, NA", "\n",
                "CrossDir -> Alpha_diversity_mean_2000,CD, NA", "\n",
                "CrossDir <-> CrossDir, C, NA", "\n",
                "CrossDir -> ", x, " ,CE, NA", "\n",
                "Sex -> ", x, " ,SE, NA", "\n",
                "Sex -> Alpha_diversity_mean_2000,SD, NA", "\n",
                "Sex -> LogMass,SM, NA", "\n",
                "Sex <-> Sex, S, NA", "\n",
                "LogMass -> ", x, " ,ME, NA", "\n",
                "LogMass -> Alpha_diversity_mean_2000,MD, NA", "\n",
                "LogMass <-> LogMass, M, NA", "\n",
                x, " <-> Alpha_diversity_mean_2000, ED, NA", "\n",
                x, " <-> worm_present, EI, NA", "\n",
                x, " <-> ", x, " , E, NA", "\n",
                "worm_present <-> Alpha_diversity_mean_2000, ID, NA", "\n",
                "worm_present <-> worm_present, I, NA", "\n",
                "Alpha_diversity_mean_2000 <-> Alpha_diversity_mean_2000,D, NA", "\n","\n","\n",
                "sem.out <- sem(RAM, S, N)", "\n", "x=summary(sem.out)", "\n",
                "y=standardizedCoefficients(sem.out)", "\n",
                'capture.output(x, y, file = ("semoutput_' , x , '.txt"))', '\n',
                "\n","\n","\n", "##Onto Next##", "\n")), file=output)
    close(output)
  }

  