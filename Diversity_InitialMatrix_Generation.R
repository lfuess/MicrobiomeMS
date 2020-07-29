##Generating a large matrix for both correlation and WGCNA analyses##

##load in packages##
library(dplyr)

##first load all the data and select what you need##
T2 = read.csv("alpha_2000.csv", header = TRUE)
FullMD = T2[c(3,4,6,12:13,18)]
names(FullMD) = c("Sample","Date", "Number", "Sex", "worm_present", "Alpha_diversity_mean_2000")
##ok now we can merge it with our RNAseq stuff##
##To get the samples we're working with, start with our normalized read output from the DESeq Analaysis with the TagSeq Paper (this has our final sample set)##
Exp=read.csv("normalizedreads_TagSeq.csv", header = FALSE)
Exp = Exp[1,]    
Exp= as.data.frame(t(Exp))  
##now that we have our list of samples, we can read in the metadata from TagSeq, which will use to merge the datasets and create universal sample names##
MetaDataExp=read.csv("Masterdata_TagSeq_SDH_23July2018.csv", header=TRUE)
names(MetaDataExp)  
MetaDataExp=MetaDataExp[c(2,20,29:30,31)]  
colnames(MetaDataExp)[2]="Date"
##now that it's polished, let's merge stuff together##
FullExp=merge(Exp, MetaDataExp,
              by.x = "1", by.y = "sample_ID", all.x = FALSE,
              sort = TRUE,  no.dups = FALSE,)
##rename the column##
colnames(FullExp)[1]="Sample"
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
##and finally remove all the rows with NA sample names 
FinalData = mergedData[!is.na(mergedData$Sample.x),] 

##Order by sample##
FinalData=FinalData[order(FinalData$Sample.x),]
##rename some columns##
colnames(FinalData)[1]="SampleRNA"
colnames(FinalData)[6]="SampleMicrobiome"
names(FinalData)
##pull out what you want for WGNCA##
FinalData=FinalData[c(1,5,7)]
write.table(FinalData, file ="DiversityData2000.csv",row.names=FALSE, sep=",") 

##ok the last thing we want to do is make our correlation matrix.##
##we need to merge our normalized reads with our Diversity figs##
  ##make a diversity fig matrix##
  Diversity = FinalData[c(1,3)]
  ##load in your readcounts##
  Reads = read.csv("normalizedreads_TagSeq.csv", check.names = FALSE, header = FALSE)
  Reads = as.data.frame(t(Reads))
  colnames(Reads) <- as.character(unlist(Reads[1,]))
  Reads = Reads[-1, ]
  ##merge them together##
  FullMatrix=merge(Reads, Diversity,
                by.x = "Sample", by.y = "SampleRNA", all.x = FALSE,
                sort = TRUE,  no.dups = FALSE,)
  ##divide them up##
  FullMatrix_Sub1 = FullMatrix[,c(2:10001,25847) ]
  FullMatrix_Sub2 = FullMatrix[,c(10002:20001,25847) ]
  FullMatrix_Sub3 = FullMatrix[,c(20002:25847) ]
  ##fix column names##
  names(FullMatrix_Sub1) <- gsub("\\.", "", names(FullMatrix_Sub1))
  names(FullMatrix_Sub2) <- gsub("\\.", "", names(FullMatrix_Sub2))
  names(FullMatrix_Sub3) <- gsub("\\.", "", names(FullMatrix_Sub3))
##and write them out##
write.csv(FullMatrix_Sub1, file = "DiversityCorrMatrix_Sub1.csv", row.names = FALSE, quote = FALSE)
write.csv(FullMatrix_Sub2, file = "DiversityCorrMatrix_Sub2.csv", row.names = FALSE, quote = FALSE)
write.csv(FullMatrix_Sub3, file = "DiversityCorrMatrix_Sub3.csv", row.names = FALSE, quote = FALSE)

  