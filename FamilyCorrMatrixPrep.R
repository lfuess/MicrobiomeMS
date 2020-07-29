library(dplyr)
library(data.table)
library(reshape2)
##begin by grabbing info needed to turn OTUs into Familys (which make more biological sense)##
  ##start by reading in the file##
  OTUdata = read.delim("Microbe_taxonomy.txt", header = TRUE, sep = " ", row.names = NULL, stringsAsFactors=FALSE)
  ##double check it read in correctly##
  OTUdata
  ##pull out just the OTU and the family##
  OTUdata_real=OTUdata[c(1,6)]
  ##check##
  names(OTUdata_real)
  ##rename##
  names(OTUdata_real) <- c("OTU", "Family")
  OTUdata=OTUdata_real
##now back up and manipulate the massive metadata file##
  MetaData = read.csv("Metadata_w_OTUcounts.csv", header = TRUE, stringsAsFactors=FALSE)
  ##there's a whole lot of junk in there!!! let's whittle it down some##
  dim(MetaData)
  names(MetaData)  
  MetaData$X
  MetaData_real=MetaData[c(1,14:11584)]
  ##fiddling with rownames so we have some reference points/column names when we transpose##
  rownames(MetaData_real) <- MetaData_real[,1]
  ##ok and transpose it so we can join##
  MetaDatat=as.data.frame(t(MetaData_real))
  MetaDatat <- tibble::rownames_to_column(MetaDatat, "Info")
  MetaDatat <- MetaDatat[2:11572,]
##now we need to combine these two data sets.. a simple merge will do!##
  Combined=merge(MetaDatat, OTUdata,
        by.x = "Info", by.y = "OTU", all.x = TRUE,
        sort = TRUE,  no.dups = FALSE,)
  Combined=Combined[c(1,695,2:694)]  
##next we work on condensing down columns to gain one value per Family!##
  ##first sort by Family##
  Combined_sorted=Combined[order(Combined$Family),]
  dim(Combined_sorted)
  ##remove that OTU colum so it won't interfere with sums##
  Combined_sorted=Combined_sorted[c(2:695)]
  Combined_sorted_num=Combined_sorted
  ##convert everything to numeric##
  cols = c(2:694);    
  Combined_sorted_num[,cols] = apply(Combined_sorted_num[,cols], 2, function(x) as.numeric(as.character(x)));
  ##next lets sum up our duplicates##
  Summed = aggregate(. ~ Family, data=Combined_sorted_num, FUN=sum)
  names(Summed)
  ##now create sample totals
  Summed_NewR = rbind(Summed, data.frame(Family="Total",t(colSums(Summed[,-1]))))
  ##transpose it back...##
  res <- recast(Summed_NewR, variable~Family, value.var='value')
  colnames(res)
  res=res[c(1:297,299:307,298)]
  ##and finally loop over the data to create new columns for porportions (each Family count divided by the Total)##
  dim(res)
  names(res)
  setDT(res)[, paste0(names(res)[2:306],"_Prop") := lapply(.SD, `/`, res$Total), .SDcols = Acetobacteraceae:Xanthomonadaceae]
  ##magic!! now select out just the names and the proportion rows
  res=as.data.frame(res)
  names(res)
  res_final=res[c(1,308:612)]
  names(res_final)
##Woohoo you got it!!! Now we have to tackle the mess of matching to our RNAseq data!##
  ##rename that first column##
  colnames(res_final)[1]="Sample"
  ##merge in data that you're going to need to help match##
  names(MetaData)
  ##select what you'll need##
  MatchInfo=MetaData[c(1,2,4,10,11)]  
  ##and merge the two##
  FamInfo=merge(res_final, MatchInfo,
                 by.x = "Sample", by.y = "X", all.x = TRUE,
                 sort = TRUE,  no.dups = FALSE,)
  ##reorder##
  names(FamInfo)
  FamInfo=FamInfo[c(1,307:310,2:306)]
##awesome, that's done.. Now we need to go back and work with our RNAseq Data to create a file that is able to be merged##
  ##start with our normalized read output from the DESeq Analaysis with the TagSeq Paper (this has our final sample set)##
  Exp=read.csv("normalizedreads_TagSeq.csv", header = FALSE)
  Exp = Exp[1,]    
  Exp= as.data.frame(t(Exp))  
  ##now that we have our list of samples, we can read in the metadata##
  MetaDataExp=read.csv("Masterdata_TagSeq_SDH_23July2018.csv", header=TRUE, stringsAsFactors = FALSE)
  names(MetaDataExp)  
  MetaDataExp=MetaDataExp[c(2,20,29:30,31)]  
  colnames(MetaDataExp)[2]="Date"
  ##now that it's polished, let's merge stuff together##
  FullExp=merge(Exp, MetaDataExp,
               by.x = "1", by.y = "sample_ID", all.x = FALSE,
               sort = TRUE,  no.dups = FALSE,)
  ##rename the column##
  colnames(FullExp)[1]="Sample"
  ##fix the dates so they're consistent across data frames (FullMD has them listed as /201X)##
  FullExp <- as.data.frame(sapply(FullExp,gsub,pattern="/15",replacement="/2015"))
  FullExp <- as.data.frame(sapply(FullExp,gsub,pattern="9/14",replacement="9/2014"))
  FullExp <- as.data.frame(sapply(FullExp,gsub,pattern="0/14",replacement="0/2014"))
  FullExp <- as.data.frame(sapply(FullExp,gsub,pattern="7/14",replacement="7/2014"))
  FullExp <- as.data.frame(sapply(FullExp,gsub,pattern="3/14",replacement="3/2014"))
  FullExp <- as.data.frame(sapply(FullExp,gsub,pattern="4/14",replacement="4/2014"))
  FullExp <- as.data.frame(sapply(FullExp,gsub,pattern="2/14",replacement="2/2014"))
  levels(FullExp$Date)
  ##perfect, now we just make sure the data classes match##
  FullExp$Number=as.integer(as.character(FullExp$Number))
  FullExp$worm_present=as.logical(as.character(FullExp$worm_present))
  FullExp$Date=as.character(as.character(FullExp$Date))
  FullExp$Sex=as.character(as.character(FullExp$Sex))
  sapply(FullExp,class)
  sapply(FamInfo,class)
  
  ##Ok, whew.. now for the tricky part: combining the FullMD and FullExp sets based on multiple matching columns##
  mergedData <- FullExp %>% right_join(FamInfo, by=c("Date","Number","worm_present", "Sex"))
  ##and finally remove all the rows with NA sample names
  FinalData = mergedData[!is.na(mergedData$Sample.x),] 
  ##Order by sample, so it will work in WGCNA##
  FinalData=FinalData[order(FinalData$Sample.x),]
  ##rename some columns##
  colnames(FinalData)[1]="SampleRNA"
  colnames(FinalData)[6]="SampleMicrobiome"
  names(FinalData)
  ##pull out what you want for WGNCA##
  FinalData=FinalData[c(1,5,7:311)]
  write.table(FinalData, file ="FamilyPropData.csv",row.names=FALSE, sep=",")
##and finally make a correlation matrix##
  ExpFull=read.csv("normalizedreads_TagSeq.csv", header = FALSE)
  ExpFull=as.data.frame(t(ExpFull))
  names(ExpFull) <- lapply(ExpFull[1, ], as.character)
  ExpFull <- ExpFull[-1,] 
  FullMatrix=merge(FinalData, ExpFull,
                by.x = "SampleRNA", by.y = "Sample", all.x = TRUE,
                sort = TRUE,  no.dups = FALSE,)
##divide into 5 subsets and write the files##
  ##divide them up##
  FullMatrix_Sub1 = FullMatrix[,c(308:5307,3:307) ]
  FullMatrix_Sub2 = FullMatrix[,c(5308:10307,3:307) ]
  FullMatrix_Sub3 = FullMatrix[,c(10308:15307,3:307) ]
  FullMatrix_Sub4 = FullMatrix[,c(15308:20307,3:307) ]
  FullMatrix_Sub5 = FullMatrix[,c(20308:25307,3:307) ]
  FullMatrix_Sub6 = FullMatrix[,c(25308:26152,3:307) ]
  ##fix column names##
  names(FullMatrix_Sub1) <- gsub("\\..*","", names(FullMatrix_Sub1))
  names(FullMatrix_Sub2) <- gsub("\\..*","", names(FullMatrix_Sub2))
  names(FullMatrix_Sub3) <- gsub("\\..*","", names(FullMatrix_Sub3))
  names(FullMatrix_Sub4) <- gsub("\\..*","", names(FullMatrix_Sub4))
  names(FullMatrix_Sub5) <- gsub("\\..*","", names(FullMatrix_Sub5))
  names(FullMatrix_Sub6) <- gsub("\\..*","", names(FullMatrix_Sub6))
  ##and write them out##
  write.csv(FullMatrix_Sub1, file = "FamilyCorrMatrix_Sub1.csv", row.names = FALSE, quote = FALSE)
  write.csv(FullMatrix_Sub2, file = "FamilyCorrMatrix_Sub2.csv", row.names = FALSE, quote = FALSE)
  write.csv(FullMatrix_Sub3, file = "FamilyCorrMatrix_Sub3.csv", row.names = FALSE, quote = FALSE)
  write.csv(FullMatrix_Sub4, file = "FamilyCorrMatrix_Sub4.csv", row.names = FALSE, quote = FALSE)
  write.csv(FullMatrix_Sub5, file = "FamilyCorrMatrix_Sub5.csv", row.names = FALSE, quote = FALSE)
  write.csv(FullMatrix_Sub6, file = "FamilyCorrMatrix_Sub6.csv", row.names = FALSE, quote = FALSE)
  

