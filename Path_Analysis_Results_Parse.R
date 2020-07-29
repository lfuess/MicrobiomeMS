##Now we're going to parse our results to figure out what % of genes are still sig correlated to diversity##
library(dplyr)

##we're going to use loops and lists to do this##
##we'll only pull the stats for ExpressionxDiversity covarrience##
myFiles <- list.files(pattern="semoutput_.*txt")
mynames<-gsub(".*semoutput_ *(.*?) *.txt.*", "\\1", myFiles)
for (i in 1:length(myFiles)) assign(mynames[i], read.table(myFiles[i], skip = 27, nrows = 1))
df_list <- mget(ls()[sapply(ls(), function(x) is.data.frame(get(x)))])
Out=bind_rows(df_list, .id = "Gene")
colnames(Out) <- c("Gene", "Factor", "Estimate", "StdError", "Zval", "Pval")

##select just the columns we really want for later manipulation##
Out_select = Out[,c(1:3,6)]

##and then select the rows with significant values##
Sig = Out_select[Out_select$Pval < 0.05, ]

dim(Sig) ##and get a number##

##1014/1929 are still significant for Diversity##
##459/639 are still sig for Betaproteobacteria##
##484/591 are still sig for Caulobacteriaceae##
##182/389 are still sig for Chitophagaceae
##210/376 are still sig for Chlamydiales##
##662/926 are still sig for Clostridiaceae_1 ##
##427/661 are still sig for Geodermatophilaceae ##
##217/505 are still sig for GP10##
##324/393 are still sig for Halomonadaceae##
##316/480 are still sig for Incertae_Sedis_XI##
##410/523 are still sig for Methylobacteriaceae##
##658/923 are still sig for Nocardiaceae##
##174/306 are still sig for Orbaceae##
##594/908 are still sig for Peptostreptococcaceae##
##187/351 are still sig for Rubrobacteraceae##
##432/488 are still sig for Spartobacteria##
##write it out so that we can redo gene ontology analyses with the subset and compare##
#write.csv(Sig, "Diversity_Sig_PostPath.csv", row.names = FALSE)


