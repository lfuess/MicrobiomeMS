##This Script takes a list of tau values from significant genes and creates GOMWU input##
##first we load our package##
library(tidyr)

##Next we load in our data##
TauRes = read.csv("Diversity_Corr_results_Tau.csv", check.names = FALSE, row.names = NULL)
dim(TauRes)
##remove the repeptitive columns; we only need first and last##
Taus=TauRes[c(1,1931)]
colnames(Taus) <- c("gene", "tau")

##Now let's write that bit out to a CSV##
write.csv(Taus, file="Diversity_SigGenes_TauValues.csv", row.names = FALSE)

##and finally let's make a GOMWU matrix!##
##import list of all genes##
genenames = read.csv("normalizedreads_TagSeq.csv", check.names = FALSE)
##fix the names for merging, but keep the part after the . for later; the pre-period portion is unique, so this works fine##
genenames=separate(data = genenames, col = Sample, into = c("gene", "postperiod"), sep = "\\.")
genenames=genenames[c(1:2)]
##and merge; make sure to keep all genes though cause you need them there for GOMWU##
Merged = merge(Taus, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names##
Merged = transform(Merged, truegenes=paste(gene, postperiod, sep="."))
##We set all nonsignificant genes to 0 (so replace NAs with 0s)##
Merged[is.na(Merged)] <- 0
dim(Merged)
names(Merged)
##grab just the full gene names and the stat##
Merged=Merged[c(4,2)]

##Now we just save it for GOMWU analyses##
write.csv(Merged, file="Diversity_Taus_GOMWU.csv", row.names = FALSE)
