##This is a script to determine ISXIbial families which are siginificantly correlated to many genes##
library(data.table)
library(tidyverse)
##input our data frames##
Sub1 = read.csv("Corr_results_pval_Sub1.csv", check.names = FALSE)
##remove the repeptitive rows##
Sub1=Sub1[c(5001:5305)]
Sub1=Sub1[c(1:5000),]
##set non-significant values to NA##
Sub1[Sub1>0.05]=NA

Sub2 = read.csv("Corr_results_pval_Sub2.csv", check.names = FALSE)
##remove the repeptitive rows##
Sub2=Sub2[c(5001:5305)]
Sub2=Sub2[c(1:5000),]
##set non-significant values to NA##
Sub2[Sub2>0.05]=NA

Sub3 = read.csv("Corr_results_pval_Sub3.csv", check.names = FALSE)
##remove the repeptitive rows##
Sub3=Sub3[c(5001:5305)]
Sub3=Sub3[c(1:5000),]
##set non-significant values to NA##
Sub3[Sub3>0.05]=NA

Sub4 = read.csv("Corr_results_pval_Sub4.csv", check.names = FALSE)
##remove the repeptitive rows##
Sub4=Sub4[c(5001:5305)]
Sub4=Sub4[c(1:5000),]
##set non-significant values to NA##
Sub4[Sub4>0.05]=NA

Sub5 = read.csv("Corr_results_pval_Sub5.csv", check.names = FALSE)
##remove the repeptitive rows##
Sub5=Sub5[c(5001:5305)]
Sub5=Sub5[c(1:5000),]
##set non-significant values to NA##
Sub5[Sub5>0.05]=NA

Sub6 = read.csv("Corr_results_pval_Sub6.csv", check.names = FALSE)
Sub6[844:850,1:10]
##remove the repeptitive rows##
Sub6=Sub6[c(846:1150)]
Sub6=Sub6[c(1:845),]
##set non-significant values to NA##
Sub6[Sub6>0.05]=NA

##ok so now we bind them all together##
All=rbind(Sub1,Sub2,Sub3,Sub4,Sub5,Sub6)
sum(!is.na(All)) ##get total of all sig correlations for MS stats##


##and now we look for those families in the top 5%##
All=as.data.frame(t(All))
##start by replacing non-sig (n/as) with 1000000)
All[is.na(All)] <- 1000000
##create a column for row sums##
All$sums = rowSums(All)
##order by sums##
All_sorted <- All[order(All$sums),] 
All_sorted$sums
##pull out the top 5% of families (roughly 15)##
All_top_sorted = All_sorted[1:15,]
row.names(All_top_sorted)
is.na(All_top_sorted) <- All_top_sorted == 1000000
##generate a list of top families for reference##
All_familes = All_top_sorted[1]
write.csv(All_familes, file="top5families.csv")

##ok now we want just the top 5% of genes that are correlated to any given family##
All_genes=rbind(Sub1,Sub2,Sub3,Sub4,Sub5,Sub6)
##and follow code from above##
  ##start by replacing non-sig (n/as) with 1000000)
  All_genes[is.na(All_genes)] <- 1000000
  ##create a column for row sums##
  All_genes$sums = rowSums(All_genes)
  ##order by sums##
  All_sorted_genes <- All_genes[order(All_genes$sums),] 
  ##pull out the top 5% of genes (roughly 1292)##
  ##we actually do 1297 because there is a large tie for genes which are correlated to the same number of genes; this is the best break##
  ##code which validates this point is shown at the end of this script##
  All_top_sorted_genes = All_sorted_genes[1:1297,]
  is.na(All_top_sorted_genes) <- All_top_sorted_genes == 1000000
  
##now parse out lists of significant genes in top 5% of both family and genes##
Caulobacteraceae = All_top_sorted_genes[,c("Caulobacteraceae_Prop","sums")]
Caulobacteraceae = Caulobacteraceae[complete.cases(Caulobacteraceae), ]
Caulobacteraceae = setDT(Caulobacteraceae, keep.rownames = TRUE)[]
Caulobacteraceae = Caulobacteraceae[,c(1)]
write.csv(Caulobacteraceae, "Caulobacteraceae_Sig_Genes.csv", row.names = FALSE)

Methylobacteriaceae = All_top_sorted_genes[,c("Methylobacteriaceae_Prop","sums")]
Methylobacteriaceae = Methylobacteriaceae[complete.cases(Methylobacteriaceae), ]
Methylobacteriaceae = setDT(Methylobacteriaceae, keep.rownames = TRUE)[]
Methylobacteriaceae = Methylobacteriaceae[,c(1)]
write.csv(Methylobacteriaceae, "Methylobacteriaceae_Sig_Genes.csv", row.names = FALSE)

Spartobacteria = All_top_sorted_genes[,c("Spartobacteria_unclassified_Prop","sums")]
Spartobacteria = Spartobacteria[complete.cases(Spartobacteria), ]
Spartobacteria = setDT(Spartobacteria, keep.rownames = TRUE)[]
Spartobacteria = Spartobacteria[,c(1)]
write.csv(Spartobacteria, "Spartobacteria_Sig_Genes.csv", row.names = FALSE)

Betaproteobacteria = All_top_sorted_genes[,c("Betaproteobacteria_unclassified_Prop","sums")]
Betaproteobacteria = Betaproteobacteria[complete.cases(Betaproteobacteria), ]
Betaproteobacteria = setDT(Betaproteobacteria, keep.rownames = TRUE)[]
Betaproteobacteria = Betaproteobacteria[,c(1)]
write.csv(Betaproteobacteria, "Betaproteobacteria_Sig_Genes.csv", row.names = FALSE)

Halomonadaceae = All_top_sorted_genes[,c("Halomonadaceae_Prop","sums")]
Halomonadaceae = Halomonadaceae[complete.cases(Halomonadaceae), ]
Halomonadaceae = setDT(Halomonadaceae, keep.rownames = TRUE)[]
Halomonadaceae = Halomonadaceae[,c(1)]
write.csv(Halomonadaceae, "Halomonadaceae_Sig_Genes.csv", row.names = FALSE)

Chitinophagaceae = All_top_sorted_genes[,c("Chitinophagaceae_Prop","sums")]
Chitinophagaceae = Chitinophagaceae[complete.cases(Chitinophagaceae), ]
Chitinophagaceae = setDT(Chitinophagaceae, keep.rownames = TRUE)[]
Chitinophagaceae = Chitinophagaceae[,c(1)]
write.csv(Chitinophagaceae, "Chitinophagaceae_Sig_Genes.csv", row.names = FALSE)

Clostridiaceae = All_top_sorted_genes[,c("Clostridiaceae_1_Prop","sums")]
Clostridiaceae = Clostridiaceae[complete.cases(Clostridiaceae), ]
Clostridiaceae = setDT(Clostridiaceae, keep.rownames = TRUE)[]
Clostridiaceae = Clostridiaceae[,c(1)]
write.csv(Clostridiaceae, "Clostridiaceae_Sig_Genes.csv", row.names = FALSE)

Orbaceae = All_top_sorted_genes[,c("Orbaceae_Prop","sums")]
Orbaceae = Orbaceae[complete.cases(Orbaceae), ]
Orbaceae = setDT(Orbaceae, keep.rownames = TRUE)[]
Orbaceae = Orbaceae[,c(1)]
write.csv(Orbaceae, "Orbaceae_Sig_Genes.csv", row.names = FALSE)

Geodermatophilaceae = All_top_sorted_genes[,c("Geodermatophilaceae_Prop","sums")]
Geodermatophilaceae = Geodermatophilaceae[complete.cases(Geodermatophilaceae), ]
Geodermatophilaceae = setDT(Geodermatophilaceae, keep.rownames = TRUE)[]
Geodermatophilaceae = Geodermatophilaceae[,c(1)]
write.csv(Geodermatophilaceae, "Geodermatophilaceae_Sig_Genes.csv", row.names = FALSE)

Nocardiaceae = All_top_sorted_genes[,c("Nocardiaceae_Prop","sums")]
Nocardiaceae = Nocardiaceae[complete.cases(Nocardiaceae), ]
Nocardiaceae = setDT(Nocardiaceae, keep.rownames = TRUE)[]
Nocardiaceae = Nocardiaceae[,c(1)]
write.csv(Nocardiaceae, "Nocardiaceae_Sig_Genes.csv", row.names = FALSE)

Gp10 = All_top_sorted_genes[,c("Gp10_unclassified_Prop","sums")]
Gp10 = Gp10[complete.cases(Gp10), ]
Gp10 = setDT(Gp10, keep.rownames = TRUE)[]
Gp10 = Gp10[,c(1)]
write.csv(Gp10, "Gp10_Sig_Genes.csv", row.names = FALSE)

Chlamydiales = All_top_sorted_genes[,c("Chlamydiales_unclassified_Prop","sums")]
Chlamydiales = Chlamydiales[complete.cases(Chlamydiales), ]
Chlamydiales = setDT(Chlamydiales, keep.rownames = TRUE)[]
Chlamydiales = Chlamydiales[,c(1)]
write.csv(Chlamydiales, "Chlamydiales_Sig_Genes.csv", row.names = FALSE)

Rubrobacteraceae = All_top_sorted_genes[,c("Rubrobacteraceae_Prop","sums")]
Rubrobacteraceae = Rubrobacteraceae[complete.cases(Rubrobacteraceae), ]
Rubrobacteraceae = setDT(Rubrobacteraceae, keep.rownames = TRUE)[]
Rubrobacteraceae = Rubrobacteraceae[,c(1)]
write.csv(Rubrobacteraceae, "Rubrobacteraceae_Sig_Genes.csv", row.names = FALSE)

Incertae_Sedis_XI = All_top_sorted_genes[,c("Incertae_Sedis_XI_Prop","sums")]
Incertae_Sedis_XI = Incertae_Sedis_XI[complete.cases(Incertae_Sedis_XI), ]
Incertae_Sedis_XI = setDT(Incertae_Sedis_XI, keep.rownames = TRUE)[]
Incertae_Sedis_XI = Incertae_Sedis_XI[,c(1)]
write.csv(Incertae_Sedis_XI, "Incertae_Sedis_XI_Sig_Genes.csv", row.names = FALSE)

Peptostreptococcaceae = All_top_sorted_genes[,c("Peptostreptococcaceae_Prop","sums")]
Peptostreptococcaceae = Peptostreptococcaceae[complete.cases(Peptostreptococcaceae), ]
Peptostreptococcaceae = setDT(Peptostreptococcaceae, keep.rownames = TRUE)[]
Peptostreptococcaceae = Peptostreptococcaceae[,c(1)]
write.csv(Peptostreptococcaceae, "Peptostreptococcaceae_Sig_Genes.csv", row.names = FALSE)

##good.. we can turn these into corrleation matrixes for generating tau values##
##read in expression data##
reads = read.csv("normalizedreads_TagSeq.csv", check.names = FALSE)
reads$Sample <- gsub("\\..*","", reads$Sample)

##read in Prop data##
Prop = read.csv("FamilyPropData.csv", check.names = FALSE)

BetaMatrix = merge(Betaproteobacteria, reads,
                   by.x = "rn", by.y = "Sample", all.x = TRUE,
                   sort = TRUE,  no.dups = FALSE,)
BetaMatrix = BetaMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
BetaMatrix = as.data.frame(t(BetaMatrix))
BetaMatrix= setDT(BetaMatrix, keep.rownames = TRUE)[]
Beta = Prop[c("SampleRNA","Betaproteobacteria_unclassified_Prop")]
BetaMatrix = merge(BetaMatrix, Beta,
                   by.x = "rn", by.y = "SampleRNA", all.x = TRUE,
                   sort = TRUE, no.dups = FALSE)
dim(BetaMatrix)
BetaMatrix = BetaMatrix[,c(2:641)]
write.csv(BetaMatrix, file = "Betaproteobacteria_Tau_Matrix.csv", row.names = FALSE)


CauloMatrix = merge(Caulobacteraceae, reads,
                    by.x = "rn", by.y = "Sample", all.x = TRUE,
                    sort = TRUE,  no.dups = FALSE,)
CauloMatrix = CauloMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
CauloMatrix = as.data.frame(t(CauloMatrix))
CauloMatrix= setDT(CauloMatrix, keep.rownames = TRUE)[]
Caulo = Prop[c("SampleRNA","Caulobacteraceae_Prop")]
CauloMatrix = merge(CauloMatrix, Caulo,
                    by.x = "rn", by.y = "SampleRNA", all.x = TRUE,
                    sort = TRUE, no.dups = FALSE)
dim(CauloMatrix)
CauloMatrix = CauloMatrix[,c(2:593)]
write.csv(CauloMatrix, file = "Caulobacteraceae_Tau_Matrix.csv", row.names = FALSE)


ChitinMatrix = merge(Chitinophagaceae, reads,
                    by.x = "rn", by.y = "Sample", all.x = TRUE,
                    sort = TRUE,  no.dups = FALSE,)
ChitinMatrix = ChitinMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
ChitinMatrix = as.data.frame(t(ChitinMatrix))
ChitinMatrix= setDT(ChitinMatrix, keep.rownames = TRUE)[]
Chitin = Prop[c("SampleRNA","Chitinophagaceae_Prop")]
ChitinMatrix = merge(ChitinMatrix, Chitin,
                    by.x = "rn", by.y = "SampleRNA", all.x = TRUE,
                    sort = TRUE, no.dups = FALSE)
dim(ChitinMatrix)
ChitinMatrix = ChitinMatrix[,c(2:391)]
write.csv(ChitinMatrix, file = "Chitinophagaceae_Tau_Matrix.csv", row.names = FALSE)


ClostMatrix = merge(Clostridiaceae, reads,
                    by.x = "rn", by.y = "Sample", all.x = TRUE,
                    sort = TRUE,  no.dups = FALSE,)
ClostMatrix = ClostMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
ClostMatrix = as.data.frame(t(ClostMatrix))
ClostMatrix= setDT(ClostMatrix, keep.rownames = TRUE)[]
Clost = Prop[c("SampleRNA","Clostridiaceae_1_Prop")]
ClostMatrix = merge(ClostMatrix, Clost,
                    by.x = "rn", by.y = "SampleRNA", all.x = TRUE,
                    sort = TRUE, no.dups = FALSE)
dim(ClostMatrix)
ClostMatrix = ClostMatrix[,c(2:928)]
write.csv(ClostMatrix, file = "Clostridiales_Tau_Matrix.csv", row.names = FALSE)


GeodMatrix = merge(Geodermatophilaceae, reads,
                   by.x = "rn", by.y = "Sample", all.x = TRUE,
                   sort = TRUE,  no.dups = FALSE,)
GeodMatrix = GeodMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
GeodMatrix = as.data.frame(t(GeodMatrix))
GeodMatrix= setDT(GeodMatrix, keep.rownames = TRUE)[]
Geod = Prop[c("SampleRNA","Geodermatophilaceae_Prop")]
GeodMatrix = merge(GeodMatrix, Geod,
                   by.x = "rn", by.y = "SampleRNA", all.x = TRUE,
                   sort = TRUE, no.dups = FALSE)
dim(GeodMatrix)
GeodMatrix = GeodMatrix[,c(2:663)]
write.csv(GeodMatrix, file = "Geodermatophilaceae_Tau_Matrix.csv", row.names = FALSE)



Gp10Matrix = merge(Gp10, reads,
                   by.x = "rn", by.y = "Sample", all.x = TRUE,
                   sort = TRUE,  no.dups = FALSE,)
Gp10Matrix = Gp10Matrix %>% remove_rownames %>% column_to_rownames(var="rn")
Gp10Matrix = as.data.frame(t(Gp10Matrix))
Gp10Matrix= setDT(Gp10Matrix, keep.rownames = TRUE)[]
Gp10 = Prop[c("SampleRNA","Gp10_unclassified_Prop")]
Gp10Matrix = merge(Gp10Matrix, Gp10,
                   by.x = "rn", by.y = "SampleRNA", all.x = TRUE,
                   sort = TRUE, no.dups = FALSE)
dim(Gp10Matrix)
Gp10Matrix = Gp10Matrix[,c(2:507)]
write.csv(Gp10Matrix, file = "Gp10_Tau_Matrix.csv", row.names = FALSE)



HaloMatrix = merge(Halomonadaceae, reads,
                   by.x = "rn", by.y = "Sample", all.x = TRUE,
                   sort = TRUE,  no.dups = FALSE,)
HaloMatrix = HaloMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
HaloMatrix = as.data.frame(t(HaloMatrix))
HaloMatrix= setDT(HaloMatrix, keep.rownames = TRUE)[]
Halo = Prop[c("SampleRNA","Halomonadaceae_Prop")]
HaloMatrix = merge(HaloMatrix, Halo,
                   by.x = "rn", by.y = "SampleRNA", all.x = TRUE,
                   sort = TRUE, no.dups = FALSE)
dim(HaloMatrix)
HaloMatrix = HaloMatrix[,c(2:395)]
write.csv(HaloMatrix, file = "Halomonadaceae_Tau_Matrix.csv", row.names = FALSE)


ISXIMatrix = merge(Incertae_Sedis_XI, reads,
                    by.x = "rn", by.y = "Sample", all.x = TRUE,
                    sort = TRUE,  no.dups = FALSE,)
ISXIMatrix = ISXIMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
ISXIMatrix = as.data.frame(t(ISXIMatrix))
ISXIMatrix= setDT(ISXIMatrix, keep.rownames = TRUE)[]
ISXI = Prop[c("SampleRNA","Incertae_Sedis_XI_Prop")]
ISXIMatrix = merge(ISXIMatrix, ISXI,
                    by.x = "rn", by.y = "SampleRNA", all.x = TRUE,
                    sort = TRUE, no.dups = FALSE)
dim(ISXIMatrix)
ISXIMatrix = ISXIMatrix[,c(2:482)]
write.csv(ISXIMatrix, file = "Incertae_Sedis_XI_Tau_Matrix.csv", row.names = FALSE)


MethylMatrix = merge(Methylobacteriaceae, reads,
                     by.x = "rn", by.y = "Sample", all.x = TRUE,
                     sort = TRUE,  no.dups = FALSE,)
MethylMatrix = MethylMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
MethylMatrix = as.data.frame(t(MethylMatrix))
MethylMatrix= setDT(MethylMatrix, keep.rownames = TRUE)[]
Methyl = Prop[c("SampleRNA","Methylobacteriaceae_Prop")]
MethylMatrix = merge(MethylMatrix, Methyl,
                     by.x = "rn", by.y = "SampleRNA", all.x = TRUE,
                     sort = TRUE, no.dups = FALSE)
dim(MethylMatrix)
MethylMatrix = MethylMatrix[,c(2:525)]
write.csv(MethylMatrix, file = "Methylobacteriaceae_Tau_Matrix.csv", row.names = FALSE)


NocardMatrix = merge(Nocardiaceae, reads,
                     by.x = "rn", by.y = "Sample", all.x = TRUE,
                     sort = TRUE,  no.dups = FALSE,)
NocardMatrix = NocardMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
NocardMatrix = as.data.frame(t(NocardMatrix))
NocardMatrix= setDT(NocardMatrix, keep.rownames = TRUE)[]
Nocard = Prop[c("SampleRNA","Nocardiaceae_Prop")]
NocardMatrix = merge(NocardMatrix, Nocard,
                     by.x = "rn", by.y = "SampleRNA", all.x = TRUE,
                     sort = TRUE, no.dups = FALSE)
dim(NocardMatrix)
NocardMatrix = NocardMatrix[,c(2:894)]
write.csv(NocardMatrix, file = "Nocardiaceae_Tau_Matrix.csv", row.names = FALSE)


OrbMatrix = merge(Orbaceae, reads,
                     by.x = "rn", by.y = "Sample", all.x = TRUE,
                     sort = TRUE,  no.dups = FALSE,)
OrbMatrix = OrbMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
OrbMatrix = as.data.frame(t(OrbMatrix))
OrbMatrix= setDT(OrbMatrix, keep.rownames = TRUE)[]
Orb = Prop[c("SampleRNA","Orbaceae_Prop")]
OrbMatrix = merge(OrbMatrix, Orb,
                     by.x = "rn", by.y = "SampleRNA", all.x = TRUE,
                     sort = TRUE, no.dups = FALSE)
dim(OrbMatrix)
OrbMatrix = OrbMatrix[,c(2:308)]
write.csv(OrbMatrix, file = "Orbaceae_Tau_Matrix.csv", row.names = FALSE)


PeptoMatrix = merge(Peptostreptococcaceae, reads,
                  by.x = "rn", by.y = "Sample", all.x = TRUE,
                  sort = TRUE,  no.dups = FALSE,)
PeptoMatrix = PeptoMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
PeptoMatrix = as.data.frame(t(PeptoMatrix))
PeptoMatrix= setDT(PeptoMatrix, keep.rownames = TRUE)[]
Pepto = Prop[c("SampleRNA","Peptostreptococcaceae_Prop")]
PeptoMatrix = merge(PeptoMatrix, Pepto,
                  by.x = "rn", by.y = "SampleRNA", all.x = TRUE,
                  sort = TRUE, no.dups = FALSE)
dim(PeptoMatrix)
PeptoMatrix = PeptoMatrix[,c(2:910)]
write.csv(PeptoMatrix, file = "Peptostreptococcaceae_Tau_Matrix.csv", row.names = FALSE)


RubroMatrix = merge(Rubrobacteraceae, reads,
                    by.x = "rn", by.y = "Sample", all.x = TRUE,
                    sort = TRUE,  no.dups = FALSE,)
RubroMatrix = RubroMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
RubroMatrix = as.data.frame(t(RubroMatrix))
RubroMatrix= setDT(RubroMatrix, keep.rownames = TRUE)[]
Rubro = Prop[c("SampleRNA","Rubrobacteraceae_Prop")]
RubroMatrix = merge(RubroMatrix, Rubro,
                    by.x = "rn", by.y = "SampleRNA", all.x = TRUE,
                    sort = TRUE, no.dups = FALSE)
dim(RubroMatrix)
RubroMatrix = RubroMatrix[,c(2:353)]
write.csv(RubroMatrix, file = "Rubrobacteraceae_Tau_Matrix.csv", row.names = FALSE)


SpartoMatrix = merge(Spartobacteria, reads,
                     by.x = "rn", by.y = "Sample", all.x = TRUE,
                     sort = TRUE,  no.dups = FALSE,)
SpartoMatrix = SpartoMatrix %>% remove_rownames %>% column_to_rownames(var="rn")
SpartoMatrix = as.data.frame(t(SpartoMatrix))
SpartoMatrix= setDT(SpartoMatrix, keep.rownames = TRUE)[]
Sparto = Prop[c("SampleRNA","Spartobacteria_unclassified_Prop")]
SpartoMatrix = merge(SpartoMatrix, Sparto,
                     by.x = "rn", by.y = "SampleRNA", all.x = TRUE,
                     sort = TRUE, no.dups = FALSE)
dim(SpartoMatrix)
SpartoMatrix = SpartoMatrix[,c(2:490)]
write.csv(SpartoMatrix, file = "Spartobacteria_Tau_Matrix.csv", row.names = FALSE)



##code to validate how we chose our top ~5% of genes##
##use this in place of existing code from lines 68-79##

##and replace sigs with 1s##
All_genes[All_genes<0.05]=1
##start by replacing non-sig (n/as) with 0)
All_genes[is.na(All_genes)] <- 0
##create a column for row sums##
All_genes$sums = rowSums(All_genes)
##order by sums##
All_sorted_genes <- All_genes[order(All_genes$sums, decreasing = TRUE),] 
##pull out the top 5% of genes (roughly 1292)##
All_top_sorted_genes = All_sorted_genes[1:1290,]
All_test = All_top_sorted_genes[,c(1,305)]
write.csv(All_test, file="ISXIbetest.csv")
is.na(All_top_sorted_genes) <- All_top_sorted_genes == 1000000
