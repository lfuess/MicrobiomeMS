library(data.table)
library(tidyr)
##load in the reads##
genenames = read.csv("normalizedreads_TagSeq.csv", header=TRUE)
genenames=separate(data = genenames, col = Sample, into = c("gene", "postperiod"), sep = "\\.")
genenames=genenames[,c(1:2)]


##and parse each file!##

Caulo = read.csv('Corr_results_tau_Caulobacteraceae.csv', check.names = FALSE)
Caulo=setDT(Caulo, keep.rownames = TRUE)[]
dim(Caulo)
Caulo=Caulo[,c(1,593)]
##make the GOMWU matrix##
colnames(Caulo) <- c("gene", "tau")
##write out the tau values for supplemental file##
write.csv(Caulo, file="Caulobacteraceae_taus.csv", row.names = FALSE)
##finish up GOMWU input##
Merged_Caulo = merge(Caulo, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names; we can do this because there are no duplicated gene names (before the period)##
Merged_Caulo = transform(Merged_Caulo, truegenes=paste(gene, postperiod, sep="."))
##replace NA with 0##
Merged_Caulo[is.na(Merged_Caulo)] <- 0
dim(Merged_Caulo)
Merged_Caulo=Merged_Caulo[,c(4,2)]
##writeitout##
write.csv(Merged_Caulo, file="Caulobacteraceae_GOMWU.csv", row.names = FALSE)

##next one##
Methyl = read.csv('Corr_results_tau_Methylobacteriaceae.csv', check.names = FALSE)
Methyl=setDT(Methyl, keep.rownames = TRUE)[]
dim(Methyl)
Methyl=Methyl[,c(1,525)]
##make the GOMWU matrix##
colnames(Methyl) <- c("gene", "tau")
##write out the tau values for supplemental file##
write.csv(Methyl, file="Methylobacteriaceae_taus.csv", row.names = FALSE)
##finish up GOMWU input##
Merged_Methyl = merge(Methyl, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names; we can do this because there are no duplicated gene names (before the period)##
Merged_Methyl = transform(Merged_Methyl, truegenes=paste(gene, postperiod, sep="."))
##replace NA with 0##
Merged_Methyl[is.na(Merged_Methyl)] <- 0
dim(Merged_Methyl)
Merged_Methyl=Merged_Methyl[,c(4,2)]
##writeitout##
write.csv(Merged_Methyl, file="Methylobacteriaceae_GOMWU.csv", row.names = FALSE)

##next one##
Spart = read.csv('Corr_results_tau_Spartobacteria.csv', check.names = FALSE)
Spart=setDT(Spart, keep.rownames = TRUE)[]
dim(Spart)
Spart=Spart[,c(1,490)]
##make the GOMWU matrix##
colnames(Spart) <- c("gene", "tau")
##write out the tau values for supplemental file##
write.csv(Spart, file="Spartobacteria_taus.csv", row.names = FALSE)
##finish up GOMWU input##
Merged_Spart = merge(Spart, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names; we can do this because there are no duplicated gene names (before the period)##
Merged_Spart = transform(Merged_Spart, truegenes=paste(gene, postperiod, sep="."))
##replace NA with 0##
Merged_Spart[is.na(Merged_Spart)] <- 0
dim(Merged_Spart)
Merged_Spart=Merged_Spart[,c(4,2)]
##writeitout##
write.csv(Merged_Spart, file="Spartobacteria_GOMWU.csv", row.names = FALSE)

##next one##
Beta = read.csv('Corr_results_tau_Betaproteobacteria.csv', check.names = FALSE)
Beta=setDT(Beta, keep.rownames = TRUE)[]
dim(Beta)
Beta=Beta[,c(1,641)]
##make the GOMWU matrix##
colnames(Beta) <- c("gene", "tau")
##write out the tau values for supplemental file##
write.csv(Beta, file="Betaproteobacteria_taus.csv", row.names = FALSE)
##finish up GOMWU input##
Merged_Beta = merge(Beta, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names; we can do this because there are no duplicated gene names (before the period)##
Merged_Beta = transform(Merged_Beta, truegenes=paste(gene, postperiod, sep="."))
##replace NA with 0##
Merged_Beta[is.na(Merged_Beta)] <- 0
dim(Merged_Beta)
Merged_Beta=Merged_Beta[,c(4,2)]
##write it out##
write.csv(Merged_Beta, file="Betaproteobacteria_GOMWU.csv", row.names = FALSE)


##next one##
Halom = read.csv('Corr_results_tau_Halomonadaceae.csv', check.names = FALSE)
Halom=setDT(Halom, keep.rownames = TRUE)[]
dim(Halom)
Halom=Halom[,c(1,395)]
##make the GOMWU matrix##
colnames(Halom) <- c("gene", "tau")
##write out the tau values for supplemental file##
write.csv(Halom, file="Halomonadaceae_taus.csv", row.names = FALSE)
##finish up GOMWU input##
Merged_Halom = merge(Halom, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names; we can do this because there are no duplicated gene names (before the period)##
Merged_Halom = transform(Merged_Halom, truegenes=paste(gene, postperiod, sep="."))
##replace NA with 0##
Merged_Halom[is.na(Merged_Halom)] <- 0
dim(Merged_Halom)
Merged_Halom=Merged_Halom[,c(4,2)]
##writeitout##
write.csv(Merged_Halom, file="Halomonadaceae_GOMWU.csv", row.names = FALSE)


##next one##
Chitin = read.csv('Corr_results_tau_Chitinophagaceae.csv', check.names = FALSE)
Chitin=setDT(Chitin, keep.rownames = TRUE)[]
dim(Chitin)
Chitin=Chitin[,c(1,391)]
##make the GOMWU matrix##
colnames(Chitin) <- c("gene", "tau")
##write out the tau values for supplemental file##
write.csv(Chitin, file="Chitinophagaceae_taus.csv", row.names = FALSE)
##finish up GOMWU input##
Merged_Chitin = merge(Chitin, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names; we can do this because there are no duplicated gene names (before the period)##
Merged_Chitin = transform(Merged_Chitin, truegenes=paste(gene, postperiod, sep="."))
##replace NA with 0##
Merged_Chitin[is.na(Merged_Chitin)] <- 0
dim(Merged_Chitin)
Merged_Chitin=Merged_Chitin[,c(4,2)]
##writeitout##
write.csv(Merged_Chitin, file="Chitinophagaceae_GOMWU.csv", row.names = FALSE)


##next one##
Chitin = read.csv('Corr_results_tau_Chitinophagaceae.csv', check.names = FALSE)
Chitin=setDT(Chitin, keep.rownames = TRUE)[]
dim(Chitin)
Chitin=Chitin[,c(1,391)]
##make the GOMWU matrix##
colnames(Chitin) <- c("gene", "tau")
##write out the tau values for supplemental file##
write.csv(Chitin, file="Chitinophagaceae_taus.csv", row.names = FALSE)
##finish up GOMWU input##
Merged_Chitin = merge(Chitin, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names; we can do this because there are no duplicated gene names (before the period)##
Merged_Chitin = transform(Merged_Chitin, truegenes=paste(gene, postperiod, sep="."))
##replace NA with 0##
Merged_Chitin[is.na(Merged_Chitin)] <- 0
dim(Merged_Chitin)
Merged_Chitin=Merged_Chitin[,c(4,2)]
##writeitout##
write.csv(Merged_Chitin, file="Chitinophagaceae_GOMWU.csv", row.names = FALSE)


##next one##
Clost = read.csv('Corr_results_tau_Clostridiaceae.csv', check.names = FALSE)
Clost=setDT(Clost, keep.rownames = TRUE)[]
dim(Clost)
Clost=Clost[,c(1,928)]
##make the GOMWU matrix##
colnames(Clost) <- c("gene", "tau")
##write out the tau values for supplemental file##
write.csv(Clost, file="Clostridiaceae_taus.csv", row.names = FALSE)
##finish up GOMWU input##
Merged_Clost = merge(Clost, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names; we can do this because there are no duplicated gene names (before the period)##
Merged_Clost = transform(Merged_Clost, truegenes=paste(gene, postperiod, sep="."))
##replace NA with 0##
Merged_Clost[is.na(Merged_Clost)] <- 0
dim(Merged_Clost)
Merged_Clost=Merged_Clost[,c(4,2)]
##writeitout##
write.csv(Merged_Clost, file="Clostridiaceae_GOMWU.csv", row.names = FALSE)


##next one##
Orb = read.csv('Corr_results_tau_Orbaceae.csv', check.names = FALSE)
Orb=setDT(Orb, keep.rownames = TRUE)[]
dim(Orb)
Orb=Orb[,c(1,308)]
##make the GOMWU matrix##
colnames(Orb) <- c("gene", "tau")
##write out the tau values for supplemental file##
write.csv(Orb, file="Orbaceae_taus.csv", row.names = FALSE)
##finish up GOMWU input##
Merged_Orb = merge(Orb, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names; we can do this because there are no duplicated gene names (before the period)##
Merged_Orb = transform(Merged_Orb, truegenes=paste(gene, postperiod, sep="."))
##replace NA with 0##
Merged_Orb[is.na(Merged_Orb)] <- 0
dim(Merged_Orb)
Merged_Orb=Merged_Orb[,c(4,2)]
##writeitout##
write.csv(Merged_Orb, file="Orbaceae_GOMWU.csv", row.names = FALSE)


##next one##
Geod = read.csv('Corr_results_tau_Geodermatophilaceae.csv', check.names = FALSE)
Geod=setDT(Geod, keep.rownames = TRUE)[]
dim(Geod)
Geod=Geod[,c(1,663)]
##make the GOMWU matrix##
colnames(Geod) <- c("gene", "tau")
##write out the tau values for supplemental file##
write.csv(Geod, file="Geodermatophilaceae_taus.csv", row.names = FALSE)
##finish up GOMWU input##
Merged_Geod = merge(Geod, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names; we can do this because there are no duplicated gene names (before the period)##
Merged_Geod = transform(Merged_Geod, truegenes=paste(gene, postperiod, sep="."))
##replace NA with 0##
Merged_Geod[is.na(Merged_Geod)] <- 0
dim(Merged_Geod)
Merged_Geod=Merged_Geod[,c(4,2)]
##writeitout##
write.csv(Merged_Geod, file="Geodermatophilaceae_GOMWU.csv", row.names = FALSE)


##next one##
Nocard = read.csv('Corr_results_tau_Nocardiaceae.csv', check.names = FALSE)
Nocard=setDT(Nocard, keep.rownames = TRUE)[]
dim(Nocard)
Nocard=Nocard[,c(1,894)]
##make the GOMWU matrix##
colnames(Nocard) <- c("gene", "tau")
##write out the tau values for supplemental file##
write.csv(Nocard, file="Nocardiaceae_taus.csv", row.names = FALSE)
##finish up GOMWU input##
Merged_Nocard = merge(Nocard, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names; we can do this because there are no duplicated gene names (before the period)##
Merged_Nocard = transform(Merged_Nocard, truegenes=paste(gene, postperiod, sep="."))
##replace NA with 0##
Merged_Nocard[is.na(Merged_Nocard)] <- 0
dim(Merged_Nocard)
Merged_Nocard=Merged_Nocard[,c(4,2)]
##writeitout##
write.csv(Merged_Nocard, file="Nocardiaceae_GOMWU.csv", row.names = FALSE)


##next one##
Gp10 = read.csv('Corr_results_tau_Gp10.csv', check.names = FALSE)
Gp10=setDT(Gp10, keep.rownames = TRUE)[]
dim(Gp10)
Gp10=Gp10[,c(1,507)]
##make the GOMWU matrix##
colnames(Gp10) <- c("gene", "tau")
##write out the tau values for supplemental file##
write.csv(Gp10, file="Gp10_taus.csv", row.names = FALSE)
##finish up GOMWU input##
Merged_Gp10 = merge(Gp10, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names; we can do this because there are no duplicated gene names (before the period)##
Merged_Gp10 = transform(Merged_Gp10, truegenes=paste(gene, postperiod, sep="."))
##replace NA with 0##
Merged_Gp10[is.na(Merged_Gp10)] <- 0
dim(Merged_Gp10)
Merged_Gp10=Merged_Gp10[,c(4,2)]
##writeitout##
write.csv(Merged_Gp10, file="Gp10_GOMWU.csv", row.names = FALSE)


##next one##
Chlam = read.csv('Corr_results_tau_Chlamydiales.csv', check.names = FALSE)
Chlam=setDT(Chlam, keep.rownames = TRUE)[]
dim(Chlam)
Chlam=Chlam[,c(1,378)]
##make the GOMWU matrix##
colnames(Chlam) <- c("gene", "tau")
##write out the tau values for supplemental file##
write.csv(Chlam, file="Chlamydiales_taus.csv", row.names = FALSE)
##finish up GOMWU input##
Merged_Chlam = merge(Chlam, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names; we can do this because there are no duplicated gene names (before the period)##
Merged_Chlam = transform(Merged_Chlam, truegenes=paste(gene, postperiod, sep="."))
##replace NA with 0##
Merged_Chlam[is.na(Merged_Chlam)] <- 0
dim(Merged_Chlam)
Merged_Chlam=Merged_Chlam[,c(4,2)]
##writeitout##
write.csv(Merged_Chlam, file="Chlamydiales_GOMWU.csv", row.names = FALSE)


##next one##
Rubro = read.csv('Corr_results_tau_Rubrobacteraceae.csv', check.names = FALSE)
Rubro=setDT(Rubro, keep.rownames = TRUE)[]
dim(Rubro)
Rubro=Rubro[,c(1,353)]
##make the GOMWU matrix##
colnames(Rubro) <- c("gene", "tau")
##write out the tau values for supplemental file##
write.csv(Rubro, file="Rubrobacteraceae_taus.csv", row.names = FALSE)
##finish up GOMWU input##
Merged_Rubro = merge(Rubro, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names; we can do this because there are no duplicated gene names (before the period)##
Merged_Rubro = transform(Merged_Rubro, truegenes=paste(gene, postperiod, sep="."))
##replace NA with 0##
Merged_Rubro[is.na(Merged_Rubro)] <- 0
dim(Merged_Rubro)
Merged_Rubro=Merged_Rubro[,c(4,2)]
##writeitout##
write.csv(Merged_Rubro, file="Rubrobacteraceae_GOMWU.csv", row.names = FALSE)


##next one##
ISXI = read.csv('Corr_results_tau_Incertae_Sedis_XI.csv', check.names = FALSE)
ISXI=setDT(ISXI, keep.rownames = TRUE)[]
dim(ISXI)
ISXI=ISXI[,c(1,482)]
##make the GOMWU matrix##
colnames(ISXI) <- c("gene", "tau")
##write out the tau values for supplemental file##
write.csv(ISXI, file="Incertae_Sedis_XI_taus.csv", row.names = FALSE)
##finish up GOMWU input##
Merged_ISXI = merge(ISXI, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names; we can do this because there are no duplicated gene names (before the period)##
Merged_ISXI = transform(Merged_ISXI, truegenes=paste(gene, postperiod, sep="."))
##replace NA with 0##
Merged_ISXI[is.na(Merged_ISXI)] <- 0
dim(Merged_ISXI)
Merged_ISXI=Merged_ISXI[,c(4,2)]
##writeitout##
write.csv(Merged_ISXI, file="Incertae_Sedis_XI_GOMWU.csv", row.names = FALSE)


##next one##
Pepto = read.csv('Corr_results_tau_Peptostreptococcaceae.csv', check.names = FALSE)
Pepto=setDT(Pepto, keep.rownames = TRUE)[]
dim(Pepto)
Pepto=Pepto[,c(1,910)]
##make the GOMWU matrix##
colnames(Pepto) <- c("gene", "tau")
##write out the tau values for supplemental file##
write.csv(Pepto, file="Peptostreptococcaceae_taus.csv", row.names = FALSE)
##finish up GOMWU input##
Merged_Pepto = merge(Pepto, genenames, by = "gene", all.y = TRUE)
##add back in the post period part to gene names; we can do this because there are no duplicated gene names (before the period)##
Merged_Pepto = transform(Merged_Pepto, truegenes=paste(gene, postperiod, sep="."))
##replace NA with 0##
Merged_Pepto[is.na(Merged_Pepto)] <- 0
dim(Merged_Pepto)
Merged_Pepto=Merged_Pepto[,c(4,2)]
##writeitout##
write.csv(Merged_Pepto, file="Peptostreptococcaceae_GOMWU.csv", row.names = FALSE)




  