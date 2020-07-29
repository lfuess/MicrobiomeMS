##This is a file to make a heatmap of all significant correlations between genes and families of interest##
library(data.table)

##load in each family's data and merge it to create a master matrix##

Beta = read.csv('Corr_results_tau_Betaproteobacteria.csv', check.names = FALSE)
Beta=setDT(Beta, keep.rownames = TRUE)[]
Beta = Beta[,c(1,641)]

Caulo = read.csv('Corr_results_tau_Caulobacteraceae.csv', check.names = FALSE)
Caulo=setDT(Caulo, keep.rownames = TRUE)[]
Caulo=Caulo[,c(1,593)]

Chito = read.csv('Corr_results_tau_Chitinophagaceae.csv', check.names = FALSE)
Chito=setDT(Chito, keep.rownames = TRUE)[]
Chito=Chito[,c(1,391)]

Chlam = read.csv('Corr_results_tau_Chlamydiales.csv', check.names = FALSE)
Chlam=setDT(Chlam, keep.rownames = TRUE)[]
Chlam=Chlam[,c(1,378)]

Clost = read.csv('Corr_results_tau_Clostridiaceae.csv', check.names = FALSE)
Clost=setDT(Clost, keep.rownames = TRUE)[]
Clost=Clost[,c(1,928)]

Geod = read.csv('Corr_results_tau_Geodermatophilaceae.csv', check.names = FALSE)
Geod=setDT(Geod, keep.rownames = TRUE)[]
Geod=Geod[,c(1,663)]

Gp10 = read.csv('Corr_results_tau_Gp10.csv', check.names = FALSE)
Gp10=setDT(Gp10, keep.rownames = TRUE)[]
Gp10=Gp10[,c(1,507)]

Halom = read.csv('Corr_results_tau_Halomonadaceae.csv', check.names = FALSE)
Halom=setDT(Halom, keep.rownames = TRUE)[]
Halom=Halom[,c(1,395)]

Incert = read.csv('Corr_results_tau_Incertae_Sedis_XI.csv', check.names = FALSE)
Incert=setDT(Incert, keep.rownames = TRUE)[]
Incert=Incert[,c(1,482)]

Methyl = read.csv('Corr_results_tau_Methylobacteriaceae.csv', check.names = FALSE)
Methyl=setDT(Methyl, keep.rownames = TRUE)[]
Methyl=Methyl[,c(1,525)]

Nocard = read.csv('Corr_results_tau_Nocardiaceae.csv', check.names = FALSE)
Nocard=setDT(Nocard, keep.rownames = TRUE)[]
Nocard=Nocard[,c(1,894)]

Orb = read.csv('Corr_results_tau_Orbaceae.csv', check.names = FALSE)
Orb=setDT(Orb, keep.rownames = TRUE)[]
Orb=Orb[,c(1,308)]

Pepto = read.csv('Corr_results_tau_Peptostreptococcaceae.csv', check.names = FALSE)
Pepto=setDT(Pepto, keep.rownames = TRUE)[]
Pepto=Pepto[,c(1,910)]

Rubro = read.csv('Corr_results_tau_Rubrobacteraceae.csv', check.names = FALSE)
Rubro=setDT(Rubro, keep.rownames = TRUE)[]
Rubro=Rubro[,c(1,353)]

Spart = read.csv('Corr_results_tau_Spartobacteria.csv', check.names = FALSE)
Spart=setDT(Spart, keep.rownames = TRUE)[]
Spart=Spart[,c(1,490)]

##merge all the things##
Master = merge(Beta, Caulo, by = "rn", all.x = TRUE, all.y = TRUE)
Master = merge(Master, Chito, by = "rn", all.x = TRUE, all.y = TRUE)
Master = merge(Master, Chlam, by = "rn", all.x = TRUE, all.y = TRUE)
Master = merge(Master, Clost, by = "rn", all.x = TRUE, all.y = TRUE)
Master = merge(Master, Geod, by = "rn", all.x = TRUE, all.y = TRUE)
Master = merge(Master, Gp10, by = "rn", all.x = TRUE, all.y = TRUE)
Master = merge(Master, Halom, by = "rn", all.x = TRUE, all.y = TRUE)
Master = merge(Master, Incert, by = "rn", all.x = TRUE, all.y = TRUE)
Master = merge(Master, Methyl, by = "rn", all.x = TRUE, all.y = TRUE)
Master = merge(Master, Nocard, by = "rn", all.x = TRUE, all.y = TRUE)
Master = merge(Master, Orb, by = "rn", all.x = TRUE, all.y = TRUE)
Master = merge(Master, Pepto, by = "rn", all.x = TRUE, all.y = TRUE)
Master = merge(Master, Rubro, by = "rn", all.x = TRUE, all.y = TRUE)
Master = merge(Master, Spart, by = "rn", all.x = TRUE, all.y = TRUE)

##remove rows you don't need##
Master = Master[c(6:1268)]

##rename columns##
names(Master) <- c("gene","Betaproteobacteria", "Caulobacteraceae", "Chitinophagaceae", "Chlamydiales", "Clostridiaceae 1",
                   "Geodermatophilaceae", "Gp10", "Halomonadaceae", "Incertae Sedis XI", "Methylobacteriaceae",
                   "Nocardiaceae", "Orbaceae", "Peptostreptococcaceae", "Rubrobacteraceae", "Spartobacteria")

##replace NAs with 0s##
Master[is.na(Master)] <- 0

##write it to a csv##
write.csv(Master, file = "TauHeatmap.csv", row.names = FALSE)

##and make the heatmap##

library(pheatmap)
library(RColorBrewer)
library(wesanderson)

data=read.csv("TauHeatmap.csv", row.names= 1)

pheatmap(data, show_rownames=F, color = colorRampPalette(c("#3B9AB2","grey","#F21A00"))(100))

