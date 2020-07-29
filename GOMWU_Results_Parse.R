##this script parases output from GOMWU analyses and pulls out a list of only those terms which are significant##

Bac = read.table("MWU_BP_new_Spartobacteria_GOMWU.csv", sep = " ", header = TRUE)
Bac = Bac[order(Bac$p.adj),]
Bac = Bac[Bac[,7]<0.05,]
write.csv(Bac, "Spartobacteria_Sig_GOs.csv")
