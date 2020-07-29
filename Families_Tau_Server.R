
##define correlation function##
cor.test.tau <- function(x){
  FUN <- function(x, y) cor.test(x, y, method="kendall", use="pairwise")[["estimate"]]
  z <- outer(
    colnames(x),
    colnames(x),
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

Caulo_Corr = read.csv("Caulobacteraceae_Tau_Matrix.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
t_Caulo=cor.test.tau(Caulo_Corr)
warnings()
write.table(t_Caulo, file ="Corr_results_tau_Caulobacteraceae.csv",sep=",")

Methyl_Corr = read.csv("Methylobacteriaceae_Tau_Matrix.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
t_Methyl=cor.test.tau(Methyl_Corr)
warnings()
write.table(t_Methyl, file ="Corr_results_tau_Methylobacteriaceae.csv",sep=",")

Sparto_Corr = read.csv("Spartobacteria_Tau_Matrix.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
t_Sparto=cor.test.tau(Sparto_Corr)
warnings()
write.table(t_Sparto, file ="Corr_results_tau_Spartobacteria.csv",sep=",")

Beta_Corr = read.csv("Betaproteobacteria_Tau_Matrix.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
t_Beta=cor.test.tau(Beta_Corr)
warnings()
write.table(t_Beta, file ="Corr_results_tau_Betaproteobacteria.csv",sep=",")

Halo_Corr = read.csv("Halomonadaceae_Tau_Matrix.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
t_Halo=cor.test.tau(Halo_Corr)
warnings()
write.table(t_Halo, file ="Corr_results_tau_Halomonadaceae.csv",sep=",")

Chitin_Corr = read.csv("Chitinophagaceae_Tau_Matrix.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
t_Chitin=cor.test.tau(Chitin_Corr)
warnings()
write.table(t_Chitin, file ="Corr_results_tau_Chitinophagaceae.csv",sep=",")

Clost_Corr = read.csv("Clostridiaceae_Tau_Matrix.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
t_Clost=cor.test.tau(Clost_Corr)
warnings()
write.table(t_Clost, file ="Corr_results_tau_Clostridiaceae.csv",sep=",")

Orb_Corr = read.csv("Orbaceae_Tau_Matrix.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
t_Orb=cor.test.tau(Orb_Corr)
warnings()
write.table(t_Orb, file ="Corr_results_tau_Orbaceae.csv",sep=",")

Geod_Corr = read.csv("Geodermatophilaceae_Tau_Matrix.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
t_Geod=cor.test.tau(Geod_Corr)
warnings()
write.table(t_Geod, file ="Corr_results_tau_Geodermatophilaceae.csv",sep=",")

Nocard_Corr = read.csv("Nocardiaceae_Tau_Matrix.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
t_Nocard=cor.test.tau(Nocard_Corr)
warnings()
write.table(t_Nocard, file ="Corr_results_tau_Nocardiaceae.csv",sep=",")

Gp10_Corr = read.csv("Gp10_Tau_Matrix.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
t_Gp10=cor.test.tau(Gp10_Corr)
warnings()
write.table(t_Gp10, file ="Corr_results_tau_Gp10.csv",sep=",")

Chlam_Corr = read.csv("Chlamydiales_Tau_Matrix.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
t_Chlam=cor.test.tau(Chlam_Corr)
warnings()
write.table(t_Chlam, file ="Corr_results_tau_Chlamydiales.csv",sep=",")

Rubro_Corr = read.csv("Rubrobacteraceae_Tau_Matrix.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
t_Rubro=cor.test.tau(Rubro_Corr)
warnings()
write.table(t_Rubro, file ="Corr_results_tau_Rubrobacteraceae.csv",sep=",")

ISXI_Corr = read.csv("Incertae_Sedis_XI_Tau_Matrix.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
t_ISXI=cor.test.tau(ISXI_Corr)
warnings()
write.table(t_ISXI, file ="Corr_results_tau_Incertae_Sedis_XI.csv",sep=",")

Pepto_Corr = read.csv("Peptostreptococcaceae_Tau_Matrix.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
t_Pepto=cor.test.tau(Pepto_Corr)
warnings()
write.table(t_Pepto, file ="Corr_results_tau_Peptostreptococcaceae.csv",sep=",")
