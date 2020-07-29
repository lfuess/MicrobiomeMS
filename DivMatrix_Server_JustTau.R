##this script will take your matrix of genes you want tau values for and run correlations, spitting out tau into a new matrix##

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

MB_Corr = read.csv("Diversity_Tau_Matrix.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
t=cor.test.tau(MB_Corr)
warnings()
write.table(t, file ="Diversity_Corr_results_Tau.csv",sep=",")
