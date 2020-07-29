##define correlation function##
cor.test.p <- function(x){
  FUN <- function(x, y) cor.test(x, y, method="kendall", use="pairwise")[["p.value"]]
  z <- outer(
    colnames(x), 
    colnames(x), 
    Vectorize(function(i,j) FUN(x[,i], x[,j]))
  )
  dimnames(z) <- list(colnames(x), colnames(x))
  z
}

##First Subset (start with smallest)##
MB_Corr = read.csv("DiversityCorrMatrix_Sub1.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
p=cor.test.p(MB_Corr)
warnings()
write.table(p, file ="Corr_results_pval_Div_Sub1.csv",sep=",")

##Next Subset##
MB_Corr = read.csv("DiversityCorrMatrix_Sub2.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
p=cor.test.p(MB_Corr)
warnings()
write.table(p, file ="Corr_results_pval_Div_Sub2.csv",sep=",")

##next subset##
MB_Corr = read.csv("DiversityCorrMatrix_Sub3.csv", check.names = FALSE)
##and run it! (will take a while.. that's a lot of pairwise correlations to make!##
p=cor.test.p(MB_Corr)
warnings()
write.table(p, file ="Corr_results_pval_Div_Sub3.csv",sep=",")

