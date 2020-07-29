##this script takes correlation data output and parses it to just recieve the significant results of interest##
##first we load our package##
library(tidyr)

##Next we load in our data, one fram at a time; we split our data into subsets so the server wouldn't choke##
CorrRes = read.csv("Corr_results_pval_Div_Sub1.csv", check.names = FALSE, row.names = NULL)
dim(CorrRes)
##remove the repeptitive columns; we only need first and last##
CorrRes_test=CorrRes[c(1,10002)]
Sub1=CorrRes_test
##set non-significant values to NA#
Sub1[2][Sub1[2]>0.05]=NA
##Remove all rows with NAs, leaving only genes with significant corrs##
Sub1_done= Sub1[complete.cases(Sub1), ]


##now we just repeat this process for the next 2 subsets##
CorrRes = read.csv("Corr_results_pval_Div_Sub2.csv", check.names = FALSE, row.names = NULL)
dim(CorrRes)
##remove the repeptitive columns; we only need first and last##
CorrRes_test=CorrRes[c(1,10002)]
Sub2=CorrRes_test
##set non-significant values to NA#
Sub2[2][Sub2[2]>0.05]=NA
##Remove all rows with NAs, leaving only genes with significant corrs##
Sub2_done= Sub2[complete.cases(Sub2), ]


CorrRes = read.csv("Corr_results_pval_Div_Sub3.csv", check.names = FALSE, row.names = NULL)
dim(CorrRes)
##remove the repeptitive columns; we only need first and last##
CorrRes_test=CorrRes[c(1,5847)]
Sub3=CorrRes_test
##set non-significant values to NA#
Sub3[2][Sub3[2]>0.05]=NA
##Remove all rows with NAs, leaving only genes with significant corrs##
Sub3_done= Sub3[complete.cases(Sub3), ]

##now let's merge them all into 1 data file##
All=rbind(Sub1_done,Sub2_done,Sub3_done)
##Perfect now let's write it to a CSV##
write.csv(All, file="SigCorrResults_withPs.csv", row.names = FALSE)
