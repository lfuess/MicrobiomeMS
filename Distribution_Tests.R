##This is a script to analyze the distribution of p-values in our data set##

  ##input our data frames##
  Sub1 = read.csv("Corr_results_pval_Sub1.csv", check.names = FALSE)
  ##remove the repeptitive rows##
  Sub1=Sub1[c(5001:5304)]
  Sub1=Sub1[c(1:5000),]
  
  Sub2 = read.csv("Corr_results_pval_Sub2.csv", check.names = FALSE)
  ##remove the repeptitive rows##
  Sub2=Sub2[c(5001:5304)]
  Sub2=Sub2[c(1:5000),]
  
  Sub3 = read.csv("Corr_results_pval_Sub3.csv", check.names = FALSE)
  ##remove the repeptitive rows##
  Sub3=Sub3[c(5001:5304)]
  Sub3=Sub3[c(1:5000),]
  
  Sub4 = read.csv("Corr_results_pval_Sub4.csv", check.names = FALSE)
  ##remove the repeptitive rows##
  Sub4=Sub4[c(5001:5304)]
  Sub4=Sub4[c(1:5000),]
  
  Sub5 = read.csv("Corr_results_pval_Sub5.csv", check.names = FALSE)
  ##remove the repeptitive rows##
  Sub5=Sub5[c(5001:5304)]
  Sub5=Sub5[c(1:5000),]
  
  Sub6 = read.csv("Corr_results_pval_Sub6.csv", check.names = FALSE)
  ##remove the repeptitive rows##
  Sub6=Sub6[c(846:1149)]
  Sub6=Sub6[c(1:845),]

##bind all the data together into a single data frame##
All=rbind(Sub1,Sub2,Sub3,Sub4,Sub5,Sub6)

##graph the distribution##
hist( All[All > 0 & All < 1] )

##It definitely doesn't look uniform, but to be sure, we're going to use a Kolmogorov-Smirnov test##

Values=unlist(All)

ks.test(Values, "punif", 0, 1)

##yeah, definitely not uniformly distributed!  Yay!

