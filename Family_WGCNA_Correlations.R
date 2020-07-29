##This is a Script to generate WGCNA Correlation Matrixes##

library(WGCNA)
library(wesanderson)
library(extrafont)
## First we load up our data from the previous WGCNA script
 load("Stickeblack-all-networkConstruction.RData")
  #make sure to check it#
  table(moduleColors)
  
## Next load the trait data ##
  traitData = read.csv("FamilyPropData.csv");
  #again make sure to check it
  dim(traitData)
  names(traitData)
  ##pick the rows you need##
  traitData=traitData[c("SampleRNA","worm_present", "Bacillales_Incertae_Sedis_XII_Prop","Betaproteobacteria_unclassified_Prop",
                        "Bradyrhizobiaceae_Prop","Caulobacteraceae_Prop", "Chlamydiales_unclassified_Prop",
                        "Clostridiales_Incertae_Sedis_XI_Prop","Firmicutes_unclassified_Prop", "Geodermatophilaceae_Prop",
                        "Gp16_unclassified_Prop", "Halomonadaceae_Prop", "Methylobacteriaceae_Prop", "Microbacteriaceae_Prop",
                        "Moraxellaceae_Prop", "Nocardiaceae_Prop", "Spartobacteria_unclassified_Prop")]


##Next we combine the sets##
  datExpr=datExpr0
  Samples = rownames(datExpr);
  traitRows = match(Samples, traitData$SampleRNA);
  datTraits = traitData[traitRows, 1:17];
  rownames(datTraits) = traitData[traitRows, 1];
  names(datTraits)
  datTraits = datTraits[c(2:17)]
  #check
  datTraits
  collectGarbage();

##Make the Matrix! This is verbatum from WGCNA tutorials##

  # Define numbers of genes and samples
  nGenes = ncol(datExpr0);
  nSamples = nrow(datExpr0);
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  MEsnames=MEs
  moduleTraitCor = bicor(MEs, datTraits, use = "pairwise.complete.obs");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  # Will display correlations and their p-values
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  #par(mar = c(6, 8.5, 3, 3));
  # Make the tombstone plot
  colors = c("#3B9AB2","#EBCC2A", "#F21A00")
  colorRampPalette(colors)
  names(MEsnames) <- substring(names(MEs), 3)
  labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEsnames),
               ySymbols = names(MEs),
               xLabelsAngle = 60,
               colorLabels = FALSE,
               colors =  colorRampPalette(colors)(100),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               yColorWidth = 100,
               cex.text = 0.5,
               cex.lab = 0.6,
               yColorOffset = 0.01,
               zlim = c(-1,1))
               #plotLegend = FALSE)

##Now we generate gene significant for infection (worm_present)
  #again taken straight from tutorials
  worm_present = as.data.frame(datTraits$worm_present);
  names(worm_present) = "worm_present"
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  geneTraitSignificance = as.data.frame(cor(datExpr, worm_present, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(worm_present), sep="");
  names(GSPvalue) = paste("p.GS.", names(worm_present), sep="");


##Finally we export all of this info (gene significant, module membership, etc) to a csv##
  # again taken straight from WGCNA tutorial
  locus = names(datExpr)
  geneInfo0 = data.frame(locus = locus, moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
  modOrder = order(-abs(cor(MEs, worm_present, use = "p")));
  for (mod in 1:ncol(geneModuleMembership))
  {
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
  }

  geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.worm_present));
  geneInfo = geneInfo0[geneOrder, ]
  write.csv(geneInfo, file = "WGCNA_Family_Corrs.csv", row.names = FALSE)

  