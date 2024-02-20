#############################################################################################################
# 
# This script performs exploratory data analysis on the reference data obtained
# from the kaggle link: https://www.kaggle.com/datasets/lachmann12/human-liver-rnaseq-gene-expression-903-samples/data
#
# Program:  init-SeqNet-sim-main.r
# Version:  RSEQREP 2.3.0
# Author:   zshahjahan
# Purpose:  Initializes important global variabels and functions. This includes knitR setup.
# Input:    
# Output:  	N/A
#############################################################################################################
suppressMessages(source(file.path("C:/Users/zayed/OneDrive/Desktop/Research_v2.0/SeqNet_playthru",
                                  "source","r","init-SeqNet-sim-main.r")))
gnexp.mat.indir = file.path(anaDir,"01")
#Run PCA for normalized vs non-normalized and standardized vs. non-standardized
#data

##perform PCA and store the results
lapply(1:length(liverPtwVect), function(i) {
  ptw.name = liverPtwVect[i]
  
  gnexp.infile.name = paste0(ptw.name,"_sim_gnexp.RDS")
  
  gnexp.mat.list = readRDS(file.path(gnexp.mat.indir,gnexp.infile.name))
  
  sim.gene.exp = gnexp.mat.list[["x"]]
  
  #Loop through TMM-normalization flags
  
  pca.res = prcomp(sim.gene.exp)
  
  
  ## store only the first three components (coefficents and rotated values)
  pca.res$rotation = pca.res$rotation[,1:3];
  pca.res$x = pca.res$x[,1:3];
  
  ## calculate percent explained variance
  pca.res$pvar=round(pca.res$sdev^2/sum(pca.res$sdev^2)*100,1);
  
})