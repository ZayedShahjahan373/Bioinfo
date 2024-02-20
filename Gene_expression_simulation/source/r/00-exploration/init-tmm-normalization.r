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
suppressMessages(source(file.path("C:/Users/zayed/OneDrive/Desktop/Research_v2.0/SeqNet_playthru","source","r","init-SeqNet-sim-main.r")))

#Load TMM-normalization factors from the init-eda script
#Load the simulated genes 
tmm.dir = file.path(anaDir,"00")
gnexp.mat.indir = file.path(anaDir,"01")

annot.df = read.table(file.path(anaDir,"annotation_processed.tab.gz"),header = T)
#list.files(tmm.dir)

lapply(1:length(liverPtwVect), function(i){
  ptw.name = liverPtwVect[i]
  

  tmm.factor.infile = paste0(ptw.name,"_tmm_factor.tab.gz")
  sim.gnexp.infile  = paste0(ptw.name,"_sim_gnexp.RDS")
  
  tmm.factor.df = read.table(file.path(tmm.dir,tmm.factor.infile))
  sim.gnexp.list = readRDS(file.path(gnexp.mat.indir,sim.gnexp.infile))
  
  gnex.mat = sim.gnexp.list[["x"]]
  
  #Subset annotation dataframe to only those in the pathway
  ptw.annot.df = annot.df %>% filter(gene_name %in% colnames(gnex.mat))
  
  
  lcpm.mat = filterUnwantedGenes(cpm(gnex.mat),clean.genes = annot.df$)
  
})