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

#First load in the RPKM adjusted gene expression counts loop by analysis variable
rpkm.indir = file.path(anaDir,"00")
#Specify output directory
gne.ntwk.outdir = file.path(anaDir,"01")

#for each liver related pathway do:
lapply(1:length(liverPtwVect), function(i){
  ptw.name = liverPtwVect[i]
  rpkm.infilename = paste0(ptw.name,"_dge_rpkm.RDS")
  
  dge.rpkm = readRDS(file.path(rpkm.indir,rpkm.infilename))
  #Get gene symbols
  gene.names = colnames(dge.rpkm$counts)
  samids     = rownames(dge.rpkm$counts)
  
  gene.network = random_network(p=length(gene.names),n_modules = NULL)

  #plot(gene.network) create a separate plotting script for gene networks
  #create heatmaps for plotting the correlation structure of the simulated
  #gene expression matrix; Compare this with the correlation structure of the
  #actual gene expression matrix that was used as the reference
  
  ##Save the resulting gene networks
  gne.ntwk.infilename = paste0(ptw.name,"_network.RDS")
  saveRDS(gene.network,file.path(gne.ntwk.outdir,gne.ntwk.infilename))
  
  generated.gene.exp.data = gen_rnaseq(n=length(samids),
                                       gene.network,
                                       dge.rpkm$counts)
  #For the simulated expression data, give genes and subjects the same
  #name as those in the reference data
  colnames(generated.gene.exp.data[["x"]])=colnames(generated.gene.exp.data[["reference"]])
  rownames(generated.gene.exp.data[["x"]])=rownames(generated.gene.exp.data[["reference"]])
  
  #Save the simulated gene expression data in
  sim.gene.infilename = paste0(ptw.name,"_sim_gnexp.RDS")
  saveRDS(generated.gene.exp.data,file.path(gne.ntwk.outdir,sim.gene.infilename))
  
  #Resulting simulated gene expression data can be found 
  msg=paste0("Resulting simulated gene expression data can be found here: ",
             file.path(gne.ntwk.outdir,sim.gene.infilename))
  msg
})
