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

lapply(1:length(liverPtwVect), function(i){
  ptw.name = liverPtwVect[i]
  
  gnexp.infile.name = paste0(ptw.name,"_sim_gnexp.RDS")
  
  gnexp.mat.list = readRDS(file.path(gnexp.mat.indir,gnexp.infile.name))
  
  corr.plot.list=lapply(1:length(refSim), function(j){
    mat.type = names(refSim)[j]
    gene.exp.mat = gnexp.mat.list[[mat.type]]
    corr.plt = get_correl_heatmap.func(gene.exp.mat,
                                       method = "pearson",
                                       sig.dig = 3)
    return(corr.plt)
  })
  #name the elements of the list 
  
  # Arrange the plots horizontally with specific legend adjustments
  plot_list_arranged = plot_grid(plotlist = corr.plot.list,
                                  align = 'h',
                                  ncol = length(corr.plot.list),
                                  rel_widths = c(1, 1),
                                 labels = c("reference","simulation"),
                                 label_size = 15)
  #Save the plots in the results directory
  plt.outfile.name = paste0(ptw.name,"_corr_htmp.png")
  
  ggsave(plt.outfile.name,
         plot_list_arranged,
         device ="png",
         path = pltDir,bg = "white",
         units = "in",
         width = 30,
         height = 20,
         dpi = "retina")
})