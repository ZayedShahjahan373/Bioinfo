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
source('../r/init-analysis.r')

## set directories
in.dir.pca  = paste(res.dir,'pca',sep='/');
in.dir.dist = paste(res.dir,'dist',sep='/');

## update specimen flag to inclue a variable across all
## types if there is more than 1 specimen type
if (length(spcFlags)>1) {
  spcFlags = c('all',spcFlags)
  spcLabls = c('All Specimen Types',spcLabls)
}

## set variables
visLabl = '';
alpha = 0.05;

## loop through specimen types
for(s in 1:length(spcFlags)) {
  spcFlag = spcFlags[s];
  spcLabl = spcLabls[s];
  
  ## subset metadata to specimen type
  if (spcFlag=='all') {
    mta.spc = mta
  } else {
    mta.spc = mta[mta$spct==spcFlag,]
  }
  
  ## plotting parameters
  par(mfrow=c(3,1)); 
  par(mar=c(3.5,3.5,3.5,0.2),oma=c(0,0,1.5,0));
  
  ## PCA, Euc Dist MDS, and Spearman MDS
  plotPcaByCell(mta,in.dir.pca,spcFlag,spcLabl,'alltp',visLabl);
  plotMdsByCell(mta,in.dir.dist,type='euc',spcFlag,spcLabl,'alltp',visLabl,text='Euclidean Distance Standardized Variables');
  plotMdsByCell(mta,in.dir.dist,type='spm',spcFlag,spcLabl,'alltp',visLabl,text='1-Spearman Cor. Distance Standardized Variables');
  
  ## add legend
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend('top',legend=unique(mta.spc$spctl),fill=unique(mta.spc$spctc),horiz=T,cex=1.1,border = FALSE,xpd=TRUE);
}