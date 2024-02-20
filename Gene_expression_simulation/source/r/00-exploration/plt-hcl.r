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

## set directories
in.dir.dist = paste(res.dir,'dist',sep='/');

## update specimen flag to inclue a variable across all
## types if there is more than 1 specimen type
if (length(spcFlags)>1) {
  spcFlags = c('all',spcFlags)
  spcLabls = c('All Specimen Types',spcLabls)
}

## loop through specimen types
for(s in 1:length(spcFlags)) {
  spcFlag = spcFlags[s];
  
  ## plotting parameters
  par(mfrow=c(1,2),oma=c(0,0,4,1));
  
  ## subset metadata to specimen type
  if (spcFlag=='all') {
    mta.spc = mta
  } else {
    mta.spc = mta[mta$spct==spcFlag,]
    mta.spc = mta.spc[!mta.spc$samid %in% outliers,]
  }
  ## loop through standardization types
  for(d in 1:length(stdFlags)) {
    stdFlag = stdFlags[d];
    stdLabl = stdLabls[d];
    
    ## import data
    filename = paste(in.dir.dist,'/',spcFlag,'_alltp_euc_dist_tmm_normalized_',stdFlag,'.RData',sep='');
    load(file=filename);
    
    ## construct HCL object and print as dendrogram
    hcl.res=hclust(dist.res,method=cluster.method);
    hcl.ids = attr(dist.res, "labels")[hcl.res$order]
    hcl.dnd=as.dendrogram(hcl.res);
    plot(colDendo(hcl.dnd,mta.spc$spctc[match(attr(dist.res, "labels"),mta.spc$samid)][hcl.res$order],hcl.ids),axes=T,horiz=T,xlab='Euclidean Distance',xlim=c(max(hcl.res$height),max(hcl.res$height)*-.1));
    if(d==1) {
      legend(max(hcl.res$height)*.75,1.12*nrow(mta.spc),legend=unique(mta.spc$spctl),fill=unique(mta.spc$spctc),cex=0.8,border = FALSE,xpd=TRUE,horiz=T);
    }
    mtext(stdLabl,cex=0.8);	
    
    ## Mark outliers
    outs.idx = which(hcl.ids %in% unlist(outliers))
    plot.edge = par()$usr[2]-(abs(par()$usr[1])+abs(par()$usr[2]))/par()$fin[1]*(par()$mai[4]*1.3)
    points(rep(plot.edge,length(outs.idx)),outs.idx,pch=1.5,col='blue',xpd=T)
  }
}

#For Kernel Density Estimates-> perform a 2-sample Kolmogorov-Smirnov Test