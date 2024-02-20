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
#retrieve the first element of the simulated gene expression dataset 
gnex.dir = file.path(anaDir,"01")

gnex.df.list = readRDS(file.path(gnex.dir,"bile_sim_gnexp.RDS"))
sim.gene.exp = gnex.df.list[["x"]]
#retrieve the metadata
mta.df = read.csv(file.path(anaDir,"project_metadata.csv"))

#retrieve the immunology dataset
imm.df = read.csv(file.path(anaDir,"immunology_outcome.csv"))

enet.mod = cv.glmnet(x=sim.gene.exp,y=imm.df[,"immune_response"])

enet.mod = glmnet(x=sim.gene.exp,y=imm.df[,"immune_response"])

plot(enet.mod)
print(enet.mod)
coef.set=coef(enet.mod,s=0.1)

cv.enet.mod = cv.glmnet(x=sim.gene.exp,y=imm.df[,"immune_response"],nfolds = 25)

plot(cv.enet.mod)

nBoot=500
#Create a subdirectory under analysis and module-03 and store the bootstrapped
#results
boot.results.dir = file.path(anaDir,"03","bootstrapped_results")

if(!dir.exists(boot.results.dir)){
  dir.create(boot.results.dir)
}


#I want to check the stability of the elastic net coefficient selectio process
#For this I will need to bootstrap the results
for (b in 1:nBoot) {
  #From the original dataset, resample a certain number of points
  #Resample 25% of the original dataset
  bootstrapped.samid = sample(imm.df$SAMID,replace=T)
  y.boot = imm.df %>% filter(SAMID %in% bootstrapped.samid)
  #Keep only the immunogenictiy variable
  y.boot = data.frame(lg10_immune=log10(y.boot[,"immune_response"]),
                      SAMID = y.boot$SAMID)
  #Do the same thing for the metadata and the design matrix
  x.boot = sim.gene.exp[rownames(sim.gene.exp) %in% bootstrapped.samid,]
  
  x.small = as.data.frame(x.boot[,c(3:9)])
  x.small$SAMID= rownames(x.small)
  mydat = merge(y.boot,x.small,by="SAMID")
  
  mydat.mta = merge(mydat,mta.df,by="SAMID",all.x=T)
  
  p=ggplot(mydat.mta,
           aes(x=ACOT8,y=lg10_immune))+
    geom_point(aes(shape=tranplType),#Shape aesthetics for points go here
               stat="identity",size=4,alpha=0.5) #No statistical transformations
  p+facet_grid(cols=vars(sexVar))+theme_bw()  
  
  base.plot=ggplot(mydat.mta,aes(x=-log2(AKR1C3),y=lg10_immune))
  
  #Melt utility
  
  
  connected_scatter=base.plot+geom_line()+geom_point()
  connected_scatter
  
  hist_immu = ggplot(mydat.mta, aes(x=lg10_immune))+
    geom_histogram()+labs("Immunogenicity Assay")+theme(text = element_text(family="berling"))
  #geom_line or geom_path for maps, geom_polygon for frequency polygons
  mydat.mta = mydat.mta %>% mutate(Day=gsub("T","",tmpt))
  #Take the log2 of the read counts before passing to the cv.glmnet function
  x.boot = log10((x.boot+0.05))
  cv.boot = cv.glmnet(data.frame(x.boot),y.boot$lg10_immune,nfolds = nFold)
  
  mydat.mta.site = mydat.mta %>% filter(treatSite=="Site A")
  mydat.mta.site %>% group_by(PATID) %>% 
    summarise(Avg_lg10 = mean(lg10_immune))%>%
  ggplot(aes(x=PATID,y=Avg_lg10,group=PATID,color=treatType))+
    geom_point()+
    geom_path()+
    scale_color_viridis(discrete=T)
  
  b.plot = ggplot(mydat.mta,aes(x=Day,y=lg10_immune))+geom_line()
  #To embed special fonts use the cairo graphics library
  
  dat = merge(y.boot,x.boot,)
  lin.mod = lm(y.boot,x.boot)
  
  plot(cv.boot)
}

#Points and smoothing function both follow the definitions in the main call to
#ggplot2

ptw2kegg = data.frame(do.call('rbind', strsplit(system('curl -q -O- http://rest.kegg.jp/link/pathway/hsa',intern=T), '\t', fixed=TRUE)),stringsAsFactors=F)

#Remedy this by specifying levels inside the geom if you don't way subsequent
#geoms to inherit the grouping behavior

#Adding overriding and removing aesthetics
#To override aesthetic=specify a different aesthetic in a downstream call
#Download gene ontology database
library(biomaRt)

ensembl = useMart("ensembl",dataset = "hsapiens_gene_ensembl")

goBP = getBM(attributes = c("ensembl_gene_id","go_id","go_term_name"),
             filters="go_domain",
             values = "biological_process",
             mart=ensembl)


