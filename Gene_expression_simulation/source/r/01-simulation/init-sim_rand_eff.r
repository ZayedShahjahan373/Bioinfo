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

#Script to simulate antibody response as a function of gene expression but for
#now do so only for the bile synthesis genes
infile.name = "bile_sim_gnexp.RDS"
gnexp.mat.list = readRDS(file.path(gnexp.mat.indir,infile.name))

sim.gene.exp = gnexp.mat.list[["x"]]

n.genes = ncol(sim.gene.exp)
n.samps = nrow(sim.gene.exp)

gene.sym = colnames(sim.gene.exp)
#Get 301 subjects, that way there will be three timepoints per subject
n.subjects = n.samps/3

plb.df = data.frame(SAMID=rownames(sim.gene.exp),
                    PATID=paste0("PAT:",sample(1:n.subjects,
                                 size = n.samps,
                                 replace = T)))

uniq.PATID = unique(plb.df$PATID)
#Assign the timepoints for which the samples are coming from 
plb.df=plb.df %>% 
  group_by(PATID) %>% 
  mutate(tmpt = paste0("T",row_number()))

#Perform subject level randomization for covariates
randomized_covar= lapply(1:length(covariateObjList), function(i){
  cvt = covariateObjList[[i]]
  cvt.name = names(covariateObjList)[i]
  
  randomized.cvt = data.frame(uniq.PATID,
                              cvt.name=sample(cvt,length(uniq.PATID),replace = TRUE))
  #Rename the columns according to the covariates; each run of this function will
  #output the same 2-column dataset rename the second column according to the 
  #name of the covariate
  colnames(randomized.cvt)[2] = cvt.name
  return(randomized.cvt)
})
patient.data = Reduce(function(x,y) merge(x,y,by="uniq.PATID"),randomized_covar)
patient.data[,"AGE"] = round(rnorm(dim(patient.data)[1],45,5))


colnames(patient.data)[colnames(patient.data)=="uniq.PATID"]="PATID"

#Construct the patient-specimen metadata
mta.df = merge(plb.df,patient.data,by="PATID",all.x=T)


#Create the outcome variable; Or even better, create a set of outcome variables
#Each using different randomly selected subset of the genes in your expression
#matrix
#First simulate one immunogenicity variable
#Specify patient and covariate specific random effects
patient_effect = rnorm(length(unique(mta.df$SAMID)),0,2)
names(patient_effect) = mta.df$SAMID
site_effect = rnorm(length(unique(mta.df$treatSite)),0,0.75)
names(site_effect)=unique(mta.df$treatSite)
#Ethnicity effect following a laplace distribution
ethn_effect = rlaplace(length(unique(mta.df$ethnVar)),0,0.2)
names(ethn_effect)=unique(mta.df$ethnVar)
#transplant_effect = 
#Specify the model matrix
design_mat = model.matrix(~AGE + sexVar + factor(treatType) + factor(ethnVar) + factor(tranplType),data = mta.df)

sim_immu.df = matrix(0,nrow = nrow(design_mat),ncol=n.genes)

#Randomly select a group genes and assign causality to them; assuming 20% of our
#genes are causally linked to the trait
cau.gene = sample(gene.sym,size=floor(0.2*length(gene.sym)))

#Now we create our immunogenicity variable
immune.vect.list = lapply(1:length(gene.sym), function(j){
  gene = gene.sym[j]
  gene.exp.vect = sim.gene.exp[,gene]
  #Add interaction effects
  interaction_effect = ethn_effect[match(mta.df$ethnVar,unique(mta.df$ethnVar))]*site_effect[match(mta.df$treatSite,unique(mta.df$treatSite))]
  names(interaction_effect)=NULL
  
  if(gene %in% cau.gene){
    causal.effect = rnorm(1,5,1)
  } else {
    causal.effect = 0 #Since only a random subset of genes will be causal
  }
  
  #This is essentially or finished product for this script
  immunogenicity = patient_effect + gene.exp.vect  + site_effect[match(mta.df$treatSite,unique(mta.df$treatSite))] + ethn_effect[match(mta.df$ethnVar,unique(mta.df$ethnVar))]+interaction_effect+causal.effect*gene.exp.vect
  names(immunogenicity) = names(patient_effect)
  return(immunogenicity)
})
names(immune.vect.list) = gene.sym

immune.df = as.data.frame(immune.vect.list)

imm.df = as.data.frame(rowSums(immune.df))
imm.df[,"SAMID"] = rownames(imm.df)
colnames(imm.df)[colnames(imm.df)=="rowSums(immune.df)"] = "immune_response"

#Save the resulting response variable
imm.df.outfilename="immunology_outcome.csv"
write.csv(imm.df,file.path(anaDir,imm.df.outfilename),row.names=F)

#Save the metadata
mta.df.infile = "project_metadata.csv"
write.csv(mta.df,file.path(anaDir,mta.df.infile))
