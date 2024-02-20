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
#read in the downloaded dataset
#We know a priori that the dataset is stored in a tsv format
#check all the files 
raw_data_files = list.files(dtaDir,pattern = "\\.tsv$")
#Use fread to read in very large rna_seq datasets
liver_rna.df = fread(file.path(dtaDir,raw_data_files))
#Convert this to a dataframe since those are easier to work with than data tables
liver_rna.df = as.data.frame(liver_rna.df)

#read in the reactome pathway list and the gene annotations dataframe, the goal
#now will be to first, We are going to select pathways related to the human liver
#for downstream simulation step. Cause we want to simulate processes occuring in
#the liver and how those can be modulated in silico
reactome.list.infile = file.path(anaDir,"reactome_pathway_list.RDS")
reactome.list = readRDS(reactome.list.infile)

pathway.names = names(reactome.list)


annot.df.infile = file.path(anaDir,"annotation_processed.tab.gz")
annot.df = read.table(annot.df.infile,header = TRUE)


#Get the sample IDs and Gene IDs
sample_IDs = colnames(liver_rna.df)[colnames(liver_rna.df)!="genes"]
gene_IDs   = liver_rna.df[,"genes"]

#set the rownames equal to the gene_IDs so that these are preserved downstream
rownames(liver_rna.df) = gene_IDs
total.genes = nrow(liver_rna.df)
total.samples = ncol(liver_rna.df)

#Transpose the subsetted gene expression matrix
liver_rna_tr.df = as.data.frame(t(liver_rna.df)[-1,])
#This is the data matrix we will  be working off of
gene.exp.df = as.data.frame(lapply(liver_rna_tr.df, as.numeric),
                            row.names = rownames(liver_rna_tr.df))

#For each liver process pathway do:
lapply(1:length(liverPtwVect), function(i){
  liver.proc = liverPtwVect[i]
  
  liver.ptw.names = pathway.names[str_detect(pathway.names,liver.proc)]
  
  liver.ptw.list = reactome.list[liver.ptw.names]
  #From the resultant liver.ptw.list object retrieve the unique ENSEMBL
  #gene IDs; There are duplicated ensembl gene names
  liver.ensembl.names = unique(unlist(liver.ptw.list))
  #Now get the ensembl gene symbols for the unique gene_IDs
  annot.liv.filtered = annot.df %>% 
    filter(gene_id %in% liver.ensembl.names & type=="gene")
  #There were 62 genes selected but the liver.ensembl names contained
  #aprroximately 400 symbols, the discrepancy comes from there being
  #other types of IDs stored in the annotation files
  
  #Subset the gene expression matrix to only keep the genes in these pathways
  keep.cols = colnames(gene.exp.df)[colnames(gene.exp.df) %in%  annot.liv.filtered$gene_name]
  gene.exp.ptw.df = gene.exp.df[,keep.cols]
  
  total.mapped.reads = colSums(gene.exp.ptw.df)
  gene.lengths = annot.liv.filtered[annot.liv.filtered$gene_name %in% keep.cols,"gene_length"]
  
  ## convert the data into an edgeR object
  dge.spc = DGEList(counts=gene.exp.ptw.df,
                    lib=colSums(gene.exp.ptw.df),
                    group=rep(1:ncol(gene.exp.ptw.df)))
  #After rpkm is performed the count matrix is no loger a dgelist
  #this should be rectified before performing tmm normalizaion
  dge.cpm = DGEList(counts=cpm(dge.spc,gene.length = gene.lengths),
                     lib=dge.spc$samples$lib.size,
                     norm.factors = dge.spc$samples$norm.factors)
  
  
  #There are a few more programming subtleties here that might require
  #more careful thinking
  #Calculate TMM-scaling factors
  dge.cpm = calcNormFactors(dge.cpm)
  
  ##Save the resulting dge objects for downstream exploratory
  ##Analysis; there should be a subdirectory corresponding to the module
  dge.outdir = file.path(anaDir,"00")
  dge.outfile.name = paste0(liver.proc,"_","dge_cpm.RDS")
  saveRDS(dge.cpm,file.path(dge.outdir,dge.outfile.name))
  
  #Save the tmm scaling factors for further analysis
  tmm.factors = data.frame(samid=rownames(dge.cpm$samples),
                           total_mapped=dge.cpm$samples$lib.size,
                           norm_factors = dge.cpm$samples$norm.factors)
  #Save this sets of tmm factors for use in 
  tmm.fac.outfile.name = paste0(liver.proc,"_","tmm_factor.tab")
  write.table(tmm.factors,file.path(dge.outdir,tmm.fac.outfile.name))
  R.utils::gzip(file.path(dge.outdir,tmm.fac.outfile.name))
  
  #TODO: For next time, calculate and save the moderated log2 counts per
  #million before and after TMM-normalization
})


