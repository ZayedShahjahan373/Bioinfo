#############################################################################################################
# 
# This script downloads the gene annotations from ENSEMBL and pre-process these
# For downstream analysis since, we don't need all of the information thats provided
# in the annotation files
#
# Program:  init-SeqNet-sim-main.r
# Version:  
# Author:   zshahjahan
# Purpose:  Initializes important global variables. This includes knitR setup.
# Input:    
# Output:  	N/A
#############################################################################################################
suppressMessages(source(file.path("C:/Users/zayed/OneDrive/Desktop/Research_v2.0/SeqNet_playthru","source","r","init-SeqNet-sim-main.r")))

ensembl.dta.infile = file.path(dtaDir,"annot","Homo_sapiens.GRCh38.104.gtf.gz")

#Check if the annotation dataset has already been downloadded
if(!file.exists(ensembl.dta.infile)){
  # Step 1: Download the Ensembl GTF file for GRCh38 human genome assembly
  print("The annotation dataset is missing in the data directory location. Downloading now...")
  ensembl_url = "ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz"
  download.file(ensembl_url, ensembl.dta.infile, mode = "wb")
} else {
print("Annotation dataset exists in the data directory, proceeding to processing step...")
}


gene.annot.df = rtracklayer::readGFF(ensembl.dta.infile)

#A priori, we know that we are working with genomic regions
#Keep 

keep.cols = c("type","start",
              "end","gene_id",
              "gene_name")

gene.annot.lite.df = dplyr::select(gene.annot.df, all_of(keep.cols))
#Create a variable substracting start from end in order to get gene length
gene.annot.lite.df = gene.annot.lite.df %>% mutate(gene_length=end-start)

#Save the resulting gene.annot.lite.df object as a compressed csv file in the 
#analysis directory for downstream use
annot.outfile = file.path(anaDir,"annotation_processed.tab")

write.table(gene.annot.lite.df,annot.outfile,quote=F,row.names=F,sep='\t')
R.utils::gzip(annot.outfile,overwrite=TRUE)

