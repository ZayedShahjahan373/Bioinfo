#############################################################################################################
# 
# This script downloads the reactome database and randomly samples a set of pathways 
# from the reactome database 
#
# Program:  init-SeqNet-sim-main.r
# Version:  
# Author:   zshahjahan
# Purpose:  Initializes important global variables. This includes knitR setup.
# Input:    
# Output:  	N/A
#############################################################################################################
suppressMessages(source(file.path("C:/Users/zayed/OneDrive/Desktop/Research_v2.0/SeqNet_playthru",
                                  "source",
                                  "r",
                                  "init-SeqNet-sim-main.r")))

##
reactome.df = read.csv('http://www.reactome.org/download/current/Ensembl2Reactome.txt',
                       header=F,
                       sep='\t',
                       stringsAsFactors=F)

##The reactome database, unfortunately, does not contain any column information
##apply column names
colnames(reactome.df) = c('ensembl_id',
                      'ptw_short',
                      'ptw_web',
                      'ptw_desc',
                      'tmp',
                      'species')

## subset to human only reactrome pathways
reactome.h_sap.df = reactome.df[reactome.df$species=='Homo sapiens',]

## remove white spaces from ptw_desc
reactome.h_sap.df$ptw_desc = gsub('^ | $','',reactome.h_sap.df$ptw_desc)

#I am expecting to get 18 hits on the string_detect and filter
#as expected
reactome.h_sap.df %>% filter(str_detect(ptw_desc,"glucuronosyltransferases")) %>% View()

homo.ptw.list = tapply(reactome.h_sap.df$ensembl_id, #Take the ensembl gene IDs
                       reactome.h_sap.df$ptw_desc,   #group these by pathway_desc
                       list)                         #Store as a list
#Store the pathway gene list as an RDS object in the analysis directory
h.ptw.list.outfile = file.path(anaDir,"reactome_pathway_list.RDS")
saveRDS(homo.ptw.list,h.ptw.list.outfile)

