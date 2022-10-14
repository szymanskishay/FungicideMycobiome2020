#Create OBJECT
#Shay Szymanski
######
#Phase 0: Load Packages (not all used immediately)
library(colorspace)
library(stringi)
library(rhdf5)
library(zlibbioc)
library(S4Vectors)
library(phyloseq)
library(Biostrings)
library(yaml)
library(colorspace)
library(ggplot2)
library(indicspecies)
library(vegan)
library(decontam)
library(data.table)
library(dplyr)
library(ggpubr)
library(tidyverse)
##### 
#Phase 1: Create Object
options(scipen = 999) # sets the number of decimals displayed on screen 
set.seed(9279) # setting the seed for reproducibility so analyses can be reproduced [9279]
setwd("C:/Users/sakuy/Desktop/PhytobiomesPaper/")
fungi_otus<- read.delim("CORE DATAS/OTU_TABLE_RELABEL.txt",
                        row.names=1)  # add otu table as a datframe
fungi_otus_phy <-otu_table(fungi_otus,
                           taxa_are_rows = TRUE) # formats the dataframe so it is interpretable by phyloseq
fungi_metadata <-read.delim("CORE DATAS/blueberry_map_fix_COPY.txt",
                            row.names=1) # read in the mapping file with your samples metadata to a dataframe
fungi_metadata_phy <-sample_data(fungi_metadata) # put into phyloseq format
fungi_taxonomy<-read.delim("CORE DATAS/constax_taxonomy_ITS_2-1-1.txt",
                           header=TRUE, 
                           row.names=1) # read in taxonomy information to dataframe
fungi_taxonomy_phy <- tax_table(as.matrix(fungi_taxonomy)) # put into phyloseq readable format
fungi_sequences <- readDNAStringSet("CORE DATAS/otus_R1.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
physeq_object_spray <- phyloseq(fungi_otus_phy, fungi_metadata_phy, fungi_taxonomy_phy, fungi_sequences)  # combine all portions into one single phyloseq object
df.mainseq<-as.data.frame(sample_data(physeq_object_spray))# Convert to data frame for exporting
write.csv(df.mainseq, file="Generated Tables/verifywhole.csv") #creates an output of the pure phyloseq object

# Check Negative Controls ---------------
df_spray <- as.data.frame(sample_data(physeq_object_spray)) # Put sample_data into a ggplot-friendly data.frame
df_spray$LibrarySize_spray <- sample_sums(physeq_object_spray) #Add a column for Library Size
df_spray <- df_spray[order(df_spray$LibrarySize_spray),] #sorts
df_spray$Index <- seq(nrow(df_spray)) #Adds another layer of order
write.csv(df_spray, file = "Generated Tables/rank_sums_spray.csv") #another useful output
#the following code allows us to abolish bad samples. The above file shows one true sample with 5 reads, so we're cutting it.
otu_table(physeq_object_spray) <- subset(otu_table(physeq_object_spray), select = -c(UW3B33)) #uw3b33 is a bad sample being cut here
otu_table(physeq_object_spray) <- subset(otu_table(physeq_object_spray), select = -c(GWB22)) #GWB22 is a bad sample being cut here (mostly plant reads)
otu_table(physeq_object_spray) <- subset(otu_table(physeq_object_spray), select = -c(GWB21)) #GWB2 is a bad sample being cut here (mostly plant reads)

sample_data(physeq_object_spray)$is.neg <- sample_data(physeq_object_spray)$Sample_or_Control == "Control Sample" #Assigns control samples to be a negative column
contamdf.prev_spray <- isContaminant(physeq_object_spray, method="prevalence", neg="is.neg") #checks what is prevalent in the controls and samples to cut potential contaminants
df.contam<-as.data.frame(sample_data(contamdf.prev_spray)) #convert to data frame for exporting
write.csv(df.contam, file="Generated Tables/contaminants.csv") #this prints us a table of contaminants in the working directory
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_spray <- transform_sample_counts(physeq_object_spray, function(abund) 1*(abund>0))
ps.pa.neg_spray <- prune_samples(sample_data(ps.pa_spray)$Sample_or_Control == "Control Sample", ps.pa_spray) #Creates an object of just Controls
ps.pa.pos_spray <- prune_samples(sample_data(ps.pa_spray)$Sample_or_Control == "True Sample", ps.pa_spray) #Creates an object of true samples
# Make data.frame of prevalence in positive and negative samples
df.pa_spray <- data.frame(pa.pos_spray=taxa_sums(ps.pa.pos_spray), pa.neg_spray=taxa_sums(ps.pa.neg_spray),
                          contaminant=contamdf.prev_spray$contaminant)
write.csv(df.pa_spray, file="Generated Tables/contaminant-sums.csv") #again just creating a table for this to export
# remove contaminants
# add OTU_3 & others to contaminants (Blueberry DNA)
contamdf.prev_spray["OTU_13",]$contaminant<-as.logical("TRUE") #Vaccinium angustifolium ITS
contamdf.prev_spray["OTU_106",]$contaminant<-as.logical("TRUE") #Vaccinium corymbosum
contamdf.prev_spray["OTU_210",]$contaminant<-as.logical("TRUE") #Panicum virgatum
contamdf.prev_spray["OTU_160",]$contaminant<-as.logical("TRUE") #Vaccinium
contamdf.prev_spray["OTU_275",]$contaminant<-as.logical("TRUE") #Vaccinium
contamdf.prev_spray["OTU_301",]$contaminant<-as.logical("TRUE") #Vaccinium
contamdf.prev_spray["OTU_304",]$contaminant<-as.logical("TRUE") #Vaccinium
contamdf.prev_spray["OTU_718",]$contaminant<-as.logical("TRUE") #Panicum virgatum
contamdf.prev_spray["OTU_723",]$contaminant<-as.logical("TRUE") #Panicum virgatum
contamdf.prev_spray["OTU_361",]$contaminant<-as.logical("TRUE") #vaccinium its
contamdf.prev_spray["OTU_901",]$contaminant<-as.logical("TRUE") #Crepis capillaris
contamdf.prev_spray["OTU_524",]$contaminant<-as.logical("TRUE") #vaccinium ITS
contamdf.prev_spray["OTU_983",]$contaminant<-as.logical("TRUE") #Triticum sp
contamdf.prev_spray["OTU_533",]$contaminant<-as.logical("TRUE") #Vaccinium ITS
contamdf.prev_spray["OTU_547",]$contaminant<-as.logical("TRUE") #vaccinium ITS
contamdf.prev_spray["OTU_1081",]$contaminant<-as.logical("TRUE") #Embryophyte/Plantago DNA
contamdf.prev_spray["OTU_181",]$contaminant<-as.logical("TRUE") #Embryophyte/Plantago lancelota
contamdf.prev_spray["OTU_24",]$contaminant<-as.logical("TRUE") #?????


ps.noncontam_spray <- prune_taxa(!(contamdf.prev_spray$contaminant), physeq_object_spray) #keeps everything thats not labeled contaminant
ps.spray<-ps.noncontam_spray #keeps original as backup
otu_table(ps.spray) <- otu_table(ps.spray)[which(rowSums(otu_table(ps.spray)) >= 1),] #only keeps OTUs with reads
ps.spray.naut<- subset_samples(ps.spray, origin%in%c("Blueberry")) #retains complete data set minus controls!
spraydata<- subset_samples(ps.spray.naut, Extraction%in%c("Kingfisher", "Chelex")) #retains only used samples
