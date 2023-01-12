#######################################
# Capscan Small Intestine Capsule Study 
# 
# Evaluation of capsule opening in  
# different locations along the small
# intestine from proximal to distal
#######################################


# load libraries
library(here)
library(ggtext)
library(ggh4x)
library(ggbeeswarm)
library(ggpubr)
library(ggvenn)
library(ggalt)
library(ggrepel)
library(ggplot2)
library(org.Hs.eg.db)
library(missMDA)
library(PCAtools)
library(limma)
library(viridis)
library(corrr)
library(corrplot)
library(vegan)
library(DESeq2)
library(magick)
library(cowplot)
library(phyloseq)
library(data.table)
# library(DivNet)
library(readxl)
library(Hmisc)
library(rstatix)
library(wesanderson)
library(RColorBrewer)
library(reshape2)
library(ropls)
library(topicmodels)
library(alto)
library(gridExtra)
library(tidyverse)


#--------------------------------------------
# source base functions  
#--------------------------------------------
source(paste0(here::here(), "/0-base-functions.R"))


#--------------------------------------------
# define directories
#--------------------------------------------
# raw_data_dir = paste0(here::here(), "/data/raw_data/")
# sherlock_dir = paste0(raw_data_dir, "fromSherlock/")
data_dir = paste0(here::here(), "/data/data/")
fig_dir = paste0(here::here(),"/4-figures/")
fig_dir_main_subpanels = paste0(here::here(),"/4-figures/1-main-figure-subpanels/")
fig_dir_ed_subpanels = paste0(here::here(),"/4-figures/2-extended-data-subpanels/")
tab_dir = paste0(here::here(),"/6-tables/")
results_dir = paste0(here::here(), "/results/")



# #--------------------------------------------
# # define raw data paths
# #--------------------------------------------
# sequencing_ids = paste0(raw_data_dir, "metagenomics_sample_list.xlsx")
# subj_capsule = paste0(raw_data_dir, "variables_all_metagenomics.xlsx")
# sample_dates = paste0(raw_data_dir, "EB-01 Data management - with DNA and pics - protected.xlsx")
# raw_reads =  paste0(raw_data_dir,"fromSherlock/raw_counts.txt")
# trim_reads =  paste0(raw_data_dir,"fromSherlock/trimmed_counts.txt")
# filt_reads =  paste0(raw_data_dir,"fromSherlock/filtered_counts.txt")
# 
# midas_species = paste0(raw_data_dir, "fromSherlock/species_profile_all.txt")

#--------------------------------------------
# define clean data paths
#--------------------------------------------
# sample_read_list = paste0(data_dir, "sampleReadList.txt")
# samples_species = paste0(data_dir, "samples_species.RDS")
sample_data = paste0(data_dir, "sample_data.RDS")

clean_phyloseq_object = paste0(data_dir, "phyloseq_alphaDiv.rds")
# clean_phyloseq_object_transformed = paste0(clean_data_dir, "phyloseq_alphaDiv_transformed.rds")
# clean_phyloseq_object_transformed_filtered = paste0(clean_data_dir, "phyloseq_transformed_filtered.rds")
phyloseq_bilesalt = paste0(data_dir, "phyloseq_bilesalt.rds")
# phyloseq_bilesalt_relabs = paste0(data_dir, "phyloseq_bilesalt_relabs.rds")

#--------------------------------------------
# Set white theme for all figures
#--------------------------------------------
theme_set(theme_minimal()+
            theme(
              panel.background = element_rect(fill = "transparent", colour = NA), # bg of the panel
              plot.background = element_rect(fill = "transparent", colour = NA), # get rid of legend bg
              legend.background = element_rect(fill = "transparent", colour = NA), 
              legend.key = element_rect(fill = "transparent", color = NA),
              legend.text=element_markdown(size=14),
              legend.title = element_text(size = 14),
              axis.text = element_text(size = 14),
              axis.title = element_text(size = 14),
              strip.text = element_text(size = 14)))


#--------------------------------------------
# Set colors for capsule types & stool
#--------------------------------------------
col<-c('#2364AA','#3DA5D9','#134611','#90D733','#743C2F')
CapTypeAndStoolColors <- col
CapType <- c(col[1:4])
CapAndStoolColors <- c(col[1],col[5])


#--------------------------------------------
# define analysis parameters 
#--------------------------------------------
seq_depth_threshold_meta = 750000
seq_depth_threshold_16s = 2500
dna_conc_threshold = 0

fig_padding = 15


##This is probably obsolete now that we have the variables hardcoded into the sample_data but should verify
#Sample to drop due to low sequencing depth 
meta_drop <- c('P17_1', 'O5_1', 'O15_1', 'N17_1', 'M4_2', 'M20_2', 'M19_1', 'M16_2', 'K1_1', 'K21_2', 'K19_2', 'K17_2', 'K11_1', 'J7_2', 'J1_2', 'J17_2', 'I13_1', 'H5_2', 'H13_2', 'F5_2', 'F23_2', 'F21_2', 'F1_2', 'E7_1', 'E21_1', 'E18_2', 'E16_2', 'E12_2', 'E10_2', 'D5_2', 'D19_1', 'D11_2', 'C3_1', 'C1_1', 'B15_1', 'A23_1', 'A21_1', 'A13_2', 'A14_2')
