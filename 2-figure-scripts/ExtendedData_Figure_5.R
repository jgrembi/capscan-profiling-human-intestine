#######################################
# Capsule Study - Gen1 
#
# Extended Data Figure 5
#
# Author: Rebecca Culver
#######################################

rm(list=ls())
# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

#######################################
# ED Fig. 5a,b - Abundance of reads mapped to CAZyme database for location and by subject
#######################################

# Import files
vars <- data.frame(sample_data(readRDS(clean_phyloseq_object))) %>% select(Main, Subject, Type, Set, meta_samplename, location, drop_meta)
cazyme_dir <- "~/Dropbox/Capsule/SECOND REVISION/repo_files/cazyme/"
caz_mapping <- read.delim(paste0(cazyme_dir, 'allCazymesFrom1kbContigsMappedToSamples.txt'), header=FALSE) %>%
  dplyr::rename(meta_samplename=V1) %>%
  mutate(cazyme_perc_mapped = as.numeric(str_split(V2, '%', simplify=TRUE)[,1])) %>% 
  select(-V2) %>%
  left_join(vars, by='meta_samplename') %>%
  mutate(location = ifelse(location == "Capsule", "Devices", as.character(location)))

df2plot <- caz_mapping %>%
  filter(drop_meta == FALSE, location %in% c('Devices','Stool'), Set %in% c('2','3','4','5','Stool'))

# Figure generation by location
(ed_fig_5a<-ggplot(df2plot, aes(x=location, y=cazyme_perc_mapped)) +
    geom_boxplot() +
    geom_point(position='jitter') + 
    labs(x = "", y = "% reads mapped to CAZyme database")) 

ggsave(paste0(fig_dir_ed_subpanels, 'ED_Fig_5a_cazyme_capsule_v_stool_compared.pdf'), ed_fig_5a, width=4, height=4)

# Figure generation by location and subject
(ed_fig_5b <- ggplot(df2plot, aes(x=reorder(Subject,cazyme_perc_mapped), y=cazyme_perc_mapped)) +
    geom_boxplot() +
    geom_point(position='jitter') +
    facet_wrap(~location, ncol=1) + 
    labs(x = "Subject", y = "% reads mapped to CAZyme database"))

ggsave(paste0(fig_dir_ed_subpanels, 'ED_Fig_5b_cazyme_supjects_compared.pdf'), ed_fig_5b, width=6, height=4)


#######################################
# ED Fig. 5c,d - ASVs/Families correlated with % reads mapped to CAZyme database
#######################################

### LOOKING AT THE ASV-LEVEL IN THE DEVICE SAMPLES ONLY


# Import taxa 16s data
taxa_devices <- readRDS(clean_phyloseq_object) %>%
  subset_samples(!drop_meta & Set %in% c("2", "3", "4", "5")) %>%
  filter_taxa(., function(x) sum(x > 3) > (0.05*length(x)), TRUE) %>% 
  transform_sample_counts(function(x) {log2(x+1)}) %>%
  psmelt() %>%
  select(ASV, Abundance, Phylum, Family, Genus, Species, meta_samplename) 

head(taxa_devices)

# Filter OTUs
capsule_asv_cazyme <- taxa_devices %>%
  left_join(df2plot, by='meta_samplename') %>%
  filter(Abundance > 0)

df_c.filt <- asv_corr(capsule_asv_cazyme, location_test = "Devices", cor_variable = 'cazyme_perc_mapped', include_0_abundance = F) %>%
  filter(adj_pval < 0.001) 

# Filter dataframe to only include significant ASVs              
df_ed5c <- capsule_asv_cazyme %>%
  filter(ASV %in% df_c.filt$ASV) %>%
  mutate(label = ifelse(!is.na(Genus), paste0(Genus, " ", Species), Family))


(ed_5c <- ggplot(df_ed5c, aes(x=cazyme_perc_mapped, y=Abundance)) + 
    geom_point(alpha = 0.6) + 
    geom_smooth(method=lm, se=FALSE) +
    labs(y=expression('log'[2]*'(ASV count)'), 
         x=expression('% reads mapped to CAZyme database')) +
    stat_cor(method = "spearman", 
             aes(label = after_stat(r.label)),
             cor.coef.name = "rho",
             r.accuracy = 0.01,
             label.x = 5.5,
             label.y = 15) +
    facet_wrap(~label+ ASV, ncol = 4) +
    theme(strip.text = element_text(size=12),
          axis.text = element_text(size = 12)))

ggsave(paste0(fig_dir_ed_subpanels, 'ED_Fig_5c_sig_corr_cazyme.pdf'), plot = ed_5c, width=8, height=8)


### LOOKING ONLY AT THE STOOL, NO SIGNIFICANT CORRELATIONS AT THE ASV LEVEL SO THIS SCRIPT GROUPS ASVS BY FAMILIES
# 
# taxa <- readRDS(clean_phyloseq_object) %>%
#     subset_samples(!drop_meta & Set %in% c("Stool")) %>% #, "Stool"
#     filter_taxa(., function(x) sum(x > 3) > (0.05*length(x)), TRUE) %>% # Gets rid of all taxa not found at a count of 3 in at least 5% of samples (that's 14 samples)
#     transform_sample_counts(function(x) {log2(x+1)}) %>%
#     psmelt() %>%
#     select(ASV, Abundance, Phylum, Family, Genus, Species, meta_samplename) %>%
#     left_join(df2plot, by='meta_samplename') %>%
# ----> #     group_by(Abundance, Phylum, Family, meta_samplename, cazyme_perc_mapped) %>%
# |  #     dplyr::summarise(Abundance=sum(Abundance)) %>%
# |  #     dplyr::rename(ASV=Family)
# |  #              
# ---Becca, above at the arrow is the group_by command that I copied below (line 135) and then removed 'Abundance'.  
#     I'm not 100% sure this is correct, but I think it is because I think you want to sum all Abundance values for the Family. 
#     Take a look and make sure I'm understanding the code correctly.
taxa_stool_family <- readRDS(clean_phyloseq_object) %>%
  subset_samples(!drop_meta & Set %in% c("Stool")) %>% #, "Stool"
  tax_glom("Family") %>%
  filter_taxa(., function(x) sum(x > 3) > (0.05*length(x)), TRUE) %>%
  transform_sample_counts(function(x) {log2(x+1)}) %>%
  psmelt() %>%
  select(Abundance, Phylum, Family, meta_samplename) 

# Filter OTUs
stool_asv_cazyme <- taxa_stool_family %>%
  left_join(df2plot, by='meta_samplename') %>%
  dplyr::rename(ASV=Family) %>%
  filter(Abundance > 0)

df_d.filt <- asv_corr(stool_asv_cazyme, location_test = "Stool", cor_variable = 'cazyme_perc_mapped', include_0_abundance = F) %>%
  filter(adj_pval < 0.01) 
# Filter dataframe to only include significant ASVs              
df_ed5d <- stool_asv_cazyme %>%
  filter(ASV %in% df_d.filt$ASV) %>%
  mutate(label = ASV)

length(unique(df_ed5d$ASV))

(ed_5d <-ggplot(df_ed5d, aes(x=cazyme_perc_mapped, y=Abundance)) + 
    geom_point()  + geom_smooth(method=lm, se=FALSE) +
    labs(y=expression('log'[2]*'(ASV count)'), x=expression('% reads mapped to CAZyme database')) +
    stat_cor(method = "spearman", 
             aes(label = after_stat(r.label)),
             cor.coef.name = "rho",
             r.accuracy = 0.01,
             label.x = 5.15,
             label.y = 13) +
    facet_wrap(~ ASV)+
    theme(strip.text = element_text(size=12),
          axis.text = element_text(size = 12)))

ggsave(paste0(fig_dir_ed_subpanels, 'ED_Fig_5d_sig_corr_family_stool.pdf'), ed_5d, width=6, height=6)

#######################################
# ED Fig. 5e - Number of CAZymes detected in isolates from the human intestinal tract (not stool) of Subject 1
#######################################

# Import RGI output files

### CHANGE THIS PATH HERE
files <- list.files(path=paste0(cazyme_dir,'cazyme_rawoutput'), pattern="*overview*", full.names=TRUE)

df<-lapply(files, function(x) read.delim(x) %>% mutate(file_name = str_split(x, '.overview.txt', simplify=TRUE)[,1]))
df<-map_dfr(df, ~.x %>% mutate(across(everything(), as.character)))

### CHANGE THIS LINE HERE BASED ON WHERE FOLDER LOCATION IS - you'll want to modify the index (currently 12) to be whichever index the file_name is in the folder path
df<-bind_rows(df) %>% 
  mutate(file_name = str_split(file_name, '/', simplify=TRUE)[,10]) %>%
  filter(X.ofTools > 1) 

# Merge in taxa information
ani<-read.table(file=paste0(cazyme_dir,'gtdbtk.bac120.summary.tsv'),sep = '\t', header = TRUE) %>%
  mutate(Sample = str_split(user_genome, '.f', simplify=T)[,1]) %>%
  mutate(Species = str_split(classification, ';', simplify=T)[,7]) %>%
  mutate(Species = ifelse(Species == 's__', str_split(classification, ';', simplify=T)[,6], Species)) %>%
  mutate(Phylum = str_split(classification, ';', simplify=T)[,2]) %>%
  mutate(Family = str_split(classification, ';', simplify=T)[,5]) %>%
  mutate(Genus = str_split(classification, ';', simplify=T)[,6]) %>% 
  mutate(Genus = str_split(Genus, '__', simplify=T)[,2]) %>%
  mutate(Family = str_split(Family, '__', simplify=T)[,2]) %>%
  mutate(Phylum = str_split(Phylum, '__', simplify=T)[,2]) %>%
  mutate(Family = ifelse(Family == 'Bacillaceae_G','Bacillaceae',Family)) %>%
  mutate(Phylum = ifelse(Phylum == 'Firmicutes_C','Firmicutes',Phylum)) %>%
  mutate(Phylum = ifelse(Phylum == 'Actinobacteriota','Actinobacteria',Phylum)) %>%
  mutate(Phylum = ifelse(Phylum == 'Bacteroidota','Bacteroidetes',Phylum)) %>%
  select(user_genome, fastani_ani, classification, Sample, Species:Genus) %>%
  dplyr::rename(file_name=Sample) %>%
  right_join(df, by="file_name") %>%
  filter(is.na(Phylum)==FALSE)

# Create an f.data where there's a cazyme column that includes the final name
tmp<-ani
tmp$cazyme<-tmp$DIAMOND
tmp<-tmp %>%
  mutate(cazyme = case_when(
    !(DIAMOND=='-') ~ DIAMOND,
    endsWith(DIAMOND, "-") ~ HMMER,
    endsWith(HMMER, "-") ~ DIAMOND
  ))

# Remove stuff in parantheses
df.final<- tmp %>%
  mutate(cazyme = gsub(" *\\(.*?\\) *",'',cazyme)) 
length(df.final$cazyme)
length(unique(df.final$cazyme))

df.final<-df.final %>%
  mutate(description = case_when(
    grepl('GH',cazyme) ~ 'Glycoside Hydrolases',
    grepl('GT',cazyme) ~ 'Glycoside Transferases',
    grepl('PL',cazyme) ~ 'Polysaccharide Lysases',
    grepl('CE',cazyme) ~ 'Carbohydrate Esterases',
    grepl('AA',cazyme) ~ 'Auxiliary Activities',
    grepl('CBM', cazyme) ~ 'Carbohydrate-Binding Molecules')) %>%
  # Remove contaiminated bifido
  filter(!user_genome %in% 'demux_S_17.filtered.scaffolds')

write.csv(df.final, paste0(tab_dir, 'Supplemental_Table_4_cazymes.csv'))

df_ed_5e <- df.final %>% 
  group_by(user_genome, Species) %>%
  dplyr::summarise(num_cazymes=n()) %>% 
  arrange(num_cazymes) %>%
  mutate(Species = gsub("_", "-", gsub("g__", "", gsub("s__", "", Species))))

(ed_5e <- ggplot(df_ed_5e, aes(y=reorder(Species, num_cazymes), x=num_cazymes)) +
  geom_point(position='jitter') +
  geom_boxplot(alpha=0) + 
  labs(y = "", x = "Number of CAZymes detected"))

ggsave(paste0(fig_dir_ed_subpanels, 'ED_Fig_5e_num_Cazymes_in_all.pdf'), plot = ed_5e, width=6, height=5)


#######################################
# ED Fig. 5f,g - Abundance of reads mapped to CARD database for location and by subject
#######################################
#--------------------------------
# resset filepath
#--------------------------------
amr_dir <- "~/Dropbox/Capsule/SECOND REVISION/repo_files/amr/"
raw_amr <- read.table(file=paste0(amr_dir, '90perc/summary.txt'), header=FALSE, sep="\t") %>%
  mutate(amr_alignment = as.numeric(str_split(V2, '%', simplify=TRUE)[,1])) %>%
  dplyr::rename(file_name=V1) %>% select(-V2)
head(raw_amr)

vars <- data.frame(sample_data(readRDS(clean_phyloseq_object))) %>% 
  mutate(file_name = meta_samplename) %>%
  select(Main, Subject, Type, Set, file_name, location)

amr.df <- raw_amr %>% 
  left_join(vars, by='file_name') %>% 
  filter(is.na(location)==FALSE) %>%
  filter(!location %in% 'Saliva')
head(amr.df)

df_5fg <- amr.df %>% 
  mutate(location = ifelse(location == "Capsule", "Devices", as.character(location)))

(ed_5f <- ggplot(df_5fg, aes(x=location, y=amr_alignment)) +
  geom_boxplot (outlier.shape=NA) +
  geom_point(position='jitter')+
  stat_compare_means() +
  labs(x = "", y = '% reads mapped to CARD'))

ggsave(paste0(fig_dir_ed_subpanels, 'ED_Fig_5f_amr_denisty_overall.pdf'), plot = ed_5f, width=4, height=4)


(ed_5g <- ggplot(df_5fg, aes(x=reorder(Subject, amr_alignment), y=amr_alignment)) +
  geom_boxplot (outlier.shape=NA) +
  geom_point(position='jitter')+
  labs(x = "", y = '% reads mapped to CARD') + 
  facet_wrap(~location, ncol=1))


ggsave(paste0(fig_dir_ed_subpanels, 'ED_Fig_5g_amr_denisty_subj.pdf'), plot = ed_5g, width=8, height=4)



#######################################
# ED Fig. 5h,i - ASVs/Families correlated with % reads mapped to CARD database
#######################################

### ASV-LEVEL CORRELATIONS FOR JUST THE DEVICE SAMPLES

# Incoporate the abundance of taxa (ie. Escherichia)
# taxa_all <- readRDS(clean_phyloseq_object) %>% 
#   subset_samples(!drop_16s & Set %in% c("2", "3", "4", "5")) %>%
#   filter_taxa(., function(x) sum(x > 3) > (0.05*length(x)), TRUE) %>% # Gets rid of all taxa not found at a count of 3 in at least 5% of samples (that's 14 samples)
#   transform_sample_counts(function(x) {log2(x + 1)}) %>%
#   psmelt() %>%
#   dplyr::rename(file_name = meta_samplename) %>%
#   select(file_name, Sample, Abundance, Phylum, Family, Genus, Species, ASV) %>% 
#   left_join(amr.df, by=c('file_name')) %>% 
#   unique() %>%
#   mutate(location = ifelse(location == "Capsule", "Devices", as.character(location)))
# head(taxa_all)

# Filter OTUs
capsule_asv_amr<- taxa_devices %>%
  dplyr::rename(file_name = meta_samplename) %>%
  left_join(amr.df, by=c('file_name')) %>% 
  unique() %>%
  mutate(location = ifelse(location == "Capsule", "Devices", as.character(location))) %>%
  filter(Abundance > 0)
## Look at distribution of Abundance and amr_alignment
ggplot(capsule_asv_amr) + 
  geom_histogram(aes(x = Abundance))

ggplot(capsule_asv_amr) + 
  geom_histogram(aes(x = amr_alignment))
##Neither are normally distributed so need to use Spearman correlations 

df_h.filt <- asv_corr(capsule_asv_amr, location_test = "Devices", cor_variable = 'amr_alignment', include_0_abundance = F) %>%
  filter(adj_pval < 0.001) 


# Filter dataframe to only include significant ASVs              
capsule_asv_amr.filt <- capsule_asv_amr %>%
  filter(ASV %in% df_h.filt$ASV) %>%
  mutate(label = ifelse(!is.na(Genus), paste0(Genus, " ", Species), Family))

length(unique(capsule_asv_amr.filt$ASV))

(ed_5h <- ggplot(capsule_asv_amr.filt, 
                 aes(x=amr_alignment, y=Abundance)) + 
    geom_point()  + 
    geom_smooth(method=lm, se=FALSE) +
    labs(y=expression('log'[2]*'(ASV count)'), x=expression('% reads mapped to CARD')) +
    stat_cor(method = "spearman", 
             aes(label = after_stat(r.label)),
             cor.coef.name = "rho",
             r.accuracy = 0.01,
             label.x.npc = "center",
             label.y.npc = "top") +
    facet_wrap(~label+ASV, ncol = 4)+
    theme(strip.text = element_text(size=12),
          axis.text = element_text(size = 12)))

ggsave(paste0(fig_dir_ed_subpanels, 'ED_Fig_5h_sig_corr_amr_device.pdf'), ed_5h, width=12, height=8)

### FAMILY-LEVEL CORRELATIONS IN JUST THE STOOL - note that there are no correlations at the ASV-level

# Filter OTUs
stool_asv_amr<- taxa_stool_family %>%
  dplyr::rename(file_name = meta_samplename) %>%
  left_join(amr.df, by=c('file_name')) %>%
  dplyr::rename(ASV=Family) %>%
  filter(Abundance > 0)

df_i.filt <- asv_corr(stool_asv_amr, location_test = "Stool", cor_variable = 'amr_alignment', include_0_abundance = F) %>%
  filter(adj_pval < 0.01) 

# Filter dataframe to only include significant ASVs              
stool_asv_amr.filt <- stool_asv_amr %>%
  filter(ASV %in% df_i.filt$ASV) %>%
  mutate(label = ASV)

head(stool_asv_amr.filt)

(ed_5i <- ggplot(stool_asv_amr.filt, aes(x=amr_alignment, y=Abundance)) + 
  geom_point()  + geom_smooth(method=lm, se=FALSE) +
  labs(y=expression('log'[2]*'(ASV count)'), x=expression('% reads mapped to CARD')) +
  stat_cor(method = "spearman", 
           aes(label = after_stat(r.label)),
           cor.coef.name = "rho",
           r.accuracy = 0.01,
           label.x.npc = "center",
           label.y.npc = "top") +
  facet_wrap(~label, scales = "free")+
  theme(strip.text = element_text(size=12),
        axis.text = element_text(size = 12)))

ggsave(paste0(fig_dir_ed_subpanels, 'ED_Fig_5i_sig_corr_amr_family_stool.pdf'), ed_5i, width=6, height=6)          

#######################################
# ED Fig. 5j - % of reads mapped to CARD database with efflux pumps removed
#######################################
raw_amr_noefflux <- read.table(file=paste0(amr_dir, '90perc_noEfflux/summary_no_efflux.txt'), header=FALSE, sep="\t") %>%
  mutate(amr_alignment = as.numeric(str_split(V2, '%', simplify=TRUE)[,1])) %>%
  dplyr::rename(file_name=V1) %>% 
  select(-V2)

amr.df_noefflux <- raw_amr_noefflux %>% 
  left_join(vars, by='file_name')  %>% 
  filter(is.na(location)==FALSE) %>%
  filter(!location %in% 'Saliva')
head(amr.df_noefflux)

df_5fj <- amr.df_noefflux %>% 
  mutate(location = ifelse(location == "Capsule", "Devices", as.character(location)))

(ed_5j <- ggplot(df_5fj, aes(x=location, y=amr_alignment)) +
    geom_boxplot (outlier.shape=NA) +
    geom_point(position='jitter')+
    stat_compare_means() +
    labs(x = "", y = '% reads mapped to CARD'))

ggsave(paste0(fig_dir_ed_subpanels, 'ED_Fig_5j_amr_denisty_no_efflux_pump.pdf'), plot = ed_5j, width=4, height=4)


#######################################
# ED Fig. 5k - Ratio of AMR to bacterial genes in MAGs
#######################################

# Import number of AMRs found in each MAG
files <- list.files(path=paste0(amr_dir, 'mags_rgi'), pattern="*.txt", full.names=TRUE)
df_amr_mags <- lapply(files, function(x) read.delim(x) %>% 
               mutate(mag=str_split(str_split(x,'mags_rgi/', simplify=TRUE)[,2],'.strict', simplify=TRUE)[,1]))
df_amr_mags <- map_dfr(df_amr_mags, ~.x %>% mutate(across(everything(), as.character)))
num_amr <- bind_rows(df_amr_mags) %>%
  group_by(mag) %>% 
  dplyr::summarise(num_amr = n())

# Import taxonomic information about mags
taxa_amr <- read.table(paste0(amr_dir,'gtdbtk.bac120.summary.tsv'), sep='\t', header=TRUE) %>%
  select(user_genome, classification) %>%
  dplyr::rename(mag=user_genome) %>%
  mutate(Family = str_split(str_split(classification ,'f__', simplify=TRUE)[,2], ';g', simplify=TRUE)[,1]) %>%
  select(-classification)

# Import number of bacterial genes identified in each MAG
mag.genes <- read.csv(paste0(data_dir, 'mags_metadata.csv')) 
# Combine all metadata together
df.final_amr_mags <- mag.genes %>% 
  left_join(taxa_amr, by='mag') %>% 
  left_join(num_amr, by='mag') %>%
  mutate(num_amr = ifelse(is.na(num_amr), 0, num_amr)) %>%
  filter(!is.na(Family))  ##This filters out 2 NA MAGs (un-annotated)

# Generate figure
order4plot <- df.final_amr_mags %>%
  mutate(ratio = log10((num_amr+0.1)/(gene_count))) %>%
  group_by(Family) %>%
  dplyr::summarise(order = median(ratio)) %>% 
  arrange(order)

(ed_5k <- ggplot(df.final_amr_mags %>% mutate(Family = factor(Family, levels=order4plot$Family)), 
       aes(x=Family, y=log10((num_amr+0.1)/(gene_count)))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(alpha=0.4) +
  labs(x = "", y = expression('log'[10]*'(AMR/bacterial genes)')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.3)))

ggsave(paste0(fig_dir_ed_subpanels, 'ED_Fig_5k_mag_amr.pdf'), plot = ed_5k, width=6, height =4)


a_b_d <- plot_grid(ed_fig_5a, ed_fig_5b, ed_5d, labels = c("a", "b", "d"), label_size = 16, ncol = 3, rel_widths = c(1, 1.8, 1))
c <- plot_grid(ed_5c, labels = c("c"), label_size = 16, ncol = 1)
e_f_g <- plot_grid(ed_5e, plot_grid(ed_5f, ed_5g, labels = c("f", "g"), label_size = 16, nrow = 2), labels = c("e", ""), label_size = 16, ncol = 2, rel_widths = c(1.5, 1))
h_i <- plot_grid(ed_5h, ed_5i, labels = c("h", "i"), label_size = 16, ncol = 2, rel_widths = c(1.8, 1))
j_k <- plot_grid(ed_5j, ed_5k, labels = c("j", "k"), label_size = 16, ncol = 2, rel_widths = c(1, 2))

plot_grid(a_b_d, c, e_f_g, h_i, j_k,
          rel_heights = c(1, 1.2, 1.8, 1, 1), 
          nrow = 5)

ggsave(paste0(fig_dir, "ED_Figure_5.pdf"), width = 18, height = 25)
