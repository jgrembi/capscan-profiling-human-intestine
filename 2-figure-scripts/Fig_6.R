#######################################
# Capsule Study - Gen1 
#
# Generate Figure 6 for manuscript
#
# Authors: Becca Culver & Jess Grembi
#######################################

rm(list=ls())
# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# Only look at data of samples from Sets 2,3,4,5 that meet bile acid cut-offs
bs <- read.csv(paste0(data_dir, "bile_salts_wAAconjugates.csv")) %>%
  dplyr::rename(Super = Sample.ID) 

# Number of samples with bile acid data
table(bs$SampleType)

df_samples <- readRDS(paste0(data_dir, "sample_data.RDS")) %>%
  mutate(Type = gsub("Capsule", "Device", Type), 
         Super = ifelse(Type == "Stool", Main, Super)) %>%
  filter(!duplicate_16s, !duplicate_meta, location != "Saliva") %>%
  select(DeviceID:Super, location, drop_16s)

df_bs <- bs %>%
  select(-(Subject.Number:Capsule.Type)) %>%
  mutate(Super = as.character(Super)) %>%
  left_join(df_samples, by = "Super") %>%
  unique() %>%
  filter(Set %in% c('2','3','4','5','Stool')) #%>% # Remove Set 1 & additional microbiota reproducibility samples

table(df_bs$location)
table(df_bs$Type)
#######################################
# Figure 6a - TCA Abundance through gastrointestinal tract
#######################################

df.ratios <- df_bs %>%
  mutate(primarySum = Cholic.acid + Chenodeoxycholic.acid) %>%
  mutate(priConjSum = Glycocholic.acid+Taurocholic.acid+Glycochenodeoxycholic.acid+Taurochenodeoxycholic.acid+Tauro.a.Muricholic.acid) %>%
  mutate(secondarySum = Deoxycholic.acid+Lithocholic.acid+Ursodeoxycholic.acid) %>%
  mutate(secConjSum = Glycodeoxycholic.acid+Taurodeoxycholic.acid+Glycolithocholic.acid+Taurolithocholic.acid+Glycoursodeoxycholic.acid+Tauroursodeoxycholic.acid+Glycohyodeoxycholic.acid) %>%
  mutate(microConj = Phenylalanocholic.acid+Leucholic.acid+Tyrosocholic.acid+Ala.TriHydroxylated.BA+Arg.Dihydroxylated.BA+Arg.TriHydroxylated.BA+Asn.Dihydroxylated.BA.1+Asn.Dihydroxylated.BA.2+Asn.TriHydroxylated.BA+Cys.Dihydroxylated.BA+Cys.TriHydroxylated.BA+Gln.Dihydroxylated.BA
         +Gln.TriHydroxylated.BA+Glu.Dihydroxylated.BA+Glu.TriHydroxylated.BA+His.TriHydroxylated.BA+Lys.Dihydroxylated.BA+Lys.TriHydroxylated.BA+Met.TriHydroxylated.BA+Phe.Dihydroxylated.BA+Ser.Dihydroxylated.BA+Ser.TriHydroxylated.BA+Trp.TriHydroxylated.BA+Tyr.Dihydroxylated.BA) %>%
  mutate(Location = ifelse(Type %in% 'Stool','Stool','Device'))


stat_test_tca <- df.ratios %>%
  wilcox_test(Taurocholic.acid ~ Type, p.adjust.method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "Type") %>%
  filter((group1 == "Device 1" & group2 %in% c("Device 4", "Stool")) | (group1 == "Device 4" & group2 == "Stool"))

(tca <- ggplot(df.ratios, aes(y=log10(Taurocholic.acid), x=Type))+
    geom_boxplot(aes(fill=Type), outlier.shape = NA, alpha = 0.9)+
    #geom_violin(alpha = 0.9) +
    geom_beeswarm(alpha=0.4)+
    scale_fill_manual(values=CapTypeAndStoolColors)+
    labs(y=expression('log'[10]*'(TCA concentration)'), x='',color='') +
    stat_pvalue_manual(stat_test_tca, 
                       label = "p.adj.signif",
                       tip.length = 0.01,
                       y.position = c(6.4, 6.7, 6.4), 
                       bracket.shorten = 0.05) +    
    theme(legend.position='none',axis.text = element_text(size = 10),axis.title.x=element_blank()))

ggsave(paste0(fig_dir_main_subpanels,'Fig_6a_subpanel_tca_concentration.pdf'), tca, width=6,height=4)


#######################################
# Figure 6b - ASV abundance correlations with TCA
#######################################
# Only look at data of samples from Sets 2,3,4,5 that meet 16S cut-offs & bile acid cut-offs
ps.bs.raw <- readRDS(phyloseq_bilesalt) %>% subset_samples(., drop_16s == F)
ps.bs<-subset_samples(ps.bs.raw, Set %in% c('2','3','4','5','Stool'))
ps.bs@sam_data$Type <- gsub("Capsule", "Device", ps.bs@sam_data$Type)

table(ps.bs@sam_data$location)
# Generate list of significantly negatively correlated ASVs in Device samples
ddff <- ps.bs %>%
  subset_samples(., !Type %in% c('Stool')) %>% # right now, just looking at device samples
  filter_taxa(., function(x) sum(x > 3) > (0.05*length(x)), TRUE) %>%
  transform_sample_counts(., function(x) {log2(x+1)}) %>%
  psmelt()

table(ddff %>% select(Sample, Type, location) %>% unique %>% pull(location))
table(ddff %>% select(Sample, Type, location) %>% unique %>% pull(Type))

# Use function to obtain ASVs correlated with Taurocholic.acid
df_b.filt <- asv_corr(ddff, location_test = "Capsule", cor_variable = 'Taurocholic.acid', include_0_abundance = F) %>%
  filter(adj_pval < 0.01) %>%
  filter(rho < 0)

ddff_stool <- ps.bs %>%
  subset_samples(., Type %in% c('Stool')) %>% # right now, just looking at device samples
  filter_taxa(., function(x) sum(x > 3) > (0.05*length(x)), TRUE) %>%
  transform_sample_counts(., function(x) {log2(x+1)}) %>%
  psmelt()
df_b.stool <- asv_corr(ddff_stool, location_test = "Stool", cor_variable = 'Taurocholic.acid', include_0_abundance = F) %>%
  filter(ASV %in% df_b.filt$ASV)

# Make a dataframe to plot correlations in Device and stool
df2plot_b <- ps.bs %>%
  prune_taxa(df_b.filt$ASV, .) %>% #Sig negatively correlated ASVs
  transform_sample_counts(., function(x) {log2(x+1)}) %>%
  psmelt() %>%
  mutate(Location = ifelse(Type %in% c('Device 1','Device 2','Device 3','Device 4'),'Devices','Stool')) %>%
  select(ASV, Main, Type, Location, Abundance, Subject, Genus, Species, Taurocholic.acid, reads_16s) %>%
  mutate(title_name = paste0(Genus, ' ', Species)) %>%
  filter(Abundance > 0) 
df_stats <- bind_rows(df_b.filt, df_b.stool) %>%
  left_join(df2plot_b %>% 
              select(ASV, title_name) %>%
              unique, by = "ASV") 
df_stats
(tca_cor<-ggplot(df2plot_b, aes(x=Abundance, y=log10(Taurocholic.acid), color=Location))+
    geom_point(size=3,alpha=0.1)+
    geom_smooth(method='lm', alpha=0.2, aes(fill = Location))+
    stat_cor(method = 'spearman',
             aes(label = after_stat(r.label)),
             cor.coef.name = "rho",
             r.accuracy = 0.01,
             label.x.npc = "center",
             label.y = 5) + 
    scale_color_manual(labels = paste("<span style='color:",
                                      CapAndStoolColors,"'>",
                                      levels(df2plot_b$Location),
                                      "</span>"),
                       values = CapAndStoolColors, 
                       guide = "none") +  
    scale_fill_manual(labels = paste("<span style='color:",
                                     CapAndStoolColors,"'>",
                                     levels(df2plot_b$Location),
                                     "</span>"),
                      values = CapAndStoolColors, 
                      guide = "none") + 
    labs(x=expression('log'[2]*'(ASV count)'), 
         y=expression('log'[10]*'(TCA concentration)'))+
    facet_wrap(~title_name + Location) +
    theme(strip.text = element_text(size = 10, face = "bold"), 
          axis.text = element_text(size = 8),
          plot.margin = margin(5,5,5,5))) 
  
ggsave(filename = paste0(fig_dir_main_subpanels, "Fig_6b_subpanel_ASV_TCA_correlations.pdf"), tca_cor, height = 10, width = 12)

#######################################
# Figure 6c & Figure 6d - % of microbially-conjugated bile acids & concentration of microbially conjugated bile acids
#######################################

## Microbially conjugated bile acids
df2plot_c <- df.ratios %>%
  mutate(pct_microConj = log10(microConj)/log10(primarySum+secondarySum+priConjSum+secConjSum+microConj)*100)

## This plot shows the concentration of microbially-conjugated bile acids
stat_test_microConj_conc <- df2plot_c %>%
  wilcox_test(microConj ~ Type, p.adjust.method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "Type") %>%
  filter((group1 == "Device 1" & group2 %in% c("Device 4", "Stool")) | (group1 == "Device 4" & group2 == "Stool"))

(abund_micro <- ggplot(df2plot_c, aes(y=log10(microConj), x=Type))+
    geom_boxplot(aes(fill=Type), outlier.shape = NA, alpha = 0.9)+
    geom_beeswarm(alpha=0.4)+
    scale_fill_manual(values=CapTypeAndStoolColors) +
    stat_pvalue_manual(stat_test_microConj_conc, 
                       label = "p.adj.signif",
                       tip.length = 0.01,
                       y.position = c(4.8, 5, 4.8), 
                       bracket.shorten = 0.05) +
    labs(y=expression('log'[10]*'(concentration of microbially conjguated bile acids)', x='')) +
    theme(legend.position='none',axis.text = element_text(size = 10),axis.title.x=element_blank()))

ggsave(filename = paste0(fig_dir_main_subpanels, "Fig_6c_subpanel_abundance_microbially_conj_bile_acids.pdf"), abund_micro)

## This plot shows the ratio of microbially-conjugated bile acids/all bile acids
stat_test_microConj_pct <- df2plot_c %>%
  wilcox_test(pct_microConj ~ Type, p.adjust.method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "Type") %>%
  filter((group1 == "Device 1" & group2 %in% c("Device 4", "Stool")) | (group1 == "Device 4" & group2 == "Stool"))

(pct_micro <- ggplot(df2plot_c, aes(y=pct_microConj, x=Type))+
    geom_boxplot(aes(fill=Type), 
                 outlier.shape = NA, 
                 alpha = 0.9)+
    geom_beeswarm(alpha=0.4)+
    scale_fill_manual(values=CapTypeAndStoolColors) +
    stat_pvalue_manual(stat_test_microConj_pct, 
                       label = "p.adj.signif",
                       tip.length = 0.01,
                       y.position = c(80, 83, 80), 
                       bracket.shorten = 0.05) + 
    labs(y=expression('% of microbially conjguated bile acids', x='')) +
    theme(legend.position='none',
          axis.text = element_text(size = 10),
          axis.title.x=element_blank()))

ggsave(filename = paste0(fig_dir_main_subpanels, "Fig_6d_subpanel_pct_microbially_conj_bile_acids.pdf"), pct_micro)


#######################################
# Figure 6e - Correlations between all bile acids in devices
#######################################          

# Create dataframe to plot
bs_capsule <- df_bs %>%
  filter(Type %in% c('Device 1','Device 2','Device 3','Device 4')) %>%
  select("Tauro.a.Muricholic.acid":'Tyr.Dihydroxylated.BA') %>%
  rename_all(~ gsub(".acid", " acid", .)) %>%
  rename_all(~ gsub(".BA", " BA", .)) %>%
  rename_all(~ gsub("BA.", "BA ", .)) %>%
  rename_all(~ gsub("\\.", "-", .)) %>%
  na.omit() %>%
  select_if(colSums(.) > 0) %>%
  mutate_all(~ log10(. + 0.000001)) %>%
  as.matrix 

res_capsule <- rcorr(bs_capsule)

# Generate figure subpanel e
pdf(height = 10, width = 14, file = paste0(fig_dir_main_subpanels,'Fig_6e_subpanel_capsule_bileacid_correlations.pdf'))
corrplot(res_capsule$r, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

dev.off()


## Abbreviated/summarised correlations (insert in main figure)
bs_capsule_summary <- df_bs %>%
  filter(Type %in% c('Device 1','Device 2','Device 3','Device 4')) %>%
  mutate(primarySum = Cholic.acid + Chenodeoxycholic.acid) %>%
  mutate(priConjSum = Glycocholic.acid+Taurocholic.acid+Glycochenodeoxycholic.acid+Taurochenodeoxycholic.acid+Tauro.a.Muricholic.acid) %>%
  mutate(secondarySum = Deoxycholic.acid+Lithocholic.acid+Ursodeoxycholic.acid) %>%
  mutate(secConjSum = Glycodeoxycholic.acid+Taurodeoxycholic.acid+Glycolithocholic.acid+Taurolithocholic.acid+Glycoursodeoxycholic.acid+Tauroursodeoxycholic.acid+Glycohyodeoxycholic.acid) %>%
  mutate(microConj = Phenylalanocholic.acid+Leucholic.acid+Tyrosocholic.acid+Ala.TriHydroxylated.BA+Arg.Dihydroxylated.BA+Arg.TriHydroxylated.BA+Asn.Dihydroxylated.BA.1+Asn.Dihydroxylated.BA.2+Asn.TriHydroxylated.BA+Cys.Dihydroxylated.BA+Cys.TriHydroxylated.BA+Gln.Dihydroxylated.BA
         +Gln.TriHydroxylated.BA+Glu.Dihydroxylated.BA+Glu.TriHydroxylated.BA+His.TriHydroxylated.BA+Lys.Dihydroxylated.BA+Lys.TriHydroxylated.BA+Met.TriHydroxylated.BA+Phe.Dihydroxylated.BA+Ser.Dihydroxylated.BA+Ser.TriHydroxylated.BA+Trp.TriHydroxylated.BA+Tyr.Dihydroxylated.BA) %>%
  select(primarySum, priConjSum, secondarySum, secConjSum, microConj) %>%
  na.omit() %>%
  select_if(colSums(.) > 0) %>%
  mutate_all(~ log10(. + 0.000001)) %>%
  as.matrix

res_capsule_summary <- rcorr(bs_capsule_summary)

# Generate the pdf of the figure          
pdf(height = 8, width = 8, file = paste0(fig_dir_main_subpanels,'Fig_6e_subpanel_insert_capsule_ba_abbr_correlations.pdf'))
corrplot(res_capsule_summary$r, type = "upper", order='alphabet',#order = "hclust", 
         tl.col = "black", tl.srt = 45)

dev.off()

#######################################
# Figure 6f - Correlations between all bile acids in stool
#######################################

# Create dataframe to plot
bs_stool <- df_bs %>%
  filter(Type %in% 'Stool') %>%
  select("Tauro.a.Muricholic.acid":'Tyr.Dihydroxylated.BA') %>%
  rename_all(~ gsub(".acid", " acid", .)) %>%
  rename_all(~ gsub(".BA", " BA", .)) %>%
  rename_all(~ gsub("BA.", "BA ", .)) %>%
  rename_all(~ gsub("\\.", "-", .)) %>%
  na.omit() %>%
  select_if(colSums(.) > 0) %>%
  mutate_all(~ log10(. + 0.000001)) %>%
  as.matrix

res_stool <- rcorr(bs_stool)

# Generate figure
pdf(height = 10, width = 14, file = paste0(fig_dir_main_subpanels,'Fig_6f_subpanel_stool_bileacid_correlations.pdf'))
corrplot(res_stool$r, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

dev.off()


## Abbreviated/summarised correlations (insert in main figure)
bs_stool_summary <- df_bs %>%
  filter(Type %in% 'Stool') %>%
  mutate(primarySum = Cholic.acid + Chenodeoxycholic.acid) %>%
  mutate(priConjSum = Glycocholic.acid+Taurocholic.acid+Glycochenodeoxycholic.acid+Taurochenodeoxycholic.acid+Tauro.a.Muricholic.acid) %>%
  mutate(secondarySum = Deoxycholic.acid+Lithocholic.acid+Ursodeoxycholic.acid) %>%
  mutate(secConjSum = Glycodeoxycholic.acid+Taurodeoxycholic.acid+Glycolithocholic.acid+Taurolithocholic.acid+Glycoursodeoxycholic.acid+Tauroursodeoxycholic.acid+Glycohyodeoxycholic.acid) %>%
  mutate(microConj = Phenylalanocholic.acid+Leucholic.acid+Tyrosocholic.acid+Ala.TriHydroxylated.BA+Arg.Dihydroxylated.BA+Arg.TriHydroxylated.BA+Asn.Dihydroxylated.BA.1+Asn.Dihydroxylated.BA.2+Asn.TriHydroxylated.BA+Cys.Dihydroxylated.BA+Cys.TriHydroxylated.BA+Gln.Dihydroxylated.BA
         +Gln.TriHydroxylated.BA+Glu.Dihydroxylated.BA+Glu.TriHydroxylated.BA+His.TriHydroxylated.BA+Lys.Dihydroxylated.BA+Lys.TriHydroxylated.BA+Met.TriHydroxylated.BA+Phe.Dihydroxylated.BA+Ser.Dihydroxylated.BA+Ser.TriHydroxylated.BA+Trp.TriHydroxylated.BA+Tyr.Dihydroxylated.BA) %>%
  select(primarySum, priConjSum, secondarySum, secConjSum, microConj) %>%
  na.omit() %>%
  select_if(colSums(.) > 0) %>%
  mutate_all(~ log10(. + 0.000001)) %>%
  as.matrix


res_stool_summary <- rcorr(bs_stool_summary)

# Generate the pdf of the figure
pdf(height = 8, width = 8, file = paste0(fig_dir_main_subpanels,'Fig_6f_subpanel_insert_stool_ba_abbr_correlations.pdf'))
corrplot(res_stool_summary$r, type = "upper", order='alphabet',#order = "hclust", 
         tl.col = "black", tl.srt = 45)

dev.off()


#######################################
# Figure 6g - Gln-Trihydroxy concentration
#######################################          

df2plot_ghi <- df_bs  %>%
  mutate(Ursodeoxycholic.acid = log10(Ursodeoxycholic.acid +1),
         Ser.TriHydroxylated.BA = log10(Ser.TriHydroxylated.BA +1)) 

stat_test_g <- df2plot_ghi %>%
  wilcox_test(Gln.TriHydroxylated.BA ~ Type, p.adjust.method = "bonferroni") %>%
  add_significance() %>%
  filter(!p.adj.signif %in% 'ns') %>%
  add_xy_position(x = "Type") %>%
  mutate(y.position = y.position - 500) 

(gln <- ggplot(df2plot_ghi, aes(y=Gln.TriHydroxylated.BA, x=Type))+
    geom_boxplot(aes(fill=Type), outlier.shape = NA, alpha = 0.9)+
    geom_beeswarm(alpha=0.4, cex=0.5)+ 
    scale_fill_manual(values=CapTypeAndStoolColors)+
    labs(y='Gln-Trihydroxy concentration', x='',color='') +
    stat_pvalue_manual(stat_test_g, 
                       label = "p.adj.signif",
                       tip.length = 0.01,
                       y.position = c(500, 600, 700, 800, 900, 1000),
                       bracket.shorten = 0.05) +    
    theme(legend.position='none',axis.text = element_text(size = 10),axis.title.x=element_blank()))

ggsave(paste0(fig_dir_main_subpanels, 'Fig_6g_subpanel_glntrihydroxy_conc.pdf'), gln, width = 4, height = 3)  

#######################################
# Figure 6h - Ursodeoxycholic acid concentration
#######################################    

stat_test_h <- df2plot_ghi %>%
  wilcox_test(Ursodeoxycholic.acid ~ Type, p.adjust.method = "bonferroni") %>%
  add_significance() %>%
  filter(!p.adj.signif %in% 'ns') %>%
  add_xy_position(x = "Type") 

(urso <- ggplot(df2plot_ghi, aes(y=Ursodeoxycholic.acid, x=Type))+
    geom_boxplot(aes(fill=Type), outlier.shape = NA, alpha = 0.9)+
    geom_beeswarm(alpha=0.4)+
    scale_fill_manual(values=CapTypeAndStoolColors)+
    labs(y=expression('log'[10]*'(Ursodeoxycholic acid concentration)'), x='',color='') +
    stat_pvalue_manual(stat_test_h, 
                       label = "p.adj.signif",
                       tip.length = 0.01,
                       y.position = c(5.2,5.7,6.2,6.7), 
                       bracket.shorten = 0.05) +    
    theme(legend.position='none',axis.text = element_text(size = 10),axis.title.x=element_blank()))

ggsave(paste0(fig_dir_main_subpanels, 'Fig_6h_subpanel_urso_conc.pdf'), urso, width = 4, height = 3)


#######################################
# Figure 6i - Ser-Trihydroxy concentration
#######################################          

stat_test_i <- df2plot_ghi %>%
  wilcox_test(Ser.TriHydroxylated.BA ~ Type, p.adjust.method = "bonferroni") %>%
  add_significance() %>%
  filter(!p.adj.signif %in% 'ns') %>%
  add_xy_position(x = "Type") 

(ser <- ggplot(df2plot_ghi, aes(y=Ser.TriHydroxylated.BA, x=Type))+
    geom_boxplot(aes(fill=Type), outlier.shape = NA, alpha = 0.9)+
    geom_beeswarm(alpha=0.4, cex=0.5)+
    scale_fill_manual(values=CapTypeAndStoolColors)+
    labs(y='log10 Ser.TriHydroxylated.BA', x='',color='') +
    stat_pvalue_manual(stat_test_i, 
                       label = "p.adj.signif",
                       tip.length = 0.01,
                       y.position = c(4,4.5,5,5.5), 
                       bracket.shorten = 0.05) +    
    theme(legend.position='none',axis.text = element_text(size = 10),axis.title.x=element_blank()))

ggsave(paste0(fig_dir_main_subpanels, 'Fig_6i_subpanel_ser_conc.pdf'), ser, width = 4, height = 3)         
