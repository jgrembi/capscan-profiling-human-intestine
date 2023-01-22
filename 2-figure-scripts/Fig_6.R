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

# Only look at data of samples from Sets 2,3,4,5 that meet 16S cut-offs & bile acid cut-offs
ps.bs.raw <- readRDS(phyloseq_bilesalt) %>% subset_samples(., drop_16s == F)
ps.bs<-subset_samples(ps.bs.raw, Set %in% c('2','3','4','5','Stool'))
ps.bs@sam_data$Type <- gsub("Capsule", "Device", ps.bs@sam_data$Type)


#######################################
# Figure 6a - TCA Abundance through gastrointestinal tract
#######################################

df.ratios <- data.frame(sample_data(ps.bs)) %>%
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
    labs(y=expression('log'[10]*' TCA concentration'), x='',color='') +
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

# Generate list of significantly negatively correlated ASVs in Device samples
ddff <- ps.bs %>%
  subset_samples(., !Type %in% c('Stool')) %>% # right now, just looking at device samples
  filter_taxa(., function(x) sum(x > 3) > (0.05*length(x)), TRUE) %>%
  transform_sample_counts(., function(x) {log2(x+1)}) %>%
  psmelt()

# Select bile acid name
bile_acid<-'Taurocholic.acid'
# Remove zeros
df_b <- ddff %>% filter(Abundance > 0)
df_b$bile_acid<-df_b[,bile_acid]
# split the data by ASV
B <- split(df_b, df_b$OTU)
# Calculate the correlation in all data.frames using lapply 
M <- lapply(B, function(x) cor.test(x$Abundance, log10(x$bile_acid), method='pearson')) 
# Get ASVs that are significant
pvals<-lapply(M, function(x) x$p.value)
pvals<-as.data.frame(pvals)
pvals.corrected<-p.adjust(pvals,method="hochberg") %>% as.data.frame() %>% filter(. < 0.01) #this is where you set the pvalue
sig<-rownames(pvals.corrected)
# Slope 
r<-lapply(M, function(x) x$estimate)       
filtered<-data.frame(r) %>% t() %>% as.data.frame() %>% filter(cor < 0) # filter for negative correlations
filtered$OTU<-rownames(filtered)
sig.r <- filtered %>% filter(OTU %in% sig)
# Filter dataframe to only include significant ASVs              
df_b.filt <- df_b %>%
  filter(OTU %in% sig.r$OTU) %>%
  select(OTU, Phylum, Family, Genus, Species) %>%
  unique() %>%
  mutate(BA = bile_acid) %>%
  left_join(.,sig.r, by='OTU')
df_b.filt


# Make a dataframe to plot correlations in Device and stool
df2plot_b <- psmelt(transform_sample_counts(ps.bs, function(x) {log2(x+1)})) %>%
  filter(OTU %in% df_b.filt$OTU) %>% #Sig negatively correlated ASVs
  mutate(Location = ifelse(Type %in% c('Device 1','Device 2','Device 3','Device 4'),'Device','Stool')) %>%
  select(OTU, Main, Type, Location, Abundance, Subject, Genus, Species, Taurocholic.acid, reads_16s) %>%
  mutate(title_name = paste0(Genus, ' ', Species)) %>%
  filter(Abundance > 0)

#formula <- 'y ~ x'

(tca_cor<-ggplot(df2plot_b, aes(x=Abundance, y=log10(Taurocholic.acid), color=Location))+
    geom_point(size=3,alpha=0.1)+
    geom_smooth(method='lm', alpha=0.2, aes(fill = Location))+
    stat_cor(method = 'pearson',p.accuracy = 0.001, r.accuracy = 0.01,
             label.x.npc = 0.5,label.y.npc = "top")+
    #ylim(0,15)+
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
    labs(x=expression('log'[2]*' Abundance'), 
         y=expression('log'[10]*' TCA concentration'),
         color = "")+
    facet_wrap(~title_name + Location)+#, scale='free_y')+
    theme(legend.position = "bottom",
          legend.margin = margin(),
          strip.text = element_text(size = 10, face = "bold"), 
          legend.spacing.x = unit(1.0, 'pt'),
          legend.text=element_markdown(size=10), 
          axis.text = element_text(size = 8),
          plot.margin = margin(5,5,5,5)))+
  guides(color = guide_legend(nrow=1, 
                              override.aes = list(label="", face=NA, shape = NA, fill=NA,linetype=NA))) ## This line removes the point on the legend for the color

ggsave(filename = paste0(fig_dir_main_subpanels, "Fig_6b_subpanel_bilophila_TCA_correlations.pdf"), tca_cor)

#######################################
# Figure 6c & Figure 6d - % of microbially-conjugated bile acids & concentration of microbially conjugated bile acids
#######################################

## Microbially conjugated bile acids
df2plot_c <- df.ratios %>%
  mutate(pct_microConj = log10(microConj)/log10(primarySum+secondarySum+priConjSum+secConjSum+microConj)*100)

# my_comparisons <- list(c('Device 1','Device 4'),c('Device 4','Stool'),c('Device 1','Stool'))
stat_test_microConj_pct <- df2plot_c %>%
  wilcox_test(pct_microConj ~ Type, p.adjust.method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "Type") %>%
  filter((group1 == "Device 1" & group2 %in% c("Device 4", "Stool")) | (group1 == "Device 4" & group2 == "Stool"))

## This plot shows the ratio of microbially-conjugated bile acids/all bile acids
(pct_micro <- ggplot(df2plot_c, aes(y=pct_microConj, x=Type))+
    geom_boxplot(aes(fill=Type), outlier.shape = NA, alpha = 0.9)+
    #geom_violin(alpha = 0.9) +
    geom_beeswarm(alpha=0.4)+
    scale_fill_manual(values=CapTypeAndStoolColors) +
    stat_pvalue_manual(stat_test_microConj_pct, 
                       label = "p.adj.signif",
                       tip.length = 0.01,
                       y.position = c(80, 83, 80), 
                       bracket.shorten = 0.05) + 
    labs(y=expression('% of microbially conjguated bile acids', x='')) +
    theme(legend.position='none',axis.text = element_text(size = 10),axis.title.x=element_blank()))

ggsave(filename = paste0(fig_dir_main_subpanels, "Fig_6c_subpanel_pct_microbially_conj_bile_acids.pdf"), pct_micro)


## This plot shows the concentration of microbially-conjugated bile acids

stat_test_microConj_conc <- df2plot_c %>%
  wilcox_test(microConj ~ Type, p.adjust.method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "Type") %>%
  filter((group1 == "Device 1" & group2 %in% c("Device 4", "Stool")) | (group1 == "Device 4" & group2 == "Stool"))

(abund_micro <- ggplot(df2plot_c, aes(y=log10(microConj), x=Type))+
    geom_boxplot(aes(fill=Type), outlier.shape = NA, alpha = 0.9)+
    #geom_violin(alpha = 0.9) +
    geom_beeswarm(alpha=0.4)+
    scale_fill_manual(values=CapTypeAndStoolColors) +
    stat_pvalue_manual(stat_test_microConj_conc, 
                       label = "p.adj.signif",
                       tip.length = 0.01,
                       y.position = c(4.8, 5, 4.8), 
                       bracket.shorten = 0.05) +
    labs(y=expression('log'[10]*'(concentration of microbially conjguated bile acids)', x='')) +
    theme(legend.position='none',axis.text = element_text(size = 10),axis.title.x=element_blank()))

ggsave(filename = paste0(fig_dir_main_subpanels, "Fig_6d_subpanel_abundance_microbially_conj_bile_acids.pdf"), abund_micro)



#######################################
# Figure 6e - Correlations between all bile acids in devices
#######################################          

# Create dataframe to plot
bs_capsule <- data.frame(sample_data(ps.bs)) %>%
  filter(Type %in% c('Device 1','Device 2','Device 3','Device 4')) %>%
  select("Tauro.a.Muricholic.acid":'Tyr.Dihydroxylated.BA') %>%
  #rename_all(~ gsub(".acid", " acid", .)) %>%
  na.omit() %>%
  select_if(colSums(.) > 0) %>%
  mutate_all(~ log10(. + 0.000001)) %>%
  as.matrix 

res_capsule <- rcorr(bs_capsule)

# Generate figure subpanel e
pdf(height = 10, width = 12, file = paste0(fig_dir_main_subpanels,'Fig_6e_subpanel_capsule_bileacid_correlations.pdf'))
corrplot(res_capsule$r, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

dev.off()


## Abbreviated/summarised correlations (insert in main figure)
bs_capsule_summary <- data.frame(sample_data(ps.bs)) %>%
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
bs_stool <- data.frame(sample_data(ps.bs)) %>%
  filter(Type %in% 'Stool') %>%
  select("Tauro.a.Muricholic.acid":'Tyr.Dihydroxylated.BA') %>%
  #rename_all(~ gsub(".acid", " acid", .)) %>%
  na.omit() %>%
  select_if(colSums(.) > 0) %>%
  mutate_all(~ log10(. + 0.000001)) %>%
  as.matrix

res_stool <- rcorr(bs_stool)

# Generate figure
pdf(height = 10, width = 12, file = paste0(fig_dir_main_subpanels,'Fig_6f_subpanel_stool_bileacid_correlations.pdf'))
corrplot(res_stool$r, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

dev.off()


## Abbreviated/summarised correlations (insert in main figure)
bs_stool_summary <- data.frame(sample_data(ps.bs)) %>%
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

df2plot_ghi <- data.frame(sample_data(ps.bs))  %>%
  mutate(Location = ifelse(Type %in% 'Stool','Stool','Device'),
         Ursodeoxycholic.acid = log10(Ursodeoxycholic.acid +1),
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
