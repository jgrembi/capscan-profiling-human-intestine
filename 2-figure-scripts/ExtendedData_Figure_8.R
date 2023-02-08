#######################################
# Capsule Study - Gen1 
#
# Extended Data Figure 8
#
# Author: Jess Grembi and Rebecca Culver
#######################################

rm(list=ls())
# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

#######################################
# ED Fig. 8a,b,c
#######################################

### sample data
df_samples <- readRDS(sample_data) 


## data from HMM search for bsh genes
bsh_hmm <- read.csv(paste0(data_dir, "0_compiled_bsh_hits.csv"), header = F, 
                    col.names = c("hit_id", "targ_id", "target_len", "hit_score", "hit_e-value", "hit_bias", "hsp_len", "hsp_score", "hsp_e-value_cond", "hsp_e-value_ind", "hsp_bias", "filepath")) %>%
  mutate(meta_samplename = gsub("workflow/out/bsh_hmmer/", "", gsub("_contigs.domtblout", "", filepath)))

bsh_hmm_filt <- bsh_hmm %>%
  filter(hsp_e.value_ind < 10e-10, hit_score > 350) %>%
  separate(targ_id, c("k", "contig", "orf"), sep = "_") %>%
  mutate(bsh = "bsh positive",
         contigName = paste0(k, "_", contig)) 


n_genes <- read.table(paste0(data_dir, "contigs_n_genes_summary.txt"), col.names = c("meta_samplename", "n_genes")) %>%
  mutate(meta_samplename = gsub("workflow/out/prodigal/contigs/", "", gsub("_contigs_out.fna", "", meta_samplename))) %>%
  left_join(df_samples, by = "meta_samplename") 

# median(n_genes$n_genes[n_genes$location %in% c("Stool", "Capsule")], na.rm = T)

df_summary <- bsh_hmm_filt %>% 
  group_by(meta_samplename) %>%
  summarise(n = n()) %>%
  left_join(n_genes %>%
              select(-n), by = "meta_samplename") %>%
  mutate(n_norm = n/n_genes) %>%
  mutate(Subject = fct_reorder(Subject, n_norm, .fun = 'median', na.rm = T)) %>%
  filter(!drop_meta, Set %in% c("2", "3", '4', "5", "Stool"))

table(df_summary$location)

(bsh_uniqueBySubj_plot <- ggplot(df_summary,
                                 aes(x = Subject, y = n_norm*100, color = location)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point(position = position_jitterdodge(0.2), alpha = 0.4) +
    # geom_beeswarm(color = "black") +
    scale_y_log10() +
    # ggsignif::geom_signif(comparisons = list(c("Small intestine", "Stool"))) +
    scale_color_manual(values = CapAndStoolColors) +
    labs(y = "% of genes identified as BSH", x = "Subject", color = "Location") +
    theme(plot.margin = margin(t = fig_padding, l = fig_padding, r = fig_padding, b = fig_padding, unit = "pt"),
          legend.position = "none"))

ggsave(filename = paste0(fig_dir_ed_subpanels, "ED_Fig_8a_summary_pct_bsh_genes_by_subj.pdf"), plot = bsh_uniqueBySubj_plot)

df_summary %>%
  group_by(location) %>%
  summarise(median = median(n_norm))

wilcox.test(x = df_summary$n_norm[df_summary$location == "Capsule"], 
            y = df_summary$n_norm[df_summary$location == "Stool"],
            alternative = "two.sided",
            conf.int = T)

## Contig depth data
df_contigs <- read.table(paste0(data_dir, "0_contig_depth_summary.txt"), header = T) %>%
  dplyr::rename(meta_samplename = A1_1)


contig_ranks <- df_contigs %>%
  group_by(meta_samplename) %>%
  mutate(rank = rank(totalAvgDepth), 
         n = n(), 
         rank = rank/n,
         contig_sample = paste0(contigName, "_", meta_samplename)) %>%
  ungroup() %>%
  left_join(bsh_hmm_filt, by = c("meta_samplename", "contigName")) %>%
  mutate(bsh = ifelse(is.na(bsh), "bsh negative", bsh)) %>%
  left_join(df_samples %>%
              select(Subject:Main, location, meta_samplename, drop_meta) %>%
              unique, by = "meta_samplename") %>%
  mutate(Type = gsub("Capsule", "Device", Type),
         location = gsub("Capsule", "Devices", location))

table(contig_ranks %>% filter(!drop_meta, Set %in% c("2", "3", "4", "5", "Stool")) %>% select(meta_samplename, Type, location) %>% unique() %>% pull(location))
table(contig_ranks$Type)

(bsh_coverage_hist <- ggplot(contig_ranks %>% filter(location %in% c("Stool", "Devices"), bsh == "bsh positive"), 
                             aes(x = rank, fill = location)) + 
    geom_histogram(binwidth = .005) + facet_wrap(~location, nrow = 2) + 
    scale_fill_manual(values = CapAndStoolColors) + 
    labs(x = "Rank coverage of BSH genes", y = "Count") + 
    theme(legend.position = "none",
          plot.margin = margin(t = fig_padding, l = fig_padding, r = fig_padding, b = fig_padding, unit = "pt")))

ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_8b_bsh_rank_distribution.pdf"), plot = bsh_coverage_hist)


(bsh_rank_boxplot <- ggplot(contig_ranks %>%
                              filter(bsh == "bsh positive") %>%
                              mutate(Type = factor(Type, levels = c("Stool", "Device 4", "Device 3", "Device 2", "Device 1"))) ,
                            aes(x = Type, y = rank, color = Type)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.4, alpha = 0.4) +
    scale_color_manual(values = rev(CapTypeAndStoolColors), guide = guide_legend(reverse = TRUE)) +
    labs(x = "", y = "Rank coverage of BSH genes", color = "") + 
    theme(legend.position = "none", 
          legend.margin=margin(),
          legend.box = "vertical",
          plot.margin = margin(t = fig_padding, l = fig_padding, r = fig_padding, b = fig_padding, unit = "pt")) + 
    coord_flip())

ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_8c_bsh_rank_coverage_by_type.pdf"), plot = bsh_rank_boxplot)

#######################################
# ED Fig. 8d - Primary bile acid concentrations
#######################################

ps.bs.raw <- readRDS(phyloseq_bilesalt) %>% subset_samples(., drop_16s == F)
ps.bs<-subset_samples(ps.bs.raw, Set %in% c('2','3','4','5','Stool'))
ps.bs@sam_data$Type <- gsub("Capsule", "Device", ps.bs@sam_data$Type)

df.ratios <- data.frame(sample_data(ps.bs)) %>%
    mutate(primarySum = Cholic.acid + Chenodeoxycholic.acid) %>%
    mutate(priConjSum = Glycocholic.acid+Taurocholic.acid+Glycochenodeoxycholic.acid+Taurochenodeoxycholic.acid+Tauro.a.Muricholic.acid) %>%
    mutate(secondarySum = Deoxycholic.acid+Lithocholic.acid+Ursodeoxycholic.acid) %>%
    mutate(secConjSum = Glycodeoxycholic.acid+Taurodeoxycholic.acid+Glycolithocholic.acid+Taurolithocholic.acid+Glycoursodeoxycholic.acid+Tauroursodeoxycholic.acid+Glycohyodeoxycholic.acid) %>%
    mutate(microConj = Phenylalanocholic.acid+Leucholic.acid+Tyrosocholic.acid+Ala.TriHydroxylated.BA+Arg.Dihydroxylated.BA+Arg.TriHydroxylated.BA+Asn.Dihydroxylated.BA.1+Asn.Dihydroxylated.BA.2+Asn.TriHydroxylated.BA+Cys.Dihydroxylated.BA+Cys.TriHydroxylated.BA+Gln.Dihydroxylated.BA
           +Gln.TriHydroxylated.BA+Glu.Dihydroxylated.BA+Glu.TriHydroxylated.BA+His.TriHydroxylated.BA+Lys.Dihydroxylated.BA+Lys.TriHydroxylated.BA+Met.TriHydroxylated.BA+Phe.Dihydroxylated.BA+Ser.Dihydroxylated.BA+Ser.TriHydroxylated.BA+Trp.TriHydroxylated.BA+Tyr.Dihydroxylated.BA)


df2plot_d <- df.ratios %>%
    mutate(Location = ifelse(Type %in% 'Stool','Stool','Device'),
           perc_priBAs = (priConjSum+primarySum)/(primarySum+secondarySum+priConjSum+secConjSum+microConj)*100)

my_comparisons <- list(c('Device 1','Device 4'),c('Device 4','Stool'),c('Device 1','Stool'))

stat_test <- df2plot_d %>%
  wilcox_test(perc_priBAs ~ Type, p.adjust.method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "Type") %>%
  filter((group1 == "Device 1" & group2 %in% c("Device 4", "Stool")) | (group1 == "Device 4" & group2 == "Stool"))

(d <- ggplot(df2plot_d, aes(x=Type, y=perc_priBAs))+
    geom_boxplot(aes(fill=Type), alpha=0.9, outlier.shape = NA) +
    #geom_violin(alpha = 0.9) +
    geom_beeswarm(alpha=0.4)+
    scale_fill_manual(values=CapTypeAndStoolColors)+
    stat_pvalue_manual(stat_test, 
                       label = "p.adj.signif",
                       tip.length = 0.01,
                       y.position = c(100,110,120), 
                       bracket.shorten = 0.05) +
    labs(title='',y='Percent of 1° bile acids',x='')+
    theme(legend.position='none'))

ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_8d_pct_primary_bile_acids_by_type.pdf"), plot = d)

#######################################
# ED Fig. 8e - GCA concentrations
#######################################

stat_test_gca <- df2plot_d %>%
  wilcox_test(Glycocholic.acid ~ Type, p.adjust.method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "Type") %>%
  filter((group1 == "Device 1" & group2 %in% c("Device 4", "Stool")) | (group1 == "Device 4" & group2 == "Stool"))

(gca <- ggplot(df2plot_d, aes(y=log10(Glycocholic.acid), x=Type))+
    geom_boxplot(aes(fill=Type), outlier.shape = NA, alpha = 0.9)+
    geom_beeswarm(alpha=0.4)+
    scale_fill_manual(values=CapTypeAndStoolColors)+
    labs(y=expression('log'[10]*'(GCA concentration)'), x='',color='') +
    stat_pvalue_manual(stat_test_gca, 
                       label = "p.adj.signif",
                       tip.length = 0.01,
                       y.position = c(6.7, 7, 6.7), 
                       bracket.shorten = 0.05) +    
    theme(legend.position='none',
          axis.title.x=element_blank()))

ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_8e_gca_conc_by_type.pdf"), plot = gca)

#######################################
# ED Fig. 8f - GCA correlations with different ASVs
#######################################

# Generate list of significantly negatively correlated ASVs (df2plot.filt) in Capsule samples
ddff <- subset_samples(ps.bs.raw, Set %in% c('2','3','4','5')) %>%
  subset_samples(., !Type %in% c('Stool')) %>%
  filter_taxa(., function(x) sum(x > 3) > (0.05*length(x)), TRUE) %>%
  transform_sample_counts(., function(x) {log2(x+1)}) %>%
  psmelt()

## FIRST LOOK AT GCA

bile_acid<-'Glycocholic.acid'
df_f <- ddff %>% filter(Abundance > 0)
df_f$bile_acid<-df_f[,bile_acid]
# split the data by ASV
B <- split(df_f, df_f$OTU)
# Calculate the correlation in all data.frames using lapply 
M <- lapply(B, function(x) cor.test(x$Abundance, log10(x$bile_acid), method='pearson')) 
# Get ASVs that are significant
pvals<-lapply(M, function(x) x$p.value)
pvals<-as.data.frame(pvals)
pvals.corrected<-p.adjust(pvals,method="hochberg") %>% as.data.frame() %>% filter(. < 0.01)
sig<-rownames(pvals.corrected)
# Slope 
r<-lapply(M, function(x) x$estimate)       
filtered<-data.frame(r) %>% t() %>% as.data.frame() %>% filter(cor < 0) # filter for negative correlations
filtered$OTU<-rownames(filtered)
sig.r <- filtered %>% filter(OTU %in% sig)
# Filter dataframe to only include significant ASVs              
df_f.filt <- df_f %>%
  filter(OTU %in% sig.r$OTU) %>%
  select(OTU, Phylum, Family, Genus, Species) %>%
  unique() %>%
  mutate(BA = bile_acid) %>%
  left_join(.,sig.r, by='OTU')
df_f.filt


# Make a dataframe to plot correlations in capsule and stool
df2plot_f<- psmelt(transform_sample_counts(ps.bs, function(x) {log2(x+1)})) %>%
  filter(OTU %in% df_f.filt$OTU) %>% #Sig negatively correlated ASVs
  mutate(Location = ifelse(Type %in% c('Device 1','Device 2','Device 3','Device 4'),'Device','Stool')) %>%
  select(OTU, Main, Type, Location, Abundance, Subject, Genus, Species, Glycocholic.acid, reads_16s) %>%
  mutate(title_name = paste0(Genus, ' ', Species)) %>%
  filter(Abundance > 0)

(gca_cor<-ggplot(df2plot_f, aes(x=Abundance, y=log10(Glycocholic.acid), color=Location))+
    geom_point(size=3,alpha=0.1)+
    geom_smooth(method='lm', alpha=0.2, aes(fill = Location))+
    stat_cor(method = 'pearson',p.accuracy = 0.001, r.accuracy = 0.01,
             label.x.npc = 0.5,label.y.npc = "top")+
    #ylim(0,15)+
    scale_color_manual(labels = paste("<span style='color:",
                                      CapAndStoolColors,"'>",
                                      levels(df2plot_f$Location),
                                      "</span>"),
                       values = CapAndStoolColors, 
                       guide = "none") +  
    scale_fill_manual(labels = paste("<span style='color:",
                                     CapAndStoolColors,"'>",
                                     levels(df2plot_f$Location),
                                     "</span>"),
                      values = CapAndStoolColors, 
                      guide = "none") + 
    labs(x=expression('log'[2]*'(ASV count)'), 
         y=expression('log'[10]*'(GCA concentration)'),
         color = "")+
    facet_wrap(~title_name + Location, nrow = 1) +
    theme(legend.position = "bottom",
          legend.margin = margin(),
          strip.text = element_text(size = 10, face = 4), 
          legend.spacing.x = unit(1.0, 'pt'),
          legend.text=element_markdown(size=10), 
          axis.text = element_text(size = 8),
          plot.margin = margin(5,5,5,5)))+
  guides(color = guide_legend(nrow=1, 
                              override.aes = list(label="", face=NA, shape = NA, fill=NA,linetype=NA))) ## This line removes the point on the legend for the color

ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_8f_gca_asv_corr.pdf"), plot = gca_cor, height = 5, width = 14)

#######################################
# ED Fig. 8g - TCDCA concentrations
####################################### 

df2plot_g <- df2plot_d %>%
    mutate(Taurochenodeoxycholic.acid = log10(Taurochenodeoxycholic.acid +1))

stat_test <- df2plot_g %>%
  wilcox_test(Taurochenodeoxycholic.acid ~ Type, p.adjust.method = "bonferroni") %>%
  add_significance() %>%
  filter(!p.adj.signif %in% 'ns') %>%
  add_xy_position(x = "Type") 

(tcdca <- ggplot(df2plot_g, aes(y=Taurochenodeoxycholic.acid, x=Type))+
    geom_boxplot(aes(fill=Type), outlier.shape = NA, alpha = 0.9)+
    geom_beeswarm(alpha=0.4, cex=0.5)+
    scale_fill_manual(values=CapTypeAndStoolColors)+
    labs(y=expression('log'[10]*'(TCDCA concentration)'), x='',color='') +
    stat_pvalue_manual(stat_test, 
                       label = "p.adj.signif",
                       tip.length = 0.01,
                       y.position = c(6,6.4,6.8,7.2,7.6), 
                       bracket.shorten = 0.05) +    
    theme(legend.position='none',axis.text = element_text(size = 10),axis.title.x=element_blank()))

ggsave(paste0(fig_dir_ed_subpanels, 'ED_Fig_8g_tcdca_conc.pdf'), tcdca, width = 4, height = 3)         

#######################################
# ED Fig. 8h,i - GCDCA, GDCA, TCDCA, TDCA correlations with different ASVs
#######################################

df_test <- ddff %>% filter(Abundance > 0)

# bile_acid_gcdca<-'Glycochenodeoxycholic.acid'
# df_test$bile_acid_gcdca<-df_f[,bile_acid_gcdca]

    

# split the data 
B <- split(df_test, df_test$OTU)


##################
##Glycochenodeoxycholic.acid
##################
# Calculate the correlation in all data.frames using lapply 
M <- lapply(B, function(x) cor.test(x$Abundance, log10(x$Glycochenodeoxycholic.acid), method='pearson'))

# Get ASVs that are significant
pvals<-lapply(M, function(x) x$p.value)
pvals<-as.data.frame(pvals)
pvals.corrected<-p.adjust(pvals,method="hochberg") %>% as.data.frame() %>% filter(. < 0.01)
sig<-rownames(pvals.corrected)

# Get negative slope for deconjugation
r<-lapply(M, function(x) x$estimate)
r<-as.data.frame(r)
filtered<-r[,r < 0]
sig.r<-colnames(filtered)             

# Filter dataframe to only include significant ASVs              
df_gcdca.filt <- df_test %>%
  filter(OTU %in% sig) %>%
  filter(OTU %in% sig.r) %>%
  filter(Glycochenodeoxycholic.acid > 1) %>%
  mutate(label = paste0 (Genus, " ", Species))
# **********************************************
### No significant ASV correlations with Glycochenodeoxycholic.acid
# **********************************************

##################
#### GLYCHODEOXYCHOLIC ACID              
##################
# Calculate the correlation in all data.frames using lapply 
M <- lapply(B, function(x) cor.test(x$Abundance, log10(x$Glycodeoxycholic.acid), method='pearson'))
            
# Get ASVs that are significant
pvals<-lapply(M, function(x) x$p.value)
pvals<-as.data.frame(pvals)
pvals.corrected<-p.adjust(pvals,method="hochberg") %>% as.data.frame() %>% filter(. < 0.01)
sig<-rownames(pvals.corrected)

# Get negative slope for deconjugation
r<-lapply(M, function(x) x$estimate)
r<-as.data.frame(r)
filtered<-r[,r < 0]
sig.r<-colnames(filtered)             

# Filter dataframe to only include significant ASVs              
df_gdca.filt <- df_test %>%
    filter(OTU %in% sig) %>%
    filter(OTU %in% sig.r) %>%
  mutate(label = paste0 (Genus, " ", Species))

# **********************************************
### No significant ASV correlations with Glycodeoxycholic.acid
# **********************************************

##################
#### TAUROCHENODEOXYCHOLIC ACID  
##################
# Calculate the correlation in all data.frames using lapply 
M <- lapply(B, function(x) cor.test(x$Abundance, log10(x$Taurochenodeoxycholic.acid), method='pearson'))

# Get ASVs that are significant
pvals<-lapply(M, function(x) x$p.value)
pvals<-as.data.frame(pvals)
pvals.corrected<-p.adjust(pvals,method="hochberg") %>% 
  as.data.frame() %>% 
  filter(. < 0.01)
sig<-rownames(pvals.corrected)

# Get negative slope for deconjugation
r<-lapply(M, function(x) x$estimate)
r<-as.data.frame(r)
filtered<-r[,r < 0]
sig.r<-colnames(filtered)             

# Filter dataframe to only include significant ASVs              
df_tcdca.filt <- df_test %>%
  filter(OTU %in% sig) %>%
  filter(OTU %in% sig.r) %>%
  mutate(label = paste0 (Genus, " ", Species))


(h <- ggplot(df_tcdca.filt, aes(x=Abundance, y=log10(Taurochenodeoxycholic.acid))) + 
    geom_point() + 
    geom_smooth(method=lm, se=FALSE) +
    labs(x=expression('log'[2]*'(ASV count)'), 
         y=expression('log'[10]*'(TCDCA concentration)')) +
    stat_cor(method = "pearson", label.x.npc = "center") +
    facet_wrap(~label) +
    theme(strip.text = element_text(size=12, face = 4),
          axis.text = element_text(size = 12)))

ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_8h_tcdca_asv_corr.pdf"), plot = h, height = 5, width = 14)

##################          
#### TAURODEOXYCHOLIC ACID              
##################
# Calculate the correlation in all data.frames using lapply 
M <- lapply(B, function(x) cor.test(x$Abundance, log10(x$Taurodeoxycholic.acid), method='pearson'))
            
# Get ASVs that are significant
pvals<-lapply(M, function(x) x$p.value)
pvals<-as.data.frame(pvals)
pvals.corrected<-p.adjust(pvals,method="hochberg") %>% as.data.frame() %>% filter(. < 0.01)
sig<-rownames(pvals.corrected)

# Get negative slope for deconjugation
r<-lapply(M, function(x) x$estimate)
r<-as.data.frame(r)
filtered<-r[,r < 0]
sig.r<-colnames(filtered)             

# Filter dataframe to only include significant ASVs              
df_tdca.filt <- df_test %>%
    filter(OTU %in% sig) %>%
    filter(OTU %in% sig.r) %>%
  mutate(label = paste0 (Genus, " ", Species))


(i <- ggplot(df_tdca.filt, aes(x=Abundance, y=log10(Taurodeoxycholic.acid))) + 
    geom_point()  + geom_smooth(method=lm, se=FALSE) +
    labs(x=expression('log'[2]*'(ASV count)'), 
         y=expression('log'[10]*'(TDCA concentration)')) +
    stat_cor(method = "pearson", label.x.npc = "center") +
    facet_wrap(~label+OTU)+
    theme(strip.text = element_text(size=12, face = 4),
    axis.text = element_text(size = 12)))

ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_8i_tdca_asv_corr.pdf"), plot = i, height = 5, width = 14)


#######################################
#  Put everything together for final figure
#######################################
abc <- plot_grid(bsh_uniqueBySubj_plot,
          plot_grid(plotlist = list(bsh_coverage_hist, bsh_rank_boxplot), ncol = 1, rel_heights = c(1, 1.2), align = "v", labels = c("b", "c"), label_size = 14),
          rel_widths = c(1.3, 1), rel_heights = c(1, 1.5), align = "v", axis = "br", labels = c("a", ""), label_size = 14)

abc_fh <- plot_grid(abc, gca_cor, h,
                    ncol = 1,
                    rel_heights = c(1.2, 0.6, 0.6),
                    labels = c("", "f", "h"),
                    label_size = 14)
degi <- plot_grid(plotlist = list(d, gca, tcdca, i), ncol = 1, rel_heights = c(1,1,1,1), align = "v", labels = c("d", "e", "g", "i"), label_size = 14)

all <- plot_grid(abc_fh, degi, ncol = 2, rel_widths = c(1, 0.6), labels = NULL)

ggsave(paste0(fig_dir, "ED_Figure_8.pdf"), width = 14, height = 14)
