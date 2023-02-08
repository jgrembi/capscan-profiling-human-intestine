#######################################
# Capsule Study - Gen1 
#
# Figure 3 and Extended Data Figure 6
#  Viromics
#
# Author: Original analysis and code written by Handuo Shi in Matlab
#         Translated to R code for inclusion in this reproducible workflow by Jess Grembi 
#######################################



rm(list=ls())
# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# df <- R.matlab::readMat(paste0(data_dir, "viral_contig_prophage_coverage_summary.mat"))

reads_threshold <- 1e6
# variables_all_meta <- read.csv(paste0(data_dir, "variables_all_metagenomics_corrected.csv")) #%>%
variables_all_meta <- read_excel(paste0(data_dir, "sample_metadata_matlab_export.xlsx"),
                                 col_types = c("text", "text", "numeric", "numeric", "numeric", "text", "text", "numeric", "text", "numeric", "numeric", "numeric", "numeric")) %>%
  mutate(Abrreviation = gsub("'", "", Abrreviation),
         Subject = factor(Subject, levels = c(1:15)),
         location = gsub("'", "", `Putative location`), 
         location = ifelse(location == "stool", "Stool", ifelse(location == "", "Control", location)), 
         Set = ifelse(location == "Stool", "Stool", ifelse(location == "Saliva", "Saliva", (gsub("Set", "", gsub("'", "", Set))))))


df_samples <- readRDS(sample_data) %>%
  filter(Set %in% c("2", "3", "4", "5", "Stool", "Saliva"))

variables_all_meta_new <- read_excel(paste0(data_dir, "sample_metadata_matlab_export.xlsx"),
                                     col_types = c("text", "text", "numeric", "numeric", "numeric", "text", "text", "numeric", "text", "numeric", "numeric", "numeric", "numeric")) %>%
  mutate(meta_samplename = gsub("'", "", Abrreviation),
         `Putative location` = ifelse(`Putative location` == "stool", "Stool", `Putative location`)) %>%
  select(-reads_humanRemoved) %>%
  left_join(df_samples %>% select(meta_samplename, reads_meta_filtered, drop_meta)) %>%
  select(-meta_samplename) %>%
  dplyr::rename(reads_humanRemoved = reads_meta_filtered)

write.csv(variables_all_meta_new, file = paste0(data_dir, "sample_metadata_matlab_export_updatedHumanReadsRemoved.csv"))

phage_names <- read_excel(paste0(data_dir, "phage_contig_names.xlsx"), col_names = F) %>%
  mutate(phage_name = gsub("'", "", `...1`))

phage_depth <- read_csv(paste0(data_dir, "phage_depth.csv"), col_names = phage_names$phage_name)
rownames(phage_depth) <- variables_all_meta$Abrreviation

phage_reads <- read_excel(paste0(data_dir, "phage_reads.xlsx"), col_names = variables_all_meta$Abrreviation)
rownames(phage_reads) <- phage_names$phage_name

phage_reads_long <- phage_reads %>%
  rownames_to_column("phage_name") %>%
  pivot_longer(cols = !phage_name, names_to = "meta_samplename", values_to = "phage_reads")
  
# canberra <- read.csv(paste0(data_dir, "canberradist.csv"))
# 
# pcoa <- read_excel(paste0(data_dir, "pcoa.xlsx"))
prophage_ind <- read_csv(paste0(data_dir, "prophage_induction_state.csv"), col_names = phage_names$phage_name)
rownames(prophage_ind) <- variables_all_meta$Abrreviation

prophage_ind_long <- prophage_ind %>%
  rownames_to_column("meta_samplename") %>%
  pivot_longer(cols = `A10_2-k141_36643`:`P15_1-k141_98150`, names_to = "phage_name", values_to = "prophage_ind")

phage_df <- df_samples %>%
  select(meta_samplename, location, Subject:Main,reads_meta_trimmed, reads_meta_filtered, drop_meta) %>%
  left_join(phage_depth %>% 
              rownames_to_column("meta_samplename"), by = "meta_samplename") %>%
  pivot_longer(cols = `A10_2-k141_36643`:`P15_1-k141_98150`, names_to = "phage_name", values_to = "phage_depth") %>%
  left_join(prophage_ind_long, by = c("meta_samplename", "phage_name")) %>%
  mutate(phage_presence = ifelse(phage_depth > 0, 1, 0)) %>%
  left_join(phage_reads_long, by = c("meta_samplename", "phage_name")) %>%
  filter(reads_meta_filtered > reads_threshold, 
         Set %in% c("Stool", "Saliva", "2", "3", "4", "5"),
         !drop_meta) %>%
  unique()
#
#-----------------------------------
# Find difference in my sampleframe and the one Handuo used
#-----------------------------------
handuo_sampleframe <- read_excel(paste0(data_dir, "capsule_samples_finalanalyses.xlsx"), col_names = F) %>%
  mutate(samplename = gsub("'", "", `...2`))

capsules_inc <- phage_df %>%
  filter(location == "Capsule") 

me_not_handuo <- setdiff(unique(capsules_inc$meta_samplename), handuo_sampleframe$samplename)

handuo_not_me <- setdiff(handuo_sampleframe$samplename, unique(capsules_inc$meta_samplename))

me <- df_samples %>%
  filter(meta_samplename %in% me_not_handuo)
me_not_handuo
# These are all Subject 3 samples and should be included.  We did not include metagenomic sequencing for Subj3 second set of samples so no need to remove anything

handuo <- df_samples %>%
  filter(meta_samplename %in% handuo_not_me)
handuo_not_me
## These are 4 samples that were sequenced twice (not sure why).  We chose to use the sample with higher shotgun sequencing depth, which were NOT these 4)
#-----------------------------------

phage_location <- phage_df %>%
  group_by(location, phage_name) %>%
  summarise(n = n(), 
            n_phage = sum(phage_depth >= 1)) %>%
  ungroup() %>%
  mutate(phage_present = ifelse(n_phage > 0, T, F)) %>%
  select(phage_name, location, phage_present) %>%
  pivot_wider(id_cols = c("phage_name"), names_from = "location", values_from = "phage_present")

## This describes the n's for this figure
phage_df %>%
  group_by(location, phage_name) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  select(location, n) %>%
  unique()

(venn_stool_capsule_phage <- plot(euler(phage_location %>% select(Capsule, Stool), quantities = T, shape = "ellipse"), 
                                 fill = CapAndStoolColors, 
                                 alpha = 0.7,
                                 quantities = T))


ggsave(paste0(fig_dir_main_subpanels, "Fig_3a_subpanel_phage_presence_devices_stool.pdf"), plot = venn_stool_capsule_phage)

#------------------------------------------
# PCoA
#------------------------------------------
sam_dat <- phage_df %>%
  select(meta_samplename, location, reads_meta_filtered, Type) %>%
  unique() %>%
  mutate(rowname = meta_samplename) %>%
  column_to_rownames("rowname")


votu <- phage_df %>%
  select(meta_samplename, phage_depth, phage_name) %>%
  pivot_wider(names_from = "phage_name", values_from = "phage_depth") %>%
  column_to_rownames("meta_samplename")


ps_votu <- phyloseq(otu_table(t(votu), taxa_are_rows = T), sample_data(sam_dat)) %>%
  transform_sample_counts(function(x) x/sum(x))
pcoa_canberra <- ordinate(ps_votu,  method = "MDS", distance = "canberra")


canberra <- distance(ps_votu, method = "canberra", type = "samples")
d = as.matrix(canberra) %>%
  data.frame %>%
  rownames_to_column("meta_samplename") %>%
  select(meta_samplename, A11_1, A11_2)
var_exp <- get_evals(pcoa_canberra)$variance_exp
scores <- get_scores_votu(pcoa_canberra, sample_data(ps_votu)) %>%
  mutate(location = factor(location, levels = c("Saliva", "Capsule", "Stool"), labels = c("Saliva", "Devices", "Stool")), 
         Type = factor(gsub("Capsule", "Device", Type), levels = c("Saliva", "Device 1", "Device 2", "Device 3", "Device 4", "Stool")))


(b <- ggplot(scores,  
             aes(x = PC1, y = PC2, color = location)) + #shape = subj_id
    geom_point(size = 3, alpha = 0.8) + 
    ## This code changes the color of the legend text so that it will mirror the points
    scale_color_manual(labels = paste("<span style='color:",
                                      c("black", CapAndStoolColors),
                                      "'>",
                                      levels(scores$location),
                                      "</span>"),
                       values = c("black", CapAndStoolColors), 
                       guide = "none") +    
    coord_fixed(var_exp[2]/var_exp[1]) +
    # scale_color_manual(values = c(CapTypeAndStoolColors, "black")) + 
    # scale_shape_manual(values = c(15, 17, 1)) +
    labs(x = paste0("Principal coordinate 1 [", round(var_exp[1],1), "% variance]"),
         y = paste0("Principal coordinate 2 [", round(var_exp[2], 1), "% variance]"),
         color = "Sample type") + 
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.title.y = element_textbox_simple(
            hjust = 0,
            orientation = "left-rotated",
            minwidth = unit(1, "in"),
            maxwidth = unit(2, "in"),
            padding = margin(4, 4, 2, 4),
            margin = margin(0, 0, 2, 0)),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.margin = margin(),
          legend.spacing.x = unit(1.0, 'pt'),
          ## This line is needed with the above (lines 151-157) to change the color of the legend text to match the points
          legend.text=element_markdown(size=14),
          legend.title = element_text(size = 14, vjust = 0.97),
          plot.margin = margin(t = fig_padding/5, l = fig_padding/5, r = fig_padding/5, b = fig_padding/5, unit = "pt")))# + 
    # guides(color = guide_legend(nrow = 6,
    #                             override.aes = list(shape = NA),
    #                             keywidth = unit(0, "pt")),
    #        shape = guide_legend(nrow = 3))) ## This line removes the point on the legend for the color

ggsave(filename = paste0(fig_dir_main_subpanels, "Fig_3b_subpanel_pcoa_canberra.pdf"), plot = b, height = 6, width = 8)

votu_01 <- phage_df %>%
  select(meta_samplename, phage_presence, phage_name) %>%
  pivot_wider(names_from = "phage_name", values_from = "phage_presence") %>%
  column_to_rownames("meta_samplename")


#------------------------------------------
# Sequencing depth
#------------------------------------------
## use medians because the raw phage_depth across all different samples is not normally distributed.
seq_depth_by_votu <- phage_df %>%
  # filter(phage_depth >= 1) %>%
  group_by(phage_name, location) %>%
  summarise(n = n(),
            median_depth = median(phage_depth/2, na.rm = T)) %>%
  pivot_wider(id_cols = phage_name, names_from = location, values_from = c("median_depth")) %>%
  mutate(both = (Stool > 0 & Capsule > 0)) %>%
  filter(both)

# How many vOTUs are used in this analysis
dim(seq_depth_by_votu)[1]

(ed_fig_6a <- ggplot(seq_depth_by_votu, 
       aes(x = log10(Stool), y = log10(Capsule))) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = c(0,0), linetype = "longdash") + 
  stat_cor(method = "pearson",
           cor.coef.name = "rho") +
  labs(x = expression('log'[10]*'(median reads in Stool)'), 
       y = expression('log'[10]*'(median reads in Devices)')))

ggsave(filename = paste0(fig_dir_ed_subpanels, "ED_Fig_6a_viralread_depth_stool_v_devices.pdf"), plot = ed_fig_6a)

seq_depth <- phage_df %>%
  filter(phage_depth >= 1) %>%
  group_by(meta_samplename, location, reads_meta_filtered, reads_meta_trimmed) %>%
  summarise(n = n(),
            viralSeq = sum(phage_reads)) %>%
  mutate(pct_readViral = viralSeq/reads_meta_filtered*100/2, 
         pct_trimmed = viralSeq/reads_meta_trimmed*100/2,
         location = factor(location, levels = c("Saliva", "Capsule", "Stool"), labels = c("Saliva", "Devices", "Stool")),
         pct_diff = pct_readViral - pct_trimmed)

summary(seq_depth$pct_diff)

table(seq_depth$location)
ggplot(seq_depth, aes(x = pct_readViral)) + 
  geom_histogram(bins = 50) + 
  facet_wrap(~location)

seq_depth_stats <- seq_depth %>%
  ungroup() %>%
  group_by(location) %>%
  summarise(median = median(pct_readViral),
            q05 = quantile(pct_readViral, probs = 0.05),
            q25 = quantile(pct_readViral, probs = 0.25),
            q75 = quantile(pct_readViral, probs = 0.75),
            q95 = quantile(pct_readViral, probs = 0.95))

(ed_fig_6b <- ggplot(seq_depth, 
                 aes(x = location, y = pct_readViral, color = location, fill = location)) + 
    geom_beeswarm(alpha = 0.4) + 
    geom_violin(alpha = 0.3, scale= 'width') + 
    geom_segment(data = seq_depth_stats, aes(x = location, xend = location, y = q25, yend = q75), inherit.aes = F, linewidth = 2, color = "gray35") + 
    geom_segment(data = seq_depth_stats, aes(x = location, xend = location, y = q05, yend = q25), inherit.aes = F, color = "gray35") + 
    geom_segment(data = seq_depth_stats, aes(x = location, xend = location, y = q75, yend = q95), inherit.aes = F, color = "gray35") + 
    geom_point(data = seq_depth_stats, aes(x = location, y = median), fill = "white", shape = 21) +
    scale_color_manual(values = c("black", CapAndStoolColors), guide = "none") + 
    scale_fill_manual(values = c("black", CapAndStoolColors), guide = "none") + 
    stat_compare_means(method = "wilcox.test", comparisons = list(c("Stool", "Devices"), c("Devices", "Saliva"), c("Stool", "Saliva"))) + 
    labs(x = "", y = "Reads mapped as viral (%)"))
  
ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_6b_percent_viral_reads.pdf"), plot = ed_fig_6b)

#------------------------------------------
# Prophage induction
#------------------------------------------
ind_prophage <- phage_df %>%
  group_by(Subject, meta_samplename, location) %>%
  summarise(n = n(),
            n_indPhage = sum(prophage_ind == 1)) %>%
  mutate(location = factor(location, levels = c("Saliva", "Capsule", "Stool"), labels = c("Saliva", "Devices", "Stool"))) %>%
  left_join(df_samples %>% 
              filter(!drop_meta) %>%
              select(pH, meta_samplename, reads_meta_filtered), by = c("meta_samplename")) 


## This describes the n's for this figure
table(ind_prophage$location)

ggplot(ind_prophage) +
  geom_histogram(aes(x = n_indPhage), bins = 100) + 
  facet_wrap(~location, scales = "free")


ind_prophage_stats <- ind_prophage %>%
  ungroup() %>%
  group_by(location) %>%
  summarise(median = median(n_indPhage),
            q05 = quantile(n_indPhage, probs = 0.05),
            q25 = quantile(n_indPhage, probs = 0.25),
            q75 = quantile(n_indPhage, probs = 0.75),
            q95 = quantile(n_indPhage, probs = 0.95))


(c <- ggplot(ind_prophage, aes(x = location, y = n_indPhage, color = location, fill = location)) + 
  geom_beeswarm(alpha = 0.4) + 
  geom_violin(alpha = 0.3, scale= 'width') + 
  geom_segment(data = ind_prophage_stats, aes(x = location, xend = location, y = q25, yend = q75), inherit.aes = F, linewidth = 2, color = "gray35") + 
  geom_segment(data = ind_prophage_stats, aes(x = location, xend = location, y = q05, yend = q25), inherit.aes = F, color = "gray35") + 
  geom_segment(data = ind_prophage_stats, aes(x = location, xend = location, y = q75, yend = q95), inherit.aes = F, color = "gray35") + 
  geom_point(data = ind_prophage_stats, aes(x = location, y = median), fill = "white", shape = 21) +
  scale_color_manual(values = c("black", CapAndStoolColors), guide = "none") + 
  scale_fill_manual(values = c("black", CapAndStoolColors), guide = "none") + 
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Stool", "Devices"), c("Devices", "Saliva"), c("Stool", "Saliva"))) + 
  labs(x = "", y = "Number of induced prophages"))
  

ggsave(paste0(fig_dir_main_subpanels, "Fig_3c_subpanel_induced_prophage_by_location.pdf"), plot = c)



(ed_fig_6c <- ggplot(ind_prophage, aes(x = Subject, y = n_indPhage, color = location)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(alpha = 0.4, position = position_jitterdodge()) + 
    scale_color_manual(values = c("black", CapAndStoolColors), guide = "none") + 
    labs(y = "Number of induced prophages"))

ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_6c_ind_phage_by_subj.pdf"), plot = ed_fig_6c)

pH_corr <- ind_prophage %>% filter(location == "Devices")
cor.test(x = pH_corr$n_indPhage, y = pH_corr$pH, method = "spearman", exact = F)

(ed_fig_6d <- ggplot(pH_corr, aes(x = pH, y = n_indPhage)) + 
    geom_point() + 
    stat_cor(method = "spearman",
             cor.coef.name = "rho") + 
    labs(y = "Number of induced prophages", x = "pH of device sample"))

ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_6d_indProphage_by_pH.pdf"), plot = ed_fig_6d)



ind_prophage_location <- phage_df %>%
  group_by(location, phage_name) %>%
  summarise(n = n(), 
            n_indProphage = sum(prophage_ind == 1)) %>%
  # ungroup() %>%
  mutate(ind_prophage_present = ifelse(n_indProphage > 0, T, F)) %>%
  select(phage_name, location, ind_prophage_present) %>%
  pivot_wider(id_cols = c("phage_name"), names_from = "location", values_from = "ind_prophage_present")


(venn_stool_capsule_indProphage <- plot(euler(ind_prophage_location %>% select(Capsule, Stool), quantities = T, shape = "ellipse"), 
                                  fill = CapAndStoolColors, 
                                  alpha = 0.7,
                                  quantities = T))
ggsave(paste0(fig_dir_main_subpanels, "Fig_3d_subpanel_phage_presence_devices_stool.pdf"), plot = venn_stool_capsule_indProphage)



#----------------------------
# Figure 3
#----------------------------

plot_grid(plot_grid(venn_stool_capsule_phage, b, labels = c("a", "b"), label_size = 16, ncol = 2, rel_widths = c(0.6, 1)), 
          plot_grid(c, venn_stool_capsule_indProphage, labels = c("c", "d"), label_size = 16, ncol = 2, rel_widths = c(1,0.6)), nrow = 2)

ggsave(paste0(fig_dir, "Figure_3.pdf"), width = 14, height = 10)




#----------------------------
# Extended Data Figure 6
#----------------------------
plot_grid(plot_grid(ed_fig_6a, ed_fig_6b, labels = c("a", "b"), label_size = 16, ncol = 2), 
          plot_grid(ed_fig_6c, ed_fig_6d, labels = c("c", "d"), label_size = 16, ncol = 2, rel_widths = c(1,0.6)), nrow = 2)

ggsave(paste0(fig_dir, "ED_Figure_6.pdf"), width = 14, height = 10)
