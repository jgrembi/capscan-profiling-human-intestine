#######################################
# Capsule Study - Gen1 
#
# Extended Data Figure 3
#
# Author: Jess Grembi and Rebecca Culver
#######################################

rm(list=ls())
# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

#######################################
# ED Fig. 3a - Abundance of taxa at the phylum-level
#######################################
ps <- readRDS(clean_phyloseq_object) %>% 
  subset_samples(., drop_16s == F) %>%
  subset_samples(., Set %in% c("2", "3", "4", "5", "Stool")) %>%
  # filter_taxa(., function(x) sum(x > 3) > (0.05*length(x)), TRUE) %>% # Gets rid of all taxa not found at a count of 3 in at least 5% of samples (that's 14 samples)
  tax_glom("Phylum") %>%
  transform_sample_counts(function(x) {log2(x + 1)})

ps@sam_data$Type <- gsub("Capsule", "Device", ps@sam_data$Type)
ps@sam_data$location <- gsub("Capsule", "Devices", ps@sam_data$location)

# Some phylum for each sample
df2plot_a <- psmelt(ps) %>%
  filter(Abundance > 0)
tmp<-df2plot_a %>% select(Main, location) %>% unique()
table(tmp$location)
(a <- ggplot(df2plot_a, aes(x=location, y=Abundance)) +
    # geom_point() +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.28,  alpha = 0.3, shape = 1) +
    stat_compare_means(label = "p.signif", method="wilcox") +
    labs(x = "", y = expression('log'[2]*'(ASV count)')) + 
    facet_wrap(~Phylum, ncol = 5))

ggsave(filename = paste0(fig_dir_ed_subpanels, "ED_Fig_3a_phylum_comparisons.pdf"), plot = a, width=10, height=6)

#######################################
# ED Fig. 3b,c - Diversity for each subjects' devices by device type and set
#######################################

## don't filter ASVs before calculating alpha diversity
(ps_noFilt <- readRDS(clean_phyloseq_object) %>%
   subset_samples(!drop_16s & Set %in% c("2", "3", "4", "5")))  

set.seed(5502)
argList <- list(s = levels(ps_noFilt@sam_data$Subject), i = 1:1000)
crossArg <- cross_df(argList)

# ***************************************
### NOTE: This takes a while to run.  At least 35 min on my laptop
# ***************************************
alpha_div_repeated <- map2_dfr(crossArg$s, crossArg$i, function(s, i) {
  out = bySubj_rarefy(s = s, i = i, phyloseq_object = ps_noFilt)
  return(out)
})

alpha_div_annotated <- alpha_div_repeated %>%
  mutate(id_16s = ifelse(level == "sample", ID, NA)) %>%
  left_join(ps_noFilt@sam_data %>% 
              data.frame %>%
              select(id_16s, Type, Set), by = "id_16s") %>%
mutate(Type = factor(ifelse(level == "sample", gsub("Capsule", "Device", as.character(Type)), ID), levels = c("Device 1", "Device 2", "Device 3", "Device 4")),
       Set = ifelse(level == "sample", as.character(Set), ifelse(level == "set", ID, NA)),
       Subject = as.numeric(Subject)) 

alpha_div_samp_summary <- alpha_div_annotated %>%
  filter(level == "sample") %>%
  group_by(Subject, id_16s, Type, level) %>%
  summarise(n = n(), 
            meanShannon = mean(Shannon, na.rm = T),
            sdShannon = sd(Shannon, na.rm = T),
            lower_ci = quantile(Shannon, probs = 0.05, na.rm = T), 
            upper_ci = quantile(Shannon, probs = 0.95, na.rm = T))

alpha_div_subj_summary <- alpha_div_annotated %>%
  filter(level == "type") %>%
  mutate(Type = factor(gsub("Capsule", "Device", ID), levels = c("Device 1", "Device 2", "Device 3", "Device 4"))) %>%
  group_by(Subject, Type, level) %>%
  summarise(n = n(), 
            meanShannon = mean(Shannon, na.rm = T),
            sdShannon = sd(Shannon, na.rm = T),
            lower_ci = quantile(Shannon, probs = 0.05, na.rm = T), 
            upper_ci = quantile(Shannon, probs = 0.95, na.rm = T)) %>%
  mutate(Subject = as.numeric(Subject))

alpha_div_summary_type <- alpha_div_samp_summary %>%
  ungroup %>%
  select(-id_16s) %>%
  bind_rows(alpha_div_subj_summary) 


(types <- ggplot(alpha_div_summary_type, 
                 aes(y = Type, x = meanShannon, color = Type, alpha = level)) + 
    geom_point() +
    geom_errorbar(aes(xmin = lower_ci, xmax = upper_ci), size = 0.8, width = 0.5) +
    scale_color_manual(labels = paste("<span style='color:",
                                      CapType,
                                      "'>",
                                      levels(alpha_div_summary_type$Type),
                                      "</span>"),
                       values = CapType) + 
    scale_alpha_manual(values = c(0.2, 0.9), guide = "none") + 
    guides(color = guide_legend(override.aes = list(shape = NA, linetype = NA),
                                keywidth = unit(0, "pt"))) +
    facet_wrap(~Subject, ncol = 5) + 
    theme(legend.position = "top", 
          legend.margin = margin(),
          legend.spacing.x = unit(0.1, 'pt')) +
    labs(x = "Shannon diversity", y = "", color = "Sample type"))


ggsave(filename = paste0(fig_dir_ed_subpanels, "ED_Fig_3b_alpha_div_by_device_type.pdf"), plot = types)

###############################
### Looking at samples vs Sets
###############################

alpha_div_samp_set <- alpha_div_annotated %>%
  filter(level == "sample") %>%
  group_by(Subject, id_16s, Set, level) %>%
  summarise(n = n(), 
            meanShannon = mean(Shannon),
            sdShannon = sd(Shannon),
            lower_ci = quantile(Shannon, probs = 0.05, na.rm = T), 
            upper_ci = quantile(Shannon, probs = 0.95, na.rm = T))

alpha_div_set_summary <- alpha_div_annotated %>%
  filter(level == "set") %>%
  group_by(Subject, Set, level) %>%
  summarise(n = n(), 
            meanShannon = mean(Shannon),
            sdShannon = sd(Shannon),
            lower_ci = quantile(Shannon, probs = 0.05, na.rm = T), 
            upper_ci = quantile(Shannon, probs = 0.95, na.rm = T))

alpha_div_summary_set <- alpha_div_samp_set %>%
  ungroup %>%
  select(-id_16s) %>%
  bind_rows(alpha_div_set_summary) %>%
  mutate(Set = factor(as.numeric(Set), levels = c("5","4","3","2")))


(sets <- ggplot(alpha_div_summary_set, 
                aes(x = Set, y = meanShannon, color = Set, alpha = level)) + 
    geom_point() +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), size = 0.8, width = 0.5) +
    coord_flip() + 
    scale_color_viridis_d(end = 0.8) + 
    scale_alpha_manual(values = c(0.2, 0.9), guide = "none") + 
    theme(legend.position = "top") +
    facet_wrap(~Subject, ncol = 5) + 
    labs(y = "Shannon diversity", x = "Set", color = "Set"))

ggsave(filename = paste0(fig_dir_ed_subpanels, "ED_Fig_3c_alpha_div_by_set.pdf"), plot = sets)

#######################################
#  Put everything together for final figure
#######################################

plot_grid(a, types, sets, 
          nrow = 3, 
          labels = c("a", "b", "c"), 
          label_size = 16, 
          rel_heights = c(0.6, 1, 1))

ggsave(paste0(fig_dir, "ED_Figure_3.pdf"), width = 10, height = 14)
