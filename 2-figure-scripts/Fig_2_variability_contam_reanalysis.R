#######################################
# Capsule Study - Gen1 
#
# Generate Figure 2 for manuscript
# Redo analysis of microbial variability and patchiness
# while removing potentially contaminated capsules
#
# Author: Jess Grembi
#######################################


rm(list=ls())
# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))


(ps <- readRDS(clean_phyloseq_object) %>%
  subset_samples(!drop_16s & (Set %in% c('2','3','4','5') & !possible_contam) | Set %in% c("Stool", "Saliva", "6","7","8","9","10","11","12","13")) %>%
  filter_taxa(., function(x) sum(x > 3) > (0.05*length(x)), TRUE) %>% # Gets rid of all taxa not found at a count of 3 in at least 1% of samples
    transform_sample_counts(function(x) {log2(x + 1)})) 


df_samples <- ps@sam_data %>%
  data.frame

df_samples_min <- df_samples %>%
  select(id_16s, Subject, Set:Volume, reads_16s, location, experiment)

method = "canberra"
dist <- phyloseq::distance(ps, method = method, type = "samples") #binary = T
dist_df <- as.matrix(dist) %>%
  data.frame() %>%
  rownames_to_column("V1") %>%
  pivot_longer(cols = starts_with("S"), names_to = "Sample2", values_to = "dist") %>%
  rename(Sample1 = V1) %>%
  mutate(Sample1 = factor(Sample1, levels = ps@sam_data$id_16s, ordered = T),
         Sample2 = factor(Sample2, levels = ps@sam_data$id_16s, ordered = T)) %>%
  filter(Sample1 != Sample2, Sample1 < Sample2) %>%
  mutate(Sample1 = as.character(Sample1),
         Sample2 = as.character(Sample2)) %>%
  ## Check this code, I think we might have a few repeated lines due to duplicate id_16s entries in df_samples
  left_join(df_samples_min, by = c("Sample1" = "id_16s")) %>%
  left_join(df_samples_min, by = c("Sample2" = "id_16s"), suffix = c("_1", "_2")) %>%
  filter(Set_1 != 1, Set_2 != 1, location_1 == location_2) %>%
  mutate(locations = paste0(location_1, "_", location_2),
         location_cat = factor(locations, levels = c("Stool_Stool", "Saliva_Saliva", "Capsule_Capsule")))

dist_long_withinSubj <- dist_df %>%
  filter(Subject_1 == Subject_2, Set_1 %in% c("2", "3", "4", "5", "Stool", "Saliva"), Set_2 %in% c("2", "3", "4", "5", "Stool", "Saliva"), location_1 == location_2) %>%
  mutate(subject = Subject_1) %>%
  group_by(subject, location_cat) %>%
  summarise(mean_dist = mean(dist))

dist_long_acrossSubj <- dist_df %>%
  filter(Subject_1 != Subject_2, Set_1 %in% c("2", "3", "4", "5", "Stool", "Saliva"), Set_2 %in% c("2", "3", "4", "5", "Stool", "Saliva"), location_1 == location_2) %>%
  mutate(subject = ifelse(as.numeric(Subject_1) < as.numeric(Subject_2), paste0(Subject_1, "_", Subject_2), paste0(Subject_2, "_", Subject_1))) %>%
  group_by(subject, location_cat) %>%
  summarise(mean_dist = mean(dist))



dist_long <- bind_rows(dist_long_withinSubj, dist_long_acrossSubj) %>%
  mutate(subj_comp = factor(ifelse(grepl("_", subject), "Across subjects", "Within subject"), levels = c("Within subject", "Across subjects")))


## Significance test looking at Stool-Stool, Saliva-Saliva, and Capsule-Capsule within subjects compared to across subjects.
# All are significantly different.
sig_test <- dist_long %>% 
  group_by(location_cat) %>%
  do(w=wilcox.test(mean_dist~subj_comp, data = ., paired = F)) %>%
  summarise(location_cat, Wilcox = w$p.value) %>%
  mutate(p_short = ifelse(Wilcox < 0.0001, "<0.0001", ifelse(Wilcox < 0.001, "<0.001", ifelse(Wilcox < 0.01, "<0.01", "NS"))), 
         p_sym = ifelse(Wilcox < 0.0001, "***", ifelse(Wilcox < 0.001, "**", ifelse(Wilcox < 0.01, "*", "NS"))))

(plot_bigPic <- ggplot(dist_long %>%
                         mutate(location_cat = factor(location_cat, levels = c("Saliva_Saliva", "Capsule_Capsule", "Stool_Stool"), labels = c("Saliva", "Capsule", "Stool"))), 
                       aes(x = location_cat, y = mean_dist, color = location_cat)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.28,  alpha = 0.3, shape = 1) +
    stat_compare_means(comparisons = list(c("Capsule", "Stool"), c("Capsule", "Saliva")), 
                       method = "wilcox.test",
                       label = "p.signif",
                       vjust = 0.5) + 
    # geom_text(data = sig_test, aes(x = location_cat, y = max(dist_long$dist) + 0.02, label = p_sym), inherit.aes = FALSE, fontface = 'bold')+
    scale_color_manual(values = c("black", CapAndStoolColors), guide = "none") + #"grey25", "#78B7C5", "#F21C00"
    # scale_shape_manual(values = c(1,2)) +
    labs(y = paste0(toupper(substr(method, 1, 1)), substr(method, 2, nchar(method)), " distance"), x = "", color = "") + 
    coord_flip() +
    facet_grid(subj_comp~.) + 
    theme(strip.text.y = element_text(size = 14)))

# ggsave(filename = paste0(fig_dir, "microbial_variability_across_locations.pdf"), plot = plot_bigPic)

medianCanberra_dist_bySubj <- dist_long %>%
  filter(subj_comp == "Within subject") %>%
  group_by(subject, location_cat) %>%
  summarise(median_canberra_dist = median(mean_dist))


##  Where is the variability in capsules coming from (type, set, time of day, etc)?
# For this analysis we will use only capsule-capsule comparisons from within the same subject
dist_capsules <- dist_df %>%
  filter(location_cat == "Capsule_Capsule", Subject_1 == Subject_2)

dist_repeat <- dist_capsules %>%
  filter(Type_1 == Type_2, Set_1 == Set_2, experiment_1 == experiment_2, experiment_1 == "repeatability") %>%
  mutate(subject = Subject_1,
         group = "Same capsule type-Same time",
         group_name = paste0(Type_1, "_", Set_1)) %>%
  group_by(subject, group_name, group) %>% #group_name
  summarise(n = n(),
            mean_dist = mean(dist)) %>%
  mutate(experiment = "Repeatability")


dist_types <- dist_capsules %>%
  filter(Type_1 == Type_2, experiment_1 == experiment_2, experiment_1 != "repeatability") %>%
  mutate(subject = Subject_1,
         group = "Same capsule type-Different time",
         group_name = Type_1) %>%
  group_by(subject, group) %>% #group_name
  summarise(n = n(),
            mean_dist = mean(dist)) %>%
  mutate(experiment = "Main study")

dist_sets <- dist_capsules %>%
  filter(Set_1 == Set_2, experiment_1 == experiment_2, experiment_1 != "repeatability") %>%
  mutate(subject = Subject_1,
         group = "Different capsule type-Same time",
         group_name = Set_1) %>%
  group_by(subject, group) %>% #group_name
  summarise(n = n(),
            mean_dist = mean(dist)) %>%
  mutate(experiment = "Main study")

dist_sets_and_types <- dist_capsules %>%
  filter(Set_1 != Set_2, Type_1 != Type_2, experiment_1 == experiment_2, experiment_1 != "repeatability") %>%
  mutate(subject = Subject_1,
         group = "Different capsule type-Different time",
         group_name = paste0("Set",Set_1,"_",Set_2,"_Type", Type_1,"_", Type_2)) %>%
  group_by(subject, group) %>%
  summarise(n = n(),
            mean_dist = mean(dist)) %>%
  mutate(experiment = "Main study")

dist_var <- dist_repeat %>%
  bind_rows(dist_types) %>%
  bind_rows(dist_sets) %>%
  bind_rows(dist_sets_and_types)

(variability <- ggplot(dist_var %>%
                         mutate(group_name = factor(ifelse(group_name %in% c("Capsule 1", "Capsule 2", "Capsule 3", "Capsule 4"), group_name, "Other"), levels = c("Capsule 1", "Capsule 2", "Capsule 3", "Capsule 4", "Other"))),
                       aes(x = group, y = mean_dist, color = experiment)) + 
   geom_boxplot(outlier.shape = NA) + 
   geom_jitter(width = 0.3, alpha = 0.4, shape = 1) +
   stat_compare_means(comparisons = list(c("Same capsule type-Same time", "Same capsule type-Different time"), c("Different capsule type-Different time", "Different capsule type-Same time")), 
                      method = "wilcox.test", 
                      label = "p.signif",
                      vjust = 0.25) + 
   coord_flip() +
   labs(y = paste0(toupper(substr(method, 1, 1)), substr(method, 2, nchar(method)), " distance"), 
        x = "", color = "Evaluation") + 
   theme(legend.position = "bottom",
         legend.margin = margin(), 
         plot.margin = margin(t = fig_padding, l = fig_padding, r = fig_padding, b = fig_padding, unit = "pt")) +
   scale_color_viridis_d(end = 0.6)) 

# ggsave(filename = paste0(fig_dir, "microbial_variability_type_vs_temporal.pdf"), plot = variability)




################################
## Alpha diversity 
################################
# Across capsule types

df_capsule <- ps@sam_data %>%
  data.frame %>%
  filter(location != "Control") %>%
  mutate(Type = factor(Type, levels = c("Saliva", "Capsule 1", "Capsule 2", "Capsule 3", "Capsule 4", "Stool")))


(alpha_div <- ggplot(df_capsule, 
      aes(x = Type, y = Shannon, fill = Type)) + 
  geom_boxplot(outlier.shape = NA) +
  ggbeeswarm::geom_beeswarm(alpha = 0.6, color = "black") +
  # geom_jitter(width = 0.3, alpha = 0.6, shape = 1) + 
  stat_compare_means(comparisons = list(c("Capsule 1", "Capsule 4"), c("Capsule 4", "Stool"), c("Saliva", "Capsule 1")),
                     method = "wilcox.test",
                     label.y = c(4.05, 4.25, 4.25),
                     # step_increase = -0.5,       
                     label = "p.signif") + 
  # ylim(0,5) +
  labs(y = "Shannon diversity", x = "") + 
  theme(legend.position = "none",
        legend.margin = margin(), 
        plot.margin = margin(t = fig_padding, l = fig_padding, r = fig_padding, b = fig_padding, unit = "pt")) +
  scale_fill_manual(values = c("black", CapTypeAndStoolColors)))

# ggsave(paste0(fig_dir, "alpha_diversity.pdf"), width = 7, height = 4)


################################
## Patchiness  
################################
## Call in phyloseq object again, but don't transform sample counts
ps_untransformed <- readRDS(clean_phyloseq_object) %>%
  subset_samples(!drop_16s & (!possible_contam | Set %in% c("Saliva", "Stool")))

top_rank <- data.frame(id_16s = rownames(ps_untransformed@sam_data), 
                       top_rank = unlist(lapply(1:dim(ps_untransformed@otu_table)[1], function(x) max(ps_untransformed@otu_table[x,])))) %>%
  left_join(ps_untransformed@sam_data %>%
              data.frame, by = "id_16s") %>%
  mutate(top_tax_abnd = top_rank/reads_16s, 
         Type = factor(Type, levels = c("Saliva", "Capsule 1", "Capsule 2", "Capsule 3", "Capsule 4", "Stool")))

(patchiness <- ggplot(top_rank, aes(top_tax_abnd, color = Type)) + 
    stat_ecdf(aes(alpha = 0.2), geom = "point") + 
    labs(x = "Relative abundance of most abundant ASV", y = "ECDF", color = "Location") +
    scale_alpha(guide = "none") + 
    # guides(colour = guide_legend(nrow = 1)) +
    scale_color_manual(labels = paste("<span style='color:",
                                      c("black", CapTypeAndStoolColors),
                                      "'>",
                                      levels(top_rank$Type),
                                      "</span>"),
                       values = c("black", CapTypeAndStoolColors)) +  
    theme(legend.position = "bottom",
          legend.margin = margin(), 
          legend.text=element_markdown(size=14),
          plot.margin = margin(t = fig_padding, l = fig_padding, r = fig_padding, b = fig_padding, unit = "pt")) + 
    guides(color = guide_legend(override.aes = list(shape = NA))))

# ggsave(paste0(fig_dir, "capsule_patchiness.pdf"), plot = patchiness, width = 7, height = 4)


plot_grid(plot_grid(plot_bigPic, variability, labels = c("A", "B"), label_size = 16), 
          plot_grid(alpha_div, patchiness, labels = c("C", "D"), label_size = 16), rel_heights = c(1.6, 1), nrow = 2)

ggsave(paste0(fig_dir, "contam_reanalysis/Figure_2_variability.pdf"), width = 10, height = 7)

