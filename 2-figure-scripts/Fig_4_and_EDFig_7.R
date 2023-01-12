#######################################
# Capsule Study - Gen1 
#
# Generate Proteomics Figures for manuscript
#  Figure 4
#  Extended Data Figure 7
#
# Authors: Peter Triet and Florian Schober, modifications by Jess Grembi
#######################################

rm(list=ls())
# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

### -- Data files

d <- read_tsv(paste0(data_dir,"proteinGroups.tsv"))
# lfq <- read_tsv(paste0(data_dir,"proteinGroups_lfq.tsv"))
# meta <- read_csv(paste0(data_dir,"meta.csv"))
meta <- readRDS(sample_data) %>%
  mutate(Pellet = ifelse(DeviceID %in% c("Saliva", "Stool"), Main, Pellet)) %>%
  filter(Set %in% c("2", "3", "4", "5", "Stool", "Saliva"), Pellet != "n/a") %>%
  select(DeviceID:Super,Swallow_date:Swallow_time, Recover_date:Recover_time, Hours_in_body, location) %>%
  unique %>%
  mutate(location = gsub("Capsule", "Device", as.character(location)))
### -- Data wrangling

# Filter data by contaminants, reverse, and only identified by site
d %>%
  filter(is.na(`Potential contaminant`)) %>%
  filter(is.na(`Reverse`)) %>%
  filter(is.na(`Only identified by site`)) -> d_filtered

d_filtered %>%
  dplyr::select(contains("Protein IDs"), contains("LFQ")) %>%
  dplyr::select(-"Majority protein IDs") %>%
  column_to_rownames("Protein IDs")-> d_lfq

#d_lfq <- read_tsv("./proteinGroups_lfq.tsv")

# Rename columns

colnames_unique <- data.frame(col_ori = colnames(d_lfq),
                              col_ori_short = str_replace_all(colnames(d_lfq), ".*T760_", ""))

colnames_unique %>%
  arrange(col_ori) %>%
  mutate(colnames_unique = make.unique(str_replace_all(col_ori_short, "_.*", ""))) %>%
  mutate(colnames_simple = str_replace_all(colnames_unique, "\\..*", "")) %>%
  arrange(desc(colnames_unique)) %>%
  distinct(colnames_simple, .keep_all = T) %>%
  pull(col_ori) -> include_samples

d_lfq %>%
  rownames_to_column("Protein") %>%
  pivot_longer(!Protein, names_to = "sample", values_to = "value") %>%
  filter(sample %in% include_samples) %>%
  mutate(sample_short =  str_replace_all(str_replace_all(sample, ".*T760_", ""), "_.*", "")) %>%
  mutate(value = ifelse(value == 0, NA, value)) -> d_red

unique(d_red$sample_short) -> colnames

meta_red <- meta %>%
  filter(Pellet %in% colnames) %>%
  arrange(Pellet) %>%
  column_to_rownames("Pellet") %>%
  mutate(Subject = as.character(Subject),
         Type = gsub("Capsule", "Device", Type))

# Sample exclusion

d_red %>%
  drop_na(value) %>%
  right_join(meta_red %>% rownames_to_column("sample_short")) %>%
  group_by(sample_short, location) %>%
  summarise(n = n()) %>%
  ungroup() -> d_stats

d_stats %>%
  group_by(location) %>%
  summarise(median = median(n)) -> stats_median

stats_sd <- sd(d_stats$n)

stats_include <- d_stats %>%
  filter(n > min(stats_median$median) - 2.5*stats_sd) %>%
  pull(sample_short)

# Wide data format

d_red %>%
  filter(sample_short %in% stats_include) %>%
  drop_na(value) %>%
  mutate(value = log2(value)) %>%
  dplyr::select(-sample) %>%
  spread(sample_short, value) %>%
  column_to_rownames("Protein") -> d_wide

# Proteins quantified in at least 70% of all samples
rowSums(!is.na(d_wide))/length(d_wide) -> NAs
names(NAs[NAs > 0.7]) -> Proteins_0.7

### -- Figures (proteomics)

## Extended Data Figure 7A, quantified proteins
d_stats %>%
  filter(sample_short %in% stats_include) %>%
  left_join(meta_red %>% rownames_to_column("sample_short")) %>%
  group_by(location) %>%
  mutate(`Sample ID` = c(1:length(unique(sample_short)))) %>%
  ungroup() %>%
  mutate(location = ifelse(sample_short %in% stats_include, as.character(location), "excluded"),
         location = factor(location, levels = c("Device", "Stool", "excluded"), labels = c("Devices", "Stool", "excluded"))) -> fig_ed7A_d

hline_df <- data.frame(location = c("Devices", "Stool"), 
                       median = c(round(stats_median$median[1], 0), round(stats_median$median[2], 0)),
                       xlab = c(180, 45), 
                       ylab = c(stats_median$median[1]-600, stats_median$median[2]-600))
(ed_7a <- ggplot(data = fig_ed7A_d %>%
                   filter(location != "excluded"), aes(x = `Sample ID`, y = n, fill = location))+
    geom_point(pch = 21, color = "grey40", size = 3) +
    scale_y_continuous(limits = c(0,10000)) +
    labs(y = "Number of proteins") +
    geom_hline(data = hline_df, aes(yintercept = median), lty = "dotted", lwd = 2, color = "grey20") +
    geom_text(data = hline_df, aes(x = xlab, y = ylab, label = median), fontface = "bold", size = 6) +# geom_hline(yintercept = stats_median$median[2], lty = "dotted", lwd = 1, color = "grey20") +
    scale_fill_manual(values = CapAndStoolColors, guide = "none") +
    facet_grid(.~location, scales = "free_x")) 

ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_7a_subpanel_n_proteins_by_stool_v_devices.pdf"), plot = ed_7a)

## Extended Data Figure 7B, Correlation ID with meta

d_stats %>%
  filter(sample_short %in% stats_include) %>%
  left_join(meta_red %>% rownames_to_column("sample_short")) %>%
  mutate(`Sample ID` = c(1:length(unique(sample_short)))) %>%
  mutate(location = ifelse(sample_short %in% stats_include, location, "excluded"),
         Type = gsub("Device", "", Type)) -> fig_ed7B_d

(ed_7b <- ggplot(data = fig_ed7B_d, aes(x = Type, y = n, group = Type, fill = Type)) +
    geom_boxplot() +
    scale_fill_manual(values = CapTypeAndStoolColors, guide = "none") +
    labs(x = "", y = "Number of proteins") +
    scale_y_continuous(limits = c(0,9000)))

ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_7b_subpanel_proteins_by_location.pdf"), plot = ed_7b)

## Extended Data Figure 7C, Correlation CV with meta
d_red %>%
  filter(sample_short %in% stats_include) %>%
  drop_na(value) %>%
  filter(Protein %in% Proteins_0.7) %>%
  left_join(meta_red %>% rownames_to_column("sample_short")) %>%
  group_by(Protein, Type) %>%
  summarise(`CV [%]` = sd(value) / mean(value) * 100) %>%
  ungroup() %>%
  group_by(Type) %>%
  arrange(`CV [%]`) %>%
  mutate(Rank = c(length(unique(Protein)):1)) %>%
  ungroup() %>%
  mutate(color = ifelse(`CV [%]` <= 100, "2", ifelse(`CV [%]` > 100 & `CV [%]` <= 200, "1", "0")),
         Type = gsub("Device", "", Type)) %>%
  drop_na(color) -> fig_ed7C_d

(ed_7c <- ggplot(data = fig_ed7C_d, aes(x = Type, y = `CV [%]`, group = Type, fill = Type))+
    geom_violin(draw_quantiles = c(0.5))+
    scale_fill_manual(values = CapTypeAndStoolColors, guide = "none") +
    labs(x = ""))

ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_7c_subpanel_proteins_by_location.pdf"), plot = ed_7c)

## Extended Data Figure 7D, Intensities by protein rank
d_red %>%
  filter(sample_short %in% stats_include) %>%
  left_join(meta_red %>% rownames_to_column("sample_short")) %>%
  drop_na(value) %>%
  filter(Protein %in% Proteins_0.7) %>%
  group_by(Protein, location) %>%
  summarise(`LFQ intensity` = mean(value)) %>%
  ungroup()%>%
  group_by(location) %>%
  arrange(`LFQ intensity`) %>%
  mutate(`Protein rank` = c(length(unique(Protein)):1)) %>%
  ungroup() %>%
  mutate(Symbol = mapIds(org.Hs.eg.db, 
                         keys=str_replace_all(Protein, ";.*", ""), 
                         column="SYMBOL", 
                         keytype="UNIPROT",
                         multiVals="first")) %>%
  left_join(d %>%
              select(`Protein IDs`, `Gene names`), by = c("Protein" = "Protein IDs")) %>%
  mutate(name = ifelse(is.na(Symbol) | Symbol == "NA", gsub(';.*', "", `Gene names`), Symbol)) -> fig_ed7D_d

(ed_7d <- ggplot(fig_ed7D_d, aes(x = `Protein rank`, y = log10(`LFQ intensity`))) +
    geom_point(pch = 21, color = "black", size = 3) +
    scale_y_log10() +
    labs(y = expression('log'[10]*' (Intensity)')) +
    geom_text_repel(data = fig_ed7D_d %>% slice_max(`Protein rank`, n = 10),
                    aes(x = `Protein rank`, y = log10(`LFQ intensity`), label = name), 
                    force_pull = 0.6,
                    # direction = "x",
                    force = 0.9,
                    point.padding = 2,
                    nudge_x = -1250
    ) +
    geom_text_repel(data = fig_ed7D_d %>% slice_min(`Protein rank`, n = 12),
                    aes(x = `Protein rank`, y = log10(`LFQ intensity`), label = name),
                    force_pull = 0.8,
                    force = 0.9,
                    point.padding = 4,
                    # nudge_y = 0.0001,
                    nudge_x = 350
    ) +
    # geom_text(data = figD_d %>% slice_min(`Protein rank`, n = 12),
    #           aes(x = 500, y = log10(`LFQ intensity`), label = name),
    #           nudge_y = 0.001
    #                 ) +
    facet_wrap(.~location) + 
    theme_classic() + 
    theme(strip.background = element_blank(),
          legend.text=element_markdown(size=14),
          legend.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 14)))

ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_7d_subpanel_protein_rank_by_intensity.pdf"), plot = ed_7d)


## Figure 4a, Protein intensity in capsule versus stool
d_red %>%
  filter(sample_short %in% stats_include) %>%
  dplyr::select(!sample) %>%
  filter(Protein %in% Proteins_0.7) %>%
  mutate(value = log10(value)) %>%
  left_join(meta_red %>% rownames_to_column("sample_short")) %>%
  group_by(Protein, location) %>%
  summarise(median = median(value, na.rm = T),
            n = n()) %>%
  ungroup() %>%
  drop_na(median) %>%
  pivot_wider(names_from = location, values_from = c(median, n)) %>%
  dplyr::rename(Device = median_Device, Stool = median_Stool) %>%
  mutate(Symbol = mapIds(org.Hs.eg.db, 
                         keys=str_replace_all(Protein, ".*;", ""), 
                         column="SYMBOL", 
                         keytype="UNIPROT",
                         multiVals="first")) %>%
  mutate(diff = Device - Stool) %>%
  left_join(d %>%
              select(`Protein IDs`, `Gene names`), by = c("Protein" = "Protein IDs")) %>%
  mutate(label = ifelse(abs(diff) > 0.5, TRUE, FALSE),
         color = ifelse(diff > 0.5, '#743C2F', ifelse(diff < 0.5, '#2364AA', NA)),
         name = ifelse(is.na(Symbol) | Symbol == "NA", gsub(';.*', "", `Gene names`), Symbol)) -> d_plot_rank_location

(main_4a <- ggplot() +
    geom_point(data = d_plot_rank_location, aes(x = Device, y = Stool), alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, lty = "dotted", size = 1, color = "grey") +
    labs(x = expression('log'[10]*' (abundance in devices)'), y = expression('log'[10]*' (abundance in stool)')) +
    geom_smooth(method = "lm") +
    geom_point(data = (d_plot_rank_location %>% filter(label == TRUE)), aes(x = Device, y = Stool, color = color), size = 3) +
    scale_color_manual(values = rev(CapAndStoolColors), guide = "none") +
    geom_text_repel(data = d_plot_rank_location %>% filter(label == TRUE), aes(x = Device, y = Stool, label = name)))

ggsave(paste0(fig_dir_main_subpanels, "Fig_4a_subpanel_abnd_stool_v_devices.pdf"), plot = main_4a)


## Figure 4d, Pearson correlations
d_red %>%
  filter(sample_short %in% stats_include) %>%
  drop_na(value) %>%
  mutate(value = log2(value)) %>%
  dplyr::select(-sample) %>%
  spread(sample_short, value) %>%
  column_to_rownames("Protein") -> d_wide

pearson <- cor(d_wide, use = "pairwise.complete.obs") 
## set lower triangle (which repeats values from upper triangle) to NA so we can easily remove later
pearson[lower.tri(pearson, diag = T)] <- NA

pearson_df <- pearson %>%
  as.data.frame() %>%
  rownames_to_column("sample_in") %>%
  gather(sample_out, value, !"sample_in") %>%
  left_join(meta_red %>% rownames_to_column("sample_in")) %>%
  dplyr::select(sample_in, sample_out, value, Subject, location) %>%
  dplyr::rename(Subject_in = Subject, Location_in = location) %>%
  left_join(meta_red %>% rownames_to_column("sample_out")) %>%
  dplyr::select(sample_in, sample_out, value, Subject_in, Location_in, Subject, location) %>%
  dplyr::rename(Subject_out = Subject, Location_out = location) %>%
  mutate(within = ifelse(Subject_in == Subject_out, "within Subjects", "between Subjects")) %>%
  filter(Location_out == Location_in, !is.na(value)) %>%
  mutate(subject_pair = ifelse(as.numeric(Subject_in) < as.numeric(Subject_out), paste0(Subject_in, "_", Subject_out), paste0(Subject_out, "_", Subject_in)))

within_subj_pearson <- pearson_df %>%
  filter(within == "within Subjects") %>%
  group_by(within, Location_in, subject_pair) %>%
  summarise(median = median(value),
            n = n()) 

btwn_subj_pearson <- pearson_df %>%
  filter(within == "between Subjects") %>%
  group_by(within, Location_in, subject_pair) %>%
  summarise(median = median(value),
            n=n())

test <- bind_rows(within_subj_pearson, btwn_subj_pearson) %>%
  mutate(Location_in = ifelse(Location_in == "Device", "Devices", Location_in), 
         within = factor(within, levels = c("within Subjects", "between Subjects"), labels = c("Within subjects", "Between subjects"))) 

(main_4d <- ggplot(data = test, aes(x = median, y = Location_in, color = Location_in))+
    geom_boxplot(alpha = 0.7, outlier.shape = NA)+
    geom_jitter(position = position_jitterdodge(jitter.width = 0.4), 
                size = 2, 
                pch = 16,
                alpha = 0.3)+
    facet_grid(within~.)+
    # scale_fill_manual(values = CapAndStoolColors, guide = "none") + 
    scale_color_manual(values = CapAndStoolColors, guide = "none") +
    labs(x = "Median Pearson Correlation", y = "") +
    stat_compare_means(method='wilcox', 
                       comparisons = list(c("Stool", "Devices")),
                       label = "p.signif")) 

ggsave(paste0(fig_dir_main_subpanels, "Fig_4d_subpanel_pearson_correlations_by_location.pdf"), plot = main_4d)

## Figure 4c, PCA stool versus capsule

d_red %>%
  filter(sample_short %in% stats_include) %>%
  filter(Protein %in% Proteins_0.7) %>%
  drop_na(value) %>%
  dplyr::select(-sample) %>%
  spread(sample_short, value) %>%
  column_to_rownames("Protein") -> d_red_0.7

imputePCA(log2(d_red_0.7)) -> d_imputed

meta %>%
  filter(Pellet %in% colnames(d_imputed$completeObs)) %>%
  column_to_rownames("Pellet") %>%
  mutate(location = ifelse(location == "Capsule", "Device", as.character(location)),
         poi = ifelse(Subject == 15 & location == "Stool", "Stool15", ifelse(Subject == 15 & location == "Device", "Device15", "else"))) -> meta_imputed

d_imputed$completeObs[,rownames(meta_imputed)] -> d_imputed_sorted

p <- pca(d_imputed_sorted, meta_imputed)

pca_res <- data.frame(p$rotated) %>%
  rownames_to_column() %>%
  left_join(data.frame(p$metadata) %>%
              rownames_to_column, by = "rowname") %>%
  mutate(poi = factor(poi, levels = c("Device15", "Stool15", "else")))
# (main_4c <- biplot(p,
#              colby = 'location',
#              hline = 0, vline = 0,
#              legendPosition = 'none',
#              lab = NULL,
#              colkey  = CapAndStoolColors,
#              xlab = expression("Principal component 1\n(4.6% variation)"),
#              xlabvjust = 1,
#              ylab = expression("Principal component 2\n(2.7% variation)")))

(main_4c <- ggplot(pca_res,
                   aes(x = PC1, y = PC2, color = location)) + 
    geom_point(size = 2, shape = 16, alpha = 0.75) + 
    scale_color_manual(values = CapAndStoolColors, guide = "none") + 
    labs(x = expression("Principal component 1 (4.6% variation)"),
         y = expression("Principal component 2\n(2.7% variation)")) +
    # labs(x = "Principal component 1 (4.6% variation)", 
    # y = "Principal component 2 (2.7% variation)") +
    # coord_fixed(ratio = 2.66/4.62) + 
    # theme_minimal() + 
    theme(plot.margin = margin(0.5, 0.5, 1, 1, "cm")))

ggsave(paste0(fig_dir_main_subpanels, "Fig_4c_subpanel_pca.pdf"), plot = main_4c)

## Figure 4e, Subject 15

# (main_4e <- biplot(p,
#              colby = 'poi',
#              hline = 0, vline = 0,
#              legendPosition = 'none',
#              lab = NULL,
#              xlab = "Principal component 1 (4.6% variation)",
#              ylab = "Principal component 2 (2.7% variation)",
#              colkey  = c('Device15'="#2364AA", 'Stool15'="#743C2F", 'else'='grey')))

(main_4e <- ggplot(pca_res,
                   aes(x = PC1, y = PC2, color = poi, alpha = poi)) + 
    geom_point(size = 2, shape = 16) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c(CapAndStoolColors, "grey"),  guide = "none") +
    scale_alpha_manual(values = c(0.9, 0.9, 0.6), guide = "none") + 
    labs(x = expression("Principal component 1 (4.6% variation)"),
         y = expression("Principal component 2\n (2.7% variation)"),
         title = "Subject 15") +
    # labs(x = "Principal component 1 (4.6% variation)", 
    # y = "Principal component 2 (2.7% variation)") +
    # coord_fixed(ratio = 2.66/4.62) + 
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = margin(0.5, 0.5, 0.5, 1, "cm")))  

ggsave(paste0(fig_dir_main_subpanels, "Fig_4e_subpanel_pca_subj15.pdf"), plot = main_4e)


## Extended Data Figure 7e, PCA on top 500 proteins
d_red %>%
  filter(sample_short %in% stats_include) %>%
  right_join(meta_imputed %>% rownames_to_column("sample_short")) %>%
  drop_na(value) %>%
  filter(Protein %in% Proteins_0.7) %>%
  group_by(Protein) %>%
  summarise(`LFQ intensity` = mean(value)) %>%
  ungroup()%>%
  arrange(`LFQ intensity`) %>%
  mutate(`Protein rank` = c(length(unique(Protein)):1)) %>%
  ungroup() %>%
  filter(`Protein rank` <= 500) %>%
  pull(Protein) -> Proteins_top500

d_red %>%
  filter(sample_short %in% stats_include) %>%
  filter(Protein %in% Proteins_top500) %>%
  drop_na(value) %>%
  dplyr::select(-sample) %>%
  spread(sample_short, value) %>%
  column_to_rownames("Protein") -> d_red_top500

imputePCA(log2(d_red_top500)) -> d_imputed_top500

meta %>%
  filter(Pellet %in% colnames(d_imputed_top500$completeObs)) %>%
  column_to_rownames("Pellet") -> meta_red_top500

d_imputed_top500$completeObs[,rownames(meta_red_top500)] -> d_imputed_top500_sorted

p_top500 <- pca(d_imputed_top500_sorted, meta_red_top500)

pca_res_top500 <- data.frame(p_top500$rotated) %>%
  rownames_to_column() %>%
  left_join(data.frame(p_top500$metadata) %>%
              rownames_to_column, by = "rowname")

(ed_7e <- ggplot(pca_res_top500,
                 aes(x = PC1, y = PC2, color = location)) + 
    geom_point(size = 2, shape = 16, alpha = 0.75) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = CapAndStoolColors, guide = "none") + 
    labs(x = expression("Principal component 1 (8.02% variation)"),
         y = expression("Principal component 2\n (4.99% variation)"), 
         title = "Top 500 proteins") +
    # labs(x = "Principal component 1 (4.6% variation)", 
    # y = "Principal component 2 (2.7% variation)") +
    # coord_fixed(ratio = 4.99/8.02) + 
    theme(plot.margin = margin(0.5, 0.5, 0.5, 1, "cm"),
          plot.title = element_text(hjust = 0.5, face = "bold")))

ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_7e_pca_top500.pdf"), plot = ed_7e)

# biplot(p_top500,
#        colby = 'location',
#        hline = 0, vline = 0,
#        legendPosition = 'right',
#        lab = NULL,
#        encircle = TRUE,
#        encircleFill = TRUE,
#        colkey  = c('Device'="#440154FF", 'Stool'="#FDE725FF"))

## Extended Data Figure 7f, PCA on not normalized data
d_filtered %>%
  dplyr::select(contains("Protein IDs"), contains("Intensity")) %>%
  dplyr::select(-"Majority protein IDs") %>%
  column_to_rownames("Protein IDs")-> d_int

d_int[,grepl(pattern = "^Intensity..*", colnames(d_int))] -> d_int

colnames_unique_int <- data.frame(col_ori = colnames(d_int),
                                  col_ori_short = str_replace_all(colnames(d_int), ".*T760_", ""))

colnames_unique_int %>%
  arrange(col_ori) %>%
  mutate(colnames_unique = make.unique(str_replace_all(col_ori_short, "_.*", ""))) %>%
  mutate(colnames_simple = str_replace_all(colnames_unique, "\\..*", "")) %>%
  arrange(desc(colnames_unique)) %>%
  distinct(colnames_simple, .keep_all = T) %>%
  pull(col_ori) -> include_samples_int

d_int %>%
  rownames_to_column("Protein") %>%
  pivot_longer(!Protein, names_to = "sample", values_to = "value") %>%
  filter(sample %in% include_samples_int) %>%
  mutate(sample_short =  str_replace_all(str_replace_all(sample, ".*T760_", ""), "_.*", "")) %>%
  mutate(value = ifelse(value == 0, NA, value)) -> d_red_nn

d_red_nn %>%
  filter(sample_short %in% stats_include) %>%
  filter(Protein %in% Proteins_0.7) %>%
  drop_na(value) %>%
  dplyr::select(-sample) %>%
  spread(sample_short, value) %>%
  column_to_rownames("Protein") -> d_red_nn_0.7

imputePCA(log2(d_red_nn_0.7)) -> d_imputed_nn

meta %>%
  filter(Pellet %in% colnames(d_imputed_nn$completeObs)) %>%
  column_to_rownames("Pellet") -> meta_red_nn

d_imputed_nn$completeObs[,rownames(meta_red_nn)] -> d_imputed_nn_sorted

p_nn <- pca(d_imputed_nn_sorted, meta_red_nn)

pca_res_nn <- data.frame(p_nn$rotated) %>%
  rownames_to_column() %>%
  left_join(data.frame(p_nn$metadata) %>%
              rownames_to_column, by = "rowname")

(ed_7f <- ggplot(pca_res_nn,
                 aes(x = PC1, y = PC2, color = location)) + 
    geom_point(size = 2, shape = 16, alpha = 0.75) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = CapAndStoolColors, guide = "none") + 
    labs(x = expression("Principal component 1 (44.69% variation)"),
         y = expression("Principal component 2\n (1.90% variation)"), 
         title = "Not normalized data") +
    # labs(x = "Principal component 1 (4.6% variation)", 
    # y = "Principal component 2 (2.7% variation)") +
    # coord_fixed(ratio = 1.9/44.69) + 
    theme(plot.margin = margin(0.5, 0.5, 0.5, 1, "cm"),
          plot.title = element_text(hjust = 0.5, 
                                    face = "bold")))

ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_7f_pca_not_normalized.pdf"), plot = ed_7f)

# biplot(p_nn,
#        colby = 'location',
#        hline = 0, vline = 0,
#        legendPosition = 'right',
#        lab = NULL,
#        encircle = TRUE,
#        encircleFill = TRUE,
#        colkey  = c('Device'="#440154FF", 'Stool'="#FDE725FF"))
# 
## Figure 4B, Volcano plot
meta_red_top500 %>%
  rownames_to_column("TubeID") -> meta_red_exp

design <- model.matrix(~ as.factor(meta_red_exp[match(colnames(d_red_0.7), meta_red_exp$TubeID),]$location))

fit <- lmFit(log2(d_red_0.7), design)
fit <- eBayes(fit)
d_limma_out<- topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
  rownames_to_column("Uniprot") %>%
  mutate(SYMBOL = mapIds(org.Hs.eg.db, 
                         keys=str_replace(Uniprot, ";.*", ""),
                         column="SYMBOL",
                         keytype="UNIPROT",
                         multiVals="first")) %>%
  mutate(Ensembl = mapIds(org.Hs.eg.db, 
                          keys=str_replace(Uniprot, ";.*", ""),
                          column="ENSEMBL",
                          keytype="UNIPROT",
                          multiVals="first"),
         color = factor(ifelse(adj.P.Val < 0.05 & logFC < -1, "devices", ifelse(adj.P.Val < 0.05 & logFC > 1, "stool", "none")), levels = c("devices", "stool", "none")))

(main_4b <- ggplot() +
    geom_point(data = d_limma_out, aes(x = logFC, y = -log10(adj.P.Val), color = color, alpha = color)) +
    # geom_point(data = d_limma_out %>% filter(P.Value < 0.05, logFC < -1), aes(x = logFC, y = -log10(adj.P.Val)), alpha = 0.9, size = 2, color = "darkblue") +
    # geom_point(data = d_limma_out %>% filter(P.Value < 0.05, logFC > 1), aes(x = logFC, y = -log10(adj.P.Val)), alpha = 0.9, size = 2, color = "brown3") +
    labs(x = expression('log'[2]*' ([stool]/[devices])'),
         y = expression('-log'[10]*" (adjusted "~italic("p")~"-value)")) + 
    scale_color_manual(values = c(CapAndStoolColors,"grey65"), guide = "none") +
    scale_alpha_manual(values = c(0.9, 0.9, 0.5), guide = "none") +
    scale_x_continuous(breaks = c(-5:3)) +
    geom_text_repel(data = d_limma_out %>% filter(adj.P.Val < 0.05, abs(logFC) > 1), aes(x = logFC, y = -log10(adj.P.Val), label = SYMBOL)) +
    geom_hline(yintercept = -log10(0.05), lty = "dotted")+
    geom_vline(xintercept = -1, lty = "dotted")+
    geom_vline(xintercept = 1, lty = "dotted"))

ggsave(paste0(fig_dir_main_subpanels, "Fig_4b_subpanel_diff_abundance_stool_v_devices.pdf"), plot = main_4b)



## Extended Data Figure 7g, Number of significant hits between locations
meta_red_top500 %>%
  mutate(Type = ifelse(location == "Stool", 5, gsub("Capsule ", "", Type))) -> meta_red_stool

for(i in c(1:5)){
  print(i)
  for(j in c(2:5)){
    
    if(j <= i){
      next
    } else{

      print(j)
      
      meta_red_stool %>%
        rownames_to_column("TubeID") %>%
        filter(Type %in% c(i, j)) -> meta_red_exp_type
      
      d_red_0.7[,meta_red_exp_type$TubeID] -> tmp
      
      design <- model.matrix(~ as.factor(meta_red_exp_type[match(colnames(tmp), meta_red_exp_type$TubeID),]$Type))
      
      fit <- lmFit(log2(tmp), design)
      fit <- eBayes(fit)
      d_limma_out <- topTable(fit, number = Inf, confint = TRUE, coef = 2, adjust.method = "fdr") %>%
        rownames_to_column("Uniprot") %>%
        mutate(SYMBOL = mapIds(org.Hs.eg.db, 
                               keys=str_replace(Uniprot, ";.*", ""),
                               column="SYMBOL",
                               keytype="UNIPROT",
                               multiVals="first")) %>%
        mutate(Ensembl = mapIds(org.Hs.eg.db, 
                                keys=str_replace(Uniprot, ";.*", ""),
                                column="ENSEMBL",
                                keytype="UNIPROT",
                                multiVals="first")) %>%
        mutate(sample_in = i, sample_out = j)
      
      if(i ==1 & j == 2){
        limma_total <- d_limma_out
      } else{
        limma_total <- full_join(limma_total, d_limma_out)
      }
    }
  }
}

ed_7g_d <- limma_total %>%
  filter(adj.P.Val < 0.05) %>%
  group_by(sample_in, sample_out) %>%
  summarise(n = n()) 

(ed_7g <- ggplot(ed_7g_d %>% 
                   filter(sample_out == 5) %>%
                   mutate(sample_in = paste0("Device ", sample_in)),
                 aes(x = sample_in, y = sample_out))+
    geom_tile(aes(fill = n)) +
    geom_text(aes(label = round(n, 1))) +
    scale_fill_viridis(option = "viridis", begin = 0.5, guide = guide_colorbar(title.vjust = 0.8), 
                       ) +
    # ylim(4.8, 5.2)
    labs(x = "", 
         y = "Stool",
         fill = "# Significant hits") + 
    coord_equal() +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 14,
                                 angle = 90,
                                 vjust = .5),
      axis.text.y = element_blank(),
      axis.title.y = element_text(size = 14,
                                  angle =0, vjust = 0.5),
      axis.ticks = element_blank(),
      legend.position='bottom',
      legend.text = element_text(size=14),
      legend.title = element_text(size = 14)))

ggsave(paste0(fig_dir_ed_subpanels, "ED_Fig_7g_significantly_different_proteins_stool_v_device_type.pdf"), plot = ed_7g)


## Figure 4f, Sample distances

canberra <- readRDS(paste0(data_dir, "canberra_pairwise_16s.RDS")) 

df_samples <- readRDS(sample_data) %>%
  filter(!drop_16s) %>%
  mutate(Pellet = ifelse(location == "Stool", Main, Pellet)) %>%
  select(Pellet, id_16s, location, Subject)
  
micro_similarity <- as.matrix(canberra) %>%
  data.frame() %>%
  rownames_to_column("V1") %>%
  pivot_longer(cols = starts_with("S"), names_to = "Sample2", values_to = "dist") %>%
  dplyr::rename(Sample1 = V1) %>%
  left_join(df_samples %>%
              select(id_16s, location, Subject, Pellet), by = c("Sample1" = "id_16s")) %>%
  left_join(df_samples %>%
              select(id_16s, location, Subject, Pellet), by = c("Sample2" = "id_16s"), suffix = c("_1", "_2")) %>%
  filter(Sample1 != Sample2, Sample1 < Sample2, location_1 == location_2) %>%
  filter(location_1 != "Saliva") %>%
  mutate(Pellet_1 = as.numeric(Pellet_1),
         Pellet_2 = as.numeric(Pellet_2),
         pair = ifelse(Pellet_1 < Pellet_2, paste0(Pellet_1, "_", Pellet_2), paste0(Pellet_2, "_", Pellet_1)))

protein_similarity <- pearson_df %>%
  select(sample_in:value) %>%
  mutate(sample_in = as.numeric(sample_in), 
         sample_out = as.numeric(sample_out), 
         pair = ifelse(sample_in < sample_out, paste0(sample_in, "_", sample_out), paste0(Pellet_2, "_", sample_out))) %>%
  left_join(micro_similarity, by = "pair") %>%
  filter(!is.na(dist)) # filter out samples without Canberra distance calculated - these are samples with poor 16S sequencing


(main_4f <- ggplot(protein_similarity, 
                   aes(x = dist, y = value, color = location_1)) + 
    geom_point(alpha = 0.6) +
    scale_x_reverse() + 
    scale_color_manual(values = CapAndStoolColors, guide = "none") + 
    geom_smooth(method = "lm", fill = "black", alpha = 1) +
    labs(x = "Microbiota Canberra distance between samples",
         y = "Correlation between host proteomes"))

ggsave(paste0(fig_dir_main_subpanels, "Fig_4f_subpanel_pearson_correlation_by_Canberra_dist.pdf"), plot = main_4f)


## Make final Figure 4

# Make final Extended Data Figure 7