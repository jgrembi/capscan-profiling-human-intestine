#######################################
# Capsule Study - Gen1 
#
# Generate Figure 1 for manuscript
# Capsule and study overview, pH, pcoa, and differential abundance
#
# Author: Jess Grembi and Rebecca Culver
#######################################


rm(list=ls())
# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

#######################################
# Figure 1A - capsule diagram and indended opening locations 
#######################################
# Read in .jpg of GI & capsules overview
a_img <- image_read(paste0(fig_dir, "gi_overview.png"))
a_ratio <- image_info(a_img)$height/image_info(a_img)$width
(a <- ggplot() + 
    coord_fixed(a_ratio) + 
    background_image(a_img)) + 
  theme(plot.margin = margin(t = fig_padding, l = fig_padding, r = fig_padding, b = fig_padding/3, unit = "pt"))
# (a <- ggdraw() + 
#     draw_image(paste0(fig_dir, "gi_overview.jpg"), height = 1.5, width = 1,) + #valign = -1.15 
#     theme(plot.margin = unit(c(5,5,5,5), "points")))

#######################################
# Figure 1B - study design
#######################################
# Read in .jpg of study overview
b_img <- image_read(paste0(fig_dir, "study_overview.png"))
b_ratio <- image_info(b_img)$height/image_info(b_img)$width
(b <- ggplot() + 
    coord_fixed(b_ratio) + 
    background_image(b_img)) + 
  theme(plot.margin = margin(t = fig_padding, l = fig_padding, r = fig_padding, b = fig_padding/3, unit = "pt"))

# (b <- ggdraw() + 
#     draw_image(paste0(fig_dir, "study_overview.jpg"), height = 1.5, width = 1) + #valign = -2 
#     theme(plot.margin = unit(c(5,5,5,5), "points")))

#######################################
# Figure 1C - pH across capsule types
#######################################

## Filter out Capsules not in Sets 2-5, and any sample without adequate 16s reads
#  This results in high-quality samples from Sets 2-5, Saliva, and Stool
#  We purposely filter out Set 1 capsules because these capsules were intended as safety tests to ensure passage through each participant.
#  Thus, they were taken at any time of day and without respect to timing before or after meals.
#  Because of the capsule design, disintegration of the capsule is pH driven, so if the capsule is taken *with* food, 
#  the location of sampling is not always as expected.  So, for analyses moving forward, we remove them.

df_samples <- readRDS(sample_data) %>%
  filter(!drop_16s & Set %in% c("2", "3", "4", "5", "Stool", "Saliva"))

# Here's the breakdown of where the remaining 297 samples are from:
table(df_samples$Type)

## First, let's look at the pH of the 210 capsule samples 
df_cap_samples <- df_samples %>%
  data.frame %>%
  filter(location == "Capsule") %>%
  mutate(Type = droplevels(Type))

summary(df_cap_samples$pH)

stat_test <- df_cap_samples %>%
  wilcox_test(pH ~ Type, p.adjust.method="bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "Type")

(c <- ggplot(df_cap_samples, 
             aes(x = Type, y = pH)) + 
    geom_boxplot(aes(fill = Type), outlier.shape = NA) +
    # geom_violin(alpha = 0.9) +
    geom_beeswarm(alpha = 0.4) +
    # geom_jitter(alpha = 0.9) +
    stat_pvalue_manual(stat_test, 
                       label = "p.adj.signif",
                       tip.length = 0.01,
                       y.position = c(8.65, 9,9.35,8.65, 9.7, 8.65), 
                       bracket.shorten = 0.05) +
    labs(x = "", color = "") + 
    scale_fill_manual(values = CapType) + 
    theme(legend.position = "none",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          plot.margin = margin(t = fig_padding, l = fig_padding, r = fig_padding, b = fig_padding, unit = "pt")))

ggsave(filename = paste0(fig_dir, "capsule_pH.pdf"), plot = c, height = 6, width = 6)

#######################################
# Figure 1D - Ordination of all capsule, stool, and saliva samples
#######################################
(ps <- readRDS(clean_phyloseq_object) %>%
   subset_samples(!drop_16s & Set %in% c("2", "3", "4", "5", "Stool", "Saliva")) %>%
   filter_taxa(., function(x) sum(x > 3) > (0.05*length(x)), TRUE) %>% # Gets rid of all taxa not found at a count of 3 in at least 5% of samples (that's 14 samples)
   transform_sample_counts(function(x) {log2(x + 1)}))

# We have a pretty large spread in library sizes/sequencing depth across samples:
summary(ps@sam_data$reads_16s)


# So, first we have to normalize the counts to library size.
# Some *easy* methods that others have done in the past include: 
#   1) rarefying the data by sub-sampling to the same sequencing depth
#   2) converting the data all to proportions
#  BUT, these methods toss really important information about the sample such as the uncertainty (sampling variance). 
#  So we prefer to normalize by sampling depth, which is easy to accomplish using variance stabilization transformation 
#  from the DESeq2 R package!  Here we use a custom function phyloseqTransformation() written by Pratheepa Jeganathan 
#   (now professor of Mathematics and Statistics at McMaster University)
#  The function description is included in the 0.base-functions.R file.

#ps_transform <- phyloseqTransformation(ps)
# ps_transform <- readRDS(clean_phyloseq_object_transformed_filtered)

# We choose canberra distance because the capsule samples are often dominated by one highly abundant ASV.  
# Canberra distance downweights these highly abundant taxa (as does our asinh-transformation above) so 
# that we're comparing the whole community and distances aren't driven primarily by the most abundant taxa
pcoa_canberra <- ordinate(ps,  method = "MDS", distance = "canberra")

var_exp <- get_evals(pcoa_canberra)$variance_exp
scores <- get_scores(pcoa_canberra, sample_data(ps))

(d <- ggplot(scores %>%
               mutate(subj_id = ifelse(Subject == 15, "15", ifelse(Subject == 10, "10", "All other"))),  
             aes(x = -PC1, y = PC2, color = Type, shape = subj_id)) +
    geom_point(size = 3, alpha = 0.8) + 
    ## This code changes the color of the legend text so that it will mirror the points
    scale_color_manual(labels = paste("<span style='color:",
                                     c("black", CapTypeAndStoolColors),
                                     "'>",
                                     levels(ps@sam_data$Type),
                                     "</span>"),
                      values = c("black", CapTypeAndStoolColors), 
                      guide = "none") +    
    coord_fixed(var_exp[2]/var_exp[1]) +
    # scale_color_manual(values = c(CapTypeAndStoolColors, "black")) + 
    scale_shape_manual(values = c(15, 17, 1)) +
    labs(x = paste0("Principal coordinate 1 [", round(var_exp[1],1), "% variance]"),
         y = paste0("Principal coordinate 2 [", round(var_exp[2], 1), "% variance]"),
         shape = "Subject", 
         color = "Sample type") + 
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.position = "top",
          legend.box = "vertical",
          legend.margin = margin(),
          legend.spacing.x = unit(1.0, 'pt'),
          ## This line is needed with the above (lines 151-157) to change the color of the legend text to match the points
          legend.text=element_markdown(size=14),
          legend.title = element_text(size = 14, vjust = 0.66),
          plot.margin = margin(t = fig_padding/5, l = fig_padding/5, r = fig_padding/5, b = fig_padding/5, unit = "pt")) + 
    guides(color = guide_legend(nrow = 1,
                                override.aes = list(shape = NA),
                                keywidth = unit(0, "pt")))) ## This line removes the point on the legend for the color

ggsave(filename = paste0(fig_dir, "pcoa_canberra.pdf"), plot = d, height = 6, width = 8)


## PERMANOVA to determine if saliva is significantly different
dist <- distance(ps, method = "canberra", type = "samples")
adonis(dist~location, data = ps@sam_data %>% data.frame)


#######################################
# Figure 1E - Differential abundance analysis
#######################################
## First look for unique ASVs and their relative abundance and get basic statistics
(ps_ra <- readRDS(clean_phyloseq_object) %>%
    subset_samples(!drop_16s & Set %in% c("2", "3", "4", "5", "Stool", "Saliva")) %>%
    transform_sample_counts(function(x) {x/ sum(x)}))# %>%

table(ps_ra@sam_data$location)

df_ps <- ps_ra %>%
  psmelt() 


###################################
#### ASVs unique to a location
###################################
df_ps_summary <- df_ps %>%
  group_by(Subject, location, ASV) %>%
  summarise(n = n(), 
            n_present = sum(Abundance > 0),
            mean_abnd = mean(Abundance),
            median_abnd = median(Abundance)) %>%
  pivot_wider(id_cols = c(Subject, ASV), names_from = location, values_from = c(n, n_present, mean_abnd, median_abnd)) %>%
  filter(n_present_Stool != 0 | n_present_Capsule != 0 | n_present_Saliva != 0)

### Number of unique ASVs per subject detected in any sample type
subj_n_asv <- df_ps %>%
  group_by(Subject) %>%
  filter(Abundance > 0) %>%
  select(Subject, ASV) %>%
  unique() %>%
  summarise(n_unique_asv = n())

### Number of ASVs per subject detected only in capsules
(capsule_only <- df_ps_summary %>%
    filter(n_present_Stool == 0 & n_present_Saliva == 0) %>%
    group_by(Subject) %>%
    summarise(n = length(unique(ASV)),
              max_median_abnd_Capsule = max(median_abnd_Capsule),
              max_mean_abnd_Capsule = max(mean_abnd_Capsule)) %>%
    left_join(subj_n_asv, by = "Subject") %>%
    mutate(pct_unique = n/ n_unique_asv))


capsule_only %>%
  ungroup() %>%
  summarise(mean_pct = mean(pct_unique)*100,
            sd_pct = sd(pct_unique)*100,
            mean_n = mean(n),
            sd_n = sd(n),
            mean_unique = mean(n_unique_asv),
            sd_unique = sd(n_unique_asv))

## For each subject, between 57-283 ASVs were detected in capsule, but not stool or saliva samples, but 
##the mean relative abundance never exceeded 6% of any capsule sample and median never exceeded 0.4%, 
# so these don't represent major constituents of the microbial community.


### Number of ASVs per subject detected only in stool
(stool_only <- df_ps_summary %>%
    filter(n_present_Capsule == 0 & n_present_Saliva == 0) %>%
    group_by(Subject) %>%
    summarise(n = length(unique(ASV)),
              max_median_abnd_stool = max(median_abnd_Stool),
              max_mean_abnd_stool = max(mean_abnd_Stool)) %>%
    left_join(subj_n_asv, by = "Subject") %>%
    mutate(pct_unique = n/ n_unique_asv))

stool_only %>%
  ungroup() %>%
  summarise(mean_pct = mean(pct_unique)*100,
            sd_pct = sd(pct_unique)*100,
            mean_n = mean(n),
            sd_n = sd(n),
            mean_unique = mean(n_unique_asv),
            sd_unique = sd(n_unique_asv))

## For each subject, 5-184 ASVs were detected in stool but not capsule or saliva samples.
## For these, the mean and median relative abundance never exceeded 0.6%, 
## except for Subj 3 who was dominated by E. coli (the mean & median relative abundance never exceeded 8%)


### Number of ASVs per subject detected only in saliva
(saliva_only <- df_ps_summary %>%
    filter(n_present_Capsule == 0 & n_present_Stool == 0) %>%
    group_by(Subject) %>%
    summarise(n = length(unique(ASV)),
              max_median_abnd_saliva = max(median_abnd_Saliva),
              max_mean_abnd_saliva = max(mean_abnd_Saliva)) %>%
    left_join(subj_n_asv, by = "Subject") %>%
    mutate(pct_unique = n/ n_unique_asv))

saliva_only %>%
  ungroup() %>%
  summarise(mean_pct = mean(pct_unique)*100,
            sd_pct = sd(pct_unique)*100,
            mean_n = mean(n),
            sd_n = sd(n),
            mean_unique = mean(n_unique_asv),
            sd_unique = sd(n_unique_asv))

## For each subject, 102-210 ASVs were detected in saliva but not capsule or stool samples.
## For these, the maximum relative abundance from any subject ranged from 2-29% 
## Thus, saliva samples harbor a community very distinct from intestine and stool samples.

###################################
#### Differential abundance between capsule and stool
#### This code is based on code published from:
#### Grembi, J.A., Nguyen, L.H., Haggerty, T.D. et al. 
#### Gut microbiota plasticity is correlated with sustained weight loss on a low-carb or 
#### low-fat dietary intervention. Sci Rep 10, 1405 (2020). 
#### https://doi.org/10.1038/s41598-020-58000-y
###################################

(ps_da <- readRDS(clean_phyloseq_object) %>%
   subset_samples(!drop_16s & Set %in% c("2", "3", "4", "5", "Stool")))
(ps_da <- subset_taxa(ps_da, taxa_sums(ps_da) > 0))


## Estimate size factors
dds <- phyloseq_to_deseq2(ps_da, design = ~ 1)
dds <- estimateSizeFactors(dds, type = "poscount")
size_factors <- sizeFactors(dds)
summary(size_factors)

# Since we do not want to artificially inflate the counts in voom transformation
# we normalize the size factors so the the smallest is equal to 1.
size_factors <- size_factors/min(size_factors)
summary(size_factors)


# Filter sequences
minSubjPrev <- .2*length(unique(ps_da@sam_data$Subject))
# Filter out rare taxa not present in at least 20% of subjects (that's 3 subjects - i.e. it needs to be generalizable to at least 3 people in our cohort)
(pstips <- filter_subject_prevalence(ps_da, thresh = minSubjPrev))
## This leaves us with 671 taxa that are present in at least 2 subjects

seqtab <- t(as(pstips@otu_table, "matrix"))
smpdata <- data.frame(pstips@sam_data) %>%
  mutate(location = factor(location, levels = c("Stool", "Capsule")))

taxtable <- data.frame(pstips@tax_table)  %>%
  rownames_to_column("SeqName") %>%
  mutate(OrgName = SeqName)

#####################################
## Differential abundance with voom
#####################################
# We apply limma + voom method on data for participants in this study.
# We choose limma because our datasets has replicates, and limma is 
# the only differential abundance testing package which can estimate
# within subject random effect. Also, we have more than a dozen
# samples available and there is no need for very sensitive methods
# as DESeq2 designed for datasets with very few number of samples.
# Using limma+voom was actually the advice we received from a developer of DESeq2,
# Michael Love.

dsgn <- ~ location

limmaRes <- limma_fit(
  seqtab, smpdata, dsgn, 
  sizefac = size_factors[rownames(smpdata)], 
  block = smpdata$Subject, 
  taxtable = taxtable,
  alpha = 1)  %>%
  left_join(data.frame(pstips@tax_table) %>%
              rownames_to_column("SeqName"))

alpha <- 0.05
limmaRes_sig <- limmaRes %>%
  filter(adj.P.Val < alpha, abs(logFC) > 0.75) %>%
  mutate(OrgName = paste0(ifelse(is.na(Genus), Family, Genus), " ", Species, " (", SeqName, ")"),
         col = factor(ifelse(logFC > 0.75, "Increased in capsule", "Increased in stool")))

(diff_abnd <- ggplot(limmaRes_sig, 
       aes(y = fct_reorder(OrgName, logFC, max), x = logFC, fill = col)) + 
  geom_col() + 
  labs(x = expression('log'[2]*' fold change'), y = "", fill = "") + 
  scale_fill_manual(labels = paste("<span style='color:",
                                   CapAndStoolColors,
                                   "'>",
                                   levels(limmaRes_sig$col),
                                   "</span>"),
                    values = CapAndStoolColors) + 
  guides(fill = guide_legend(nrow = 1,
                             reverse = T,
                             override.aes = list(shape = NA, fill=NA))) + 
  theme(legend.position = "bottom",
        legend.box.margin = margin(-2, unit = "pt"),
        legend.margin = margin(t=-2, unit = "pt"),
        axis.text.y = element_text(size = 13),
        # legend.margin = margin(l = 10, unit = "pt"),
        plot.margin = margin(r = -0.8, unit = "pt")))


### This looks for the prevalence of all differentially abundant ASVs in both capsule and stool samples
asv_labels = limmaRes_sig %>%
  ungroup() %>%
  arrange((logFC)) %>%
  pull(OrgName)
  # mutate(OrgName = fct_reorder(OrgName, logFC, max))
prev.cap <- data.frame(df_ps) %>%
  filter(ASV %in% limmaRes_sig$ASV, location %in% c("Capsule","Stool")) %>%
  mutate(location = factor(location, levels = c("Capsule", "Stool"))) %>%
  group_by(ASV, location) %>%
  summarise(n = n(), 
            n_samp = sum(Abundance > 0),
            prev_pct = n_samp/n) %>%
  select(location, ASV, prev_pct) %>%
  filter(location %in% 'Capsule')

prev.stool <- data.frame(df_ps) %>%
  filter(ASV %in% limmaRes_sig$ASV, location %in% c("Capsule","Stool")) %>%
  mutate(location = factor(location, levels = c("Capsule", "Stool"))) %>%
  group_by(ASV, location) %>%
  summarise(n = n(), 
            n_samp = sum(Abundance > 0),
            prev_pct = n_samp/n) %>%
  select(location, ASV, prev_pct) %>%
  filter(location %in% 'Stool')

prev<-rbind(prev.cap,prev.stool) %>%
  left_join(.,limmaRes_sig,by='ASV') %>%
  mutate(location = factor(location, levels = c("Capsule", "Stool"), labels = c("C", "S")),
         OrgName = factor(OrgName, levels = asv_labels))

(prev_hm <- ggplot(prev, aes(x = location, y = OrgName, fill = log10(prev_pct*100))) +
  geom_tile(color = "white") +
  labs(x = "", y = "", fill = "Prev (%)") +
  coord_equal() + 
  theme(legend.position='right',
        legend.margin = margin(t = -2, b = -2, unit = "pt"),
        axis.text.y = element_blank(),
        plot.margin=margin(l=0, r = -1, unit="pt")) +
  scale_fill_gradient(low = "light blue", high = "dark blue", 
                      breaks = c(1, 1.25, 1.5, 1.75, 2), 
                      labels = c(10, 18, 32, 56, 100)))
## Comparing mean of log10 abundance
ggplot(df_ps %>%
         filter(ASV %in% limmaRes_sig$ASV, location %in% c("Capsule", "Stool")), 
       aes(x = log10(Abundance*100))) + 
  geom_histogram(bins = 25) + 
  facet_wrap(ASV~location, scales= "free")


#relative abundance of ASVs
rel_abnd_heatmap <- df_ps %>%
  filter(ASV %in% limmaRes_sig$ASV, location %in% c("Capsule", "Stool")) %>%
  group_by(ASV, location) %>%
  summarise(mean_relAbnd = mean(Abundance),
            median_relAbnd = median(Abundance)) %>%
  ungroup() %>%
  left_join(limmaRes_sig %>%
              select(ASV, OrgName)) %>%
  mutate(location = factor(location, levels = c("Capsule", "Stool"), labels = c("C", "S")),
         OrgName = factor(OrgName, levels = asv_labels))

(relabnd_hm <- ggplot(rel_abnd_heatmap, aes(x = location, y = OrgName, fill = log10(mean_relAbnd*100))) +
  geom_tile(color = "white") +
  labs(x = "", y = "", fill = str_wrap('Relative abundance (%)', width = 10)) +
  coord_equal() + 
  theme(legend.position='right',
        legend.margin = margin(t = -2, b = -2, unit = "pt"),
        axis.text.y = element_blank(),
        plot.margin=margin(l=-1, r = 0, unit="pt")) +
  scale_fill_gradient(low = "yellow", high = "red"))

(prev_abnd_hm <- plot_grid(prev_hm + theme(legend.position = "none", 
                                            plot.margin = margin(l = -1, r = -1, unit = "pt")), 
                            relabnd_hm + theme(legend.position = "none", 
                                               plot.margin = margin(l = -1, r = -1, unit = "pt")), nrow = 1, ncol = 2, rel_widths = c(1,1), align = "vh", axis = "l"))
(legend_plot <- plot_grid(get_legend(prev_hm), 
                     get_legend(relabnd_hm), ncol = 1, axis = "b", align = "v"))
(heatmaps <- plot_grid(prev_abnd_hm, legend_plot, nrow = 1, axis = "l"))
(e <- plot_grid(diff_abnd, 
                heatmaps,
          nrow = 1, rel_widths = c(7, 2), axis = "t"))

ggsave(paste0(fig_dir, "diff_abundance.pdf"), height = 5.8)



### Plot final figure for manuscript
 
lay <- rbind(c(1,2),c(1,3), c(4,4),c(4,4),c(5,5), c(5,5))

margin <- theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

g<-list(a,b,c,d,e)
p.final<-arrangeGrob(grobs=lapply(g,"+",margin), layout_matrix = lay)

p <- as_ggplot(p.final) + # transform to a ggplot
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 18,x = c(0, 0.5, 0.5, 0, 0), y = c(1,1,0.75,0.75,0.25))
ggsave(filename = paste0(fig_dir, "Figure_1.pdf"),p, width = 10, height = 14)
 

