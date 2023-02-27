#######################################
# Capsule Study - Gen1 
#
# Extended Data Figure 1
#
# Author: Jess Grembi
#######################################

rm(list=ls())
# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

### Get phyloseq object, subset to only blooming samples with read counts over 2000 (one sample had 2465 so we included it even though our cutoff for other analyses was 2500)
(ps <- readRDS(clean_phyloseq_object) %>%
    subset_samples(!is.na(bloomingID) & reads_16s > 2000))  
    

df_ps <- ps %>%
  #Gets rid of all taxa not found at a count of 15 in at least 1% of samples (that's 1 sample)
  filter_taxa(., function(x) sum(x) > 15, TRUE) %>% # ## 99% of reads from each sample were retained using this filtering
  transform_sample_counts(., function(x) x / sum(x)) %>%
  psmelt() %>%
  select(-(DeviceID:reads_meta_filtered), -(location:experiment)) %>%
  mutate(bloomingID = gsub("bloom.", "", bloomingID)) %>%
  separate(bloomingID, into = c("Type","h_incubate"), sep = "[.]") %>%
  mutate(h_incubate = factor(h_incubate, levels = c("32h", "58h", "87h")), 
         Type = paste0("Device ", Type),
         Genus = factor(ifelse(is.na(Genus), ifelse(!is.na(Family), paste0("(Unclassified ", Family, ")"), "(Unclassified family)"), Genus)))

genus_list <- levels(df_ps$Genus)

df_ps <- df_ps %>%
  mutate(Genus = factor(Genus, levels = c(genus_list[4:29], genus_list[1:3])))

label.df <- data.frame(h_incubate = c("58h", "87h"), Type = c("Device 2", "Device 3"), label = rep("Library prep failed", 2))
(blooming <- ggplot(df_ps,
                      # mutate(Genus = fct_explicit_na(Genus, na_level = "(Unclassified)")), 
                    aes(x = h_incubate, y = Abundance)) + 
    geom_col(aes(fill = Genus)) + 
    geom_text(data = label.df,
              aes(label = label, x = h_incubate, y = 0.25, angle = 90), size = 4) +
    facet_wrap(~Type,  nrow = 1) + 
    labs(y = "Relative abundance", x = "Sampling timepoint") + 
    theme(strip.text.x = element_text(size = 14),
          legend.position = "right") +  
    scale_fill_manual(values = c(
      "#FF7F00", 
      "green4",
      "#E31A1C", # red
      "#6A3D9A", # purple
      "dodgerblue2",
      "coral", # orange
      "black", "gold1",
      "skyblue2", "#FB9A99", # lt pink
      "palegreen2",
      "deeppink1", 
      "#FDBF6F", # lt orange
      "gray70", "khaki2",
      "maroon", "orchid1", "#CAB2D6", # lt purple
      "blue1", "steelblue4",
      "darkturquoise", "green1", "yellow4", "yellow3",
      "darkorange4", "brown", "gray90", "gray20", "gray50"
    )))

ggsave(filename = paste0(fig_dir, "ED_Figure_1.pdf"), blooming, width = 14)

########################
## Below are additional analyses that produces the values reported in this section of the manuscript. 
## None of the plots below appear in the manuscript as main, extended data, or supplementary figures.
########################

##############################################
# Coefficient of variation for relative abundance from 32h to 58h or 87h
##############################################

top_30 <- df_ps %>%
  group_by(Type, h_incubate) %>%
  slice_max(order_by = Abundance, n = 30) %>%
  arrange(-Abundance) %>%
  mutate(rank = rank(-Abundance))

ranks_32h <- df_ps %>%
  filter(h_incubate == "32h") %>%
  group_by(Type) %>%
  mutate(rank = rank(-Abundance, ties.method = "first")) %>%
  mutate(rank = ifelse(Type == "Capsule 1" & rank > 64, NA, 
                       ifelse(Type == "Capsule 2" & rank > 94, NA, 
                              ifelse(Type == "Capsule 3" & rank > 77, NA, 
                                     ifelse(Type == "Capsule 4" & rank > 92, NA, rank)))))

nodetect_32h <- df_ps %>%
  filter(h_incubate == "32h", Abundance == 0) %>%
  mutate(no_32h_detect = "no 32h detect")
  
detect_after32h <- df_ps %>%
  filter(h_incubate != "32h", Abundance > 0) %>%
  mutate(detect_after32 = "detect after 32h")


top_changes <- df_ps %>%
  left_join(ranks_32h %>%
              select(Type, ASV, rank), by = c("Type", "ASV")) %>%
  left_join(nodetect_32h %>%
              select(Type, ASV, no_32h_detect), by = c("Type", "ASV")) %>%
  left_join(detect_after32h %>%
              select(Type, ASV, detect_after32), by = c("Type", "ASV")) %>%
  unique() %>%
  arrange(Type, h_incubate, -Abundance)

(plot <- ggplot(top_changes %>%
                  group_by(Type) %>%
                  filter(!is.na(rank) | (!is.na(detect_after32) & Abundance > 0.1)) %>%
                  ungroup(),
                aes(y = tidytext::reorder_within(ASV, -rank, within = Type), 
                    x = Abundance, fill = Genus)) + 
    geom_col() + 
    facet_grid(Type~h_incubate, scales = "free_y", drop=TRUE) + 
    labs(y = "ASV", x = "Relative abundance") +
    # scale_color_manual(limits = c(NA, "no 32h detect"), values = c(NA, "red")) +
    theme(axis.text.y = element_blank(),
          legend.position = "right") +
    scale_fill_manual(values = c(
      "#FF7F00", 
      "green4",
      "#E31A1C", # red
      "#6A3D9A", # purple
      "dodgerblue2",
      "coral", # orange
      "black", "gold1",
      "skyblue2", "#FB9A99", # lt pink
      "palegreen2",
      "deeppink1", 
      "#FDBF6F", # lt orange
      "gray70", "khaki2",
      "maroon", "orchid1", "#CAB2D6", # lt purple
      "blue1", "steelblue4",
      "darkturquoise", "green1", "yellow4", "yellow3",
      "darkorange4", "brown", "gray90", "gray20", "gray50"
    )))

##############################################
### Change in rank of the top 30 ASVs from 32h to 58h and 87h 
##############################################
df_ranks <- ps %>%
  transform_sample_counts(., rank) %>%
  psmelt %>%
  select(-(DeviceID:reads_meta_filtered), -(location:experiment)) %>%
  mutate(bloomingID = gsub("bloom.", "", bloomingID)) %>%
  separate(bloomingID, into = c("Type","h_incubate"), sep = "[.]") %>%
  mutate(h_incubate = factor(h_incubate, levels = c("32h", "58h", "87h")), 
         Type = paste0("Capsule ", Type)) %>%
  group_by(Type, h_incubate) %>%
  arrange(-Abundance) %>%
  mutate(rank = rank(-Abundance))



df_top30ranks <- df_ranks %>%
  filter(h_incubate == "32h") %>%
  group_by(Type) %>%
  slice_max(order_by = Abundance, n = 30) %>%
  mutate(top30 = "top 30 ASVs")

df_ranks_top30 <- df_ranks %>%
  left_join(df_top30ranks %>% select(Type, ASV, top30)) %>%
  filter(!is.na(top30)) %>%
  pivot_wider(id_cols = c(Type, ASV, Kingdom:Species), names_from = h_incubate, values_from = rank) %>%
  mutate(s_m = abs(`32h`- `58h`),
         s_e = abs(`32h` - `87h`))


df_ranks_top30 %>% 
  group_by(Type) %>% 
  summarise(median_sm = median(s_m),
            max_sm = max(s_m),
            median_se = median(s_e),
            max_se = max(s_e))

##############################################
### New ASVs not detected at 32h but detected 
### at 58h, 87h, or both 
##############################################
ps_asinh = ps %>%
    transform_sample_counts(., function(x) {asinh(x)})

asinh_sums <- sample_sums(ps_asinh)
df_asinh <- ps_asinh %>%
  psmelt() %>%
  select(-(DeviceID:reads_meta_filtered), -(location:experiment)) %>%
  left_join(data.frame(Sample = names(asinh_sums), asinh_sum = asinh_sums), by = "Sample") %>%
  mutate(bloomingID = gsub("bloom.", "", bloomingID),
         rel_abnd = Abundance/asinh_sum) %>%
  separate(bloomingID, into = c("Type","h_incubate"), sep = "[.]") %>%
  mutate(h_incubate = factor(h_incubate, levels = c("32h", "58h", "87h")),
         Type = paste0("Capsule ", Type))

asinh_sums <- df_asinh %>%
  group_by(Sample) %>%
  summarise(total_asinh_abnd = sum(Abundance, na.rm = T))

new_ASVs <- df_asinh %>%
  left_join(asinh_sums %>%
              select(Sample, total_asinh_abnd), by = "Sample") %>%
  mutate(rel_abnd = Abundance/ total_asinh_abnd) %>%
  pivot_wider(id_cols = c(Type, ASV, Kingdom:Species), names_from = h_incubate, values_from = rel_abnd) %>%
  select(Type:`32h`, `58h`,`87h`) %>%
  filter(`32h` == 0 & (`58h` > 0 | `87h` > 0))


df_ranks_top30 %>% 
  group_by(Type) %>% 
  summarise(mean = mean(s_m),
            median = median(s_m),
            max = max(s_m))

## Summary statistics presented in the manuscript
new_ASVs %>%
  group_by(Type) %>%
  summarise(n = n(),
            max_abnd_58 = max(`58h`, na.rm = T),
            max_abnd_87h = max(`87h`, na.rm = T),
            sum_58 = sum(`58h`,na.rm = T),
            sum_87 = sum(`87h`, na.rm = T))
