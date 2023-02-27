#######################################
# Capsule Study - Gen1 
#
# Extended Data Figure 2
#
# Author: Jess Grembi
#######################################

rm(list=ls())
# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

#######################################
# ED Fig. 2a - Schematic overview generated in adobe illustrator
#######################################
a_img <- image_read(paste0(fig_dir_ed_subpanels, "ED_Fig_2a_sample_flowchart.png")) 
a_ratio <- image_info(a_img)$height/image_info(a_img)$width
(a <- ggplot() + 
    coord_fixed(a_ratio) + 
    background_image(a_img)) + 
  theme(plot.margin = margin(t = fig_padding, l = fig_padding, r = fig_padding, b = fig_padding/3, unit = "pt"))

#######################################
# ED Fig. 2b - Transit time of devices by subject
#######################################

#Read in sample data
df_samples <- readRDS(sample_data) %>%
  dplyr::rename(transit_time = Hours_in_body) %>%
  mutate(Subject = factor(Subject)) %>%
  filter(Subject != "Control", location == "Capsule", !is.na(Time_of_day), Set %in% c("2","3","4","5"), duplicate_16s == F, duplicate_meta == F)


df_samples_plot <- df_samples %>%
  mutate(Subject = fct_reorder(Subject, transit_time, .fun = 'median', na.rm = T),
         Time_of_day = factor(Time_of_day, levels = c("After lunch", "After dinner")))


## Summary statistics of transit time for capsules in Sets 2-5, taken 3h post-meal as directed.
## We do not include capsules from Set 1 because these were taken at any time of day that participants desired, including with a meal.
## Therefore, these samples are different from those evaluated in the main study in Sets 2-5.
summary(df_samples_plot$transit_time)
cols <- scales::viridis_pal(end = 0.6)(length(unique(df_samples_plot$Time_of_day)))

(transit_time_bySubj <- ggplot(df_samples_plot, 
                               aes(x = Subject, y = transit_time)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(alpha = 0.5, width = 0.3, aes(color = Time_of_day), size = 2.5) + 
    ylim(5,70) +
    labs(x = "Subject number", y = "Device gut transit time (h)", color = "Time device swallowed") + 
    theme(legend.position = "top",
          legend.spacing.x = unit(1.0, 'pt')) + 
    # scale_color_viridis_d(end = 0.6) + 
    scale_color_manual(labels = paste("<span style='color:",
                                      cols,
                                      "'>",
                                      levels(df_samples_plot$Time_of_day),
                                      "</span>"),
                       values = cols) +
    guides(color = guide_legend(nrow = 1,
                                override.aes = list(shape = NA),
                                keywidth = unit(0, "pt"))))

ggsave(filename = paste0(fig_dir_ed_subpanels, "ED_Fig_2b_subject_transit_time.pdf"), plot = transit_time_bySubj)

#######################################
# ED Fig. 2c,d - Transit time of devices by diet and time
#######################################

#Read in dietary recall data
diet.df <- read_excel(paste0(data_dir, "FoodQuestionnaire.xlsx")) %>%
  select(Subject, Capsule_set, TimeOfDay_capsule, Food_pre_post,Dairy:Coffee) %>%
  mutate(Coffee = as.numeric(ifelse(Coffee == "1 (decaf)", 0.5, Coffee)), 
         Subject = factor(Subject), 
         Capsule_set = as.character(Capsule_set)) %>%
  dplyr::rename(Grains = Bread_cereal_rice_pasta, Set = Capsule_set)


df <- df_samples %>%
  left_join(diet.df, by = c("Subject", "Set")) %>%
  pivot_longer(cols = Dairy:Coffee, names_to = "Food_type", values_to = "Est_servings") %>%
  mutate(Est_servings = ifelse(is.na(Est_servings), 0, Est_servings),
         Food_01 = factor(ifelse(Est_servings > 0, 1, 0), labels = c("No", "Yes")),
         Food_type = gsub("_", "/", ifelse(Food_type == "Nuts_Beans_Eggs", "Nuts_beans_eggs", Food_type)), 
         Food_type_sig_pre = ifelse(Food_type %in% c("Dairy", "Fruit", "Meat/fish","Vegetables"), "sig", "ns"),
         Food_type_sig_post = ifelse(Food_type %in% c("Coffee", "Meat/fish", "Nuts/beans/eggs", "Vegetables"), "ns", "sig"))


(pre_capsule_foods <- ggplot(df %>% 
                               filter(!location %in% c("Stool", "Saliva", "Control"), 
                                      !is.na(Food_pre_post), 
                                      Set %in% c("2","3","4","5"),
                                      Food_pre_post == "Pre-capsule"), 
                             aes(x = Food_01, y = transit_time)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes(alpha = Food_type_sig_pre, color = Time_of_day), width = 0.2, shape = 21) + 
    scale_alpha_manual(values = c(0.4, 1), guide = "none") +
    scale_color_manual(labels = paste("<span style='color:",
                                      cols,
                                      "'>",
                                      levels(df_samples_plot$Time_of_day),
                                      "</span>"),
                       values = cols) +
    stat_compare_means(comparisons = list(c('No','Yes')), 
                       method = "wilcox.test",
                       label = "p.signif",
                       label.y = 65) + 
    ylim(5,70) +
    facet_wrap(~Food_type, nrow = 2) + 
    labs(x = "Food reported eaten in meal prior to device ingestion", y = "Device gut transit time (h)") + 
    theme(strip.text.x = element_text(size = 14, face = "bold"),
          legend.text=element_markdown(size=14),
          legend.position = "none") + 
    guides(color = guide_legend(nrow = 1,
                                override.aes = list(shape = NA),
                                keywidth = unit(0, "pt"))))


ggsave(filename = paste0(fig_dir_ed_subpanels, "ED_Fig_2c_transit_time_by_pre_foods.pdf"), plot = pre_capsule_foods)

(post_capsule_foods <- ggplot(df %>% 
                                filter(!location %in% c("Stool", "Saliva", "Control"), 
                                       !is.na(Food_pre_post), 
                                       Set %in% c("2","3","4","5"),
                                       Food_pre_post == "Post-capsule"), 
                              aes(x = Food_01, y = transit_time)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(aes(alpha = Food_type_sig_post, color = Time_of_day), width = 0.2, shape = 21) + 
    scale_alpha_manual(values = c(0.4, 1), guide = "none") +
    stat_compare_means(comparisons = list(c('No','Yes')), 
                       method = "wilcox.test",
                       label = "p.signif",
                       label.y = 65) +
    scale_color_manual(labels = paste("<span style='color:",
                                      cols,
                                      "'>",
                                      levels(df_samples_plot$Time_of_day),
                                      "</span>"),
                       values = cols) +
    ylim(5,70) +
    facet_wrap(~Food_type, nrow = 2) + 
    labs(x = "Food reported eaten in meal after device ingestion", y = "Device gut transit time (h)") + 
    theme(strip.text.x = element_text(size = 14, face = "bold"),
          legend.text=element_markdown(size=14), 
          legend.position = "none") +
    guides(color = guide_legend(nrow = 1,
                                override.aes = list(shape = NA),
                                keywidth = unit(0, "pt"))))


ggsave(filename = paste0(fig_dir_ed_subpanels, "ED_Fig_2d_transit_time_by_post_foods.pdf"), plot = post_capsule_foods)


plot_grid(a, transit_time_bySubj, pre_capsule_foods, post_capsule_foods, labels = "auto", label_size = 16, nrow = 3, rel_heights = c(1, 1.15, 1.15))
ggsave(filename = paste0(fig_dir, "ED_Figure_2.pdf"), width = 16, height = 14)
