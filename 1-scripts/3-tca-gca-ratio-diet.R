#######################################
# Capsule Study - Gen1 
#
# Diet vs Bile acids
#  
#
# Author: Original analysis and code written by Becca Culver, updated by Jess Grembi 
#######################################


rm(list=ls())
# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))


# Read in metabolomics data
bs <- read.csv(paste0(data_dir, "bile_salts_wAAconjugates.csv")) %>%
  dplyr::rename(Super = Sample.ID) 

# Number of samples with bile acid data
table(bs$SampleType)

# Read in sample data
df_samples <- readRDS(paste0(data_dir, "sample_data.RDS")) %>%
  mutate(Type = gsub("Capsule", "Device", Type), 
         Super = ifelse(Type == "Stool", Main, Super)) %>%
  filter(!duplicate_16s, !duplicate_meta, location != "Saliva") %>%
  select(DeviceID:Super, location, drop_16s)

# Read in diet data
diet.df <- readxl::read_xlsx(paste0(data_dir, "FoodQuestionnaire.xlsx")) %>%
  select(Subject, Capsule_set, TimeOfDay_capsule, Food_pre_post,Dairy:Coffee) %>%
  mutate(Coffee = as.numeric(ifelse(Coffee == "1 (decaf)", 0.5, Coffee)), 
         Subject = factor(Subject), 
         Capsule_set = as.character(Capsule_set)) %>%
  dplyr::rename(Grains = Bread_cereal_rice_pasta, Set = Capsule_set)

# Only look at data of samples from Sets 2,3,4,5 that meet bile acid cut-offs
df <- bs %>%
  select(-(Subject.Number:Capsule.Type)) %>%
  mutate(Super = as.character(Super)) %>%
  left_join(df_samples, by = "Super") %>%
  unique() %>%
  filter(Set %in% c('2','3','4','5')) %>% # Remove Set 1 & additional microbiota reproducibility samples
  left_join(diet.df, by = c("Subject", "Set")) %>%
  pivot_longer(cols = Dairy:Coffee, names_to = "Food_type", values_to = "Est_servings") %>%
  mutate(Est_servings = ifelse(is.na(Est_servings), 0, Est_servings),
         Food_01 = factor(ifelse(Est_servings > 0, 1, 0), labels = c("No", "Yes")),
         Food_type = gsub("_", "/", ifelse(Food_type == "Nuts_Beans_Eggs", "Nuts_beans_eggs", Food_type)),
         Food_type_time = paste0(Food_type,'_',Food_pre_post,'uni',Food_01),
         base=paste0(Food_type,'_',Food_pre_post),
         Food_type_01 = paste0(Food_type, '_',Food_01))

# Create a data frame that contains the ratios for gca and tca abundances
# GCA and TCA only
df.ratio <- df %>%
  mutate(gca_ratio = Glycocholic.acid/(Glycocholic.acid+Taurocholic.acid),
         tca_ratio = Taurocholic.acid/(Glycocholic.acid+Taurocholic.acid))
# can also define ratio ias tca_ratio = Taurocholic.acid/Glycocholic.acid and all statistical associations remain significant

# Plot distributions (not normal)
ggplot(df.ratio, aes(gca_ratio)) +
  geom_histogram()

ggplot(df.ratio, aes(tca_ratio)) +
  geom_histogram()

# Use a corrected, wilcox test to evaluate whether the ratios differed across food types

# Checked both though they should be identical
#TCA ratio
stat_test <- df.ratio %>%
  wilcox_test(tca_ratio ~ Food_type_01, p.adjust.method = "bonferroni") %>%
  mutate(base1 = str_split(group1, '_', simplify=TRUE)[,1],
         base2 = str_split(group2, '_', simplify=TRUE)[,1]) %>%
  filter(base1==base2) %>%
  add_significance() #%>%
  #add_xy_position(x = "Food_type_time") %>%
  # filter(!p.adj.signif=='ns')

stat_test


#GCA ratio - IDENTICAL TO TCA RATIO
#stat_test <- df.ratio %>%
#    wilcox_test(gca_ratio ~ Food_type_01, p.adjust.method = "bonferroni") %>%
#    mutate(base1 = str_split(group1, '_', simplify=TRUE)[,1],
#        base2 = str_split(group2, '_', simplify=TRUE)[,1]) %>%
#    filter(base1==base2) %>%
#    add_significance() %>%
#add_xy_position(x = "Food_type_time") %>%
#    filter(!p.adj.signif=='ns')


#NOTE: the stat_test does not output which group is in the largest abundance, plot to see when ratio is higher
### Sanity check plots; figures not in manuscript - pvalues in figures are not corrected

ggplot(df.ratio %>% filter(Food_type %in% c('Dairy', 'Vegetables')), aes(x=Food_01, y=tca_ratio)) +
  geom_boxplot() +
  facet_wrap(~Food_type) 

# Looking at absolute concentrations
ggplot(df %>% filter(Food_type %in% c('Dairy', 'Vegetables')),
       aes(y=log10(Glycocholic.acid), x=Food_01)) +
  geom_boxplot() +
  theme(axis.text.x=element_text(angle = 45)) +
  facet_wrap(~Food_type)

ggplot(df %>% filter(Food_type %in% c('Dairy', 'Vegetables')),
       aes(y=log10(Taurocholic.acid), x=Food_01)) +
  geom_boxplot() +
  theme(axis.text.x=element_text(angle = 45)) +
  facet_wrap(~Food_type)
