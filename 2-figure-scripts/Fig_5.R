#######################################
# Capsule Study - Gen1 
#
# Generate Figure 5 for manuscript
#
# Authors: Becca Culver & Jess Grembi
#######################################

rm(list=ls())
# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# Only look at data of samples from Sets 2,3,4,5 that meet 16S cut-offs & bile acid cut-offs
# ps.bs.raw <- readRDS(phyloseq_bilesalt) %>% subset_samples(., drop_16s == F)
# ps.bs<-subset_samples(ps.bs.raw, Set %in% c('2','3','4','5','Stool'))
# ps.bs@sam_data$Type <- gsub("Capsule", "Device", ps.bs@sam_data$Type)

bs <- read.csv(paste0(data_dir, "bile_salts_wAAconjugates.csv")) %>%
  dplyr::rename(Super = Sample.ID) %>%
  select(-(Subject.Number:Capsule.Type)) %>%
  mutate(Super = as.character(Super))
  
df_samples <- readRDS(paste0(data_dir, "sample_data.RDS")) %>%
  mutate(Type = gsub("Capsule", "Device", Type), 
         Super = ifelse(Type == "Stool", Main, Super))

df_bs <- bs %>%
  left_join(df_samples %>%
              filter(!duplicate_16s, !duplicate_meta) %>%
              select(DeviceID:Super, location, drop_16s), by = "Super") %>%
  unique() %>%
  filter(Set %in% c('2','3','4','5','Stool')) #%>% # Remove Set 1 & additional microbiota reproducibility samples


table(df_bs$location)
table(df_bs$Type)
#######################################
# Figure 5a - Bile acid deconjugation/dehyrdroxylation overview
#######################################

padding=5
a_img <- image_read(paste0(fig_dir_main_subpanels, "Fig_5a_subpanel_ba_overview.png"))
a_ratio <- image_info(a_img)$height/image_info(a_img)$width
(ba.overview <- ggplot() + 
    coord_fixed(a_ratio) + 
    background_image(a_img)) + 
  theme(plot.margin = margin( t = padding, l = padding, r = padding, b = padding, unit = "pt"))

#######################################
# Figure 5b - Concentration of bile acids
#######################################

# Primary: 'Cholic.acid','Chenodeoxycholic.acid'
# Primary Conjugated: 'Glycocholic.acid','Taurocholic.acid','Glycochenodeoxycholic.acid','Taurochenodeoxycholic.acid','Tauro.a.Muricholic.acid','Phenylalanocholic.acid','Leucholic.acid','Tyrosocholic.acid'
# Secondary: 'Deoxycholic.acid','Lithocholic.acid','Ursodeoxycholic.acid'
# Secondary Conjugated:'Glycodeoxycholic.acid','Taurodeoxycholic.acid','Glycolithocholic.acid','Taurolithocholic.acid','Glycoursodeoxycholic.acid','Tauroursodeoxycholic.acid','Glycohyodeoxycholic.acid'

df.ratios <- df_bs %>%
  mutate(Location = ifelse(Type %in% 'Stool','Stool','Device')) %>%
  mutate(primarySum = Cholic.acid + Chenodeoxycholic.acid) %>%
  mutate(priConjSum = Glycocholic.acid+Taurocholic.acid+Glycochenodeoxycholic.acid+Taurochenodeoxycholic.acid+Tauro.a.Muricholic.acid) %>%
  mutate(secondarySum = Deoxycholic.acid+Lithocholic.acid+Ursodeoxycholic.acid) %>%
  mutate(secConjSum = Glycodeoxycholic.acid+Taurodeoxycholic.acid+Glycolithocholic.acid+Taurolithocholic.acid+Glycoursodeoxycholic.acid+Tauroursodeoxycholic.acid+Glycohyodeoxycholic.acid) %>%
  mutate(microConj = Phenylalanocholic.acid+Leucholic.acid+Tyrosocholic.acid+Ala.TriHydroxylated.BA+Arg.Dihydroxylated.BA+Arg.TriHydroxylated.BA+Asn.Dihydroxylated.BA.1+Asn.Dihydroxylated.BA.2+Asn.TriHydroxylated.BA+Cys.Dihydroxylated.BA+Cys.TriHydroxylated.BA+Gln.Dihydroxylated.BA
         +Gln.TriHydroxylated.BA+Glu.Dihydroxylated.BA+Glu.TriHydroxylated.BA+His.TriHydroxylated.BA+Lys.Dihydroxylated.BA+Lys.TriHydroxylated.BA+Met.TriHydroxylated.BA+Phe.Dihydroxylated.BA+Ser.Dihydroxylated.BA+Ser.TriHydroxylated.BA+Trp.TriHydroxylated.BA+Tyr.Dihydroxylated.BA,
         total_abundances =(priConjSum+secConjSum+primarySum+secondarySum+microConj))

df2plot_b <- df.ratios %>%
  mutate(total_bileAcids = log10(total_abundances))

my_comparisons_overall <- list(c('Device','Stool'))
my_comparisons <- list(c('Device 1','Device 4'),c('Device 4','Stool'),c('Device 1','Stool'))

stat_test <- df2plot_b %>%
  wilcox_test(total_bileAcids ~ Type, p.adjust.method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "Type") %>%
  filter((group1 == "Device 1" & group2 %in% c("Device 4", "Stool")) | (group1 == "Device 4" & group2 == "Stool"))

(b <- ggplot(df2plot_b, aes(y=total_bileAcids, x=Type))+
    geom_boxplot(aes(fill=Type), alpha=0.9, outlier.shape = NA)+
    geom_beeswarm(alpha=0.4)+
    scale_fill_manual(values=CapTypeAndStoolColors)+
    stat_pvalue_manual(stat_test, 
                       label = "p.adj.signif",
                       tip.length = 0.01,
                       y.position = c(7.5, 7.7, 7.5), 
                       bracket.shorten = 0.05) +
    labs(y=expression('log'[10]*'(concentration of bile acids)'),x='')+
    theme(legend.position='none'))

ggsave(paste0(fig_dir_main_subpanels,'Fig_5b_subpanel_boxplot_abundance_allTypes.pdf'), b, width=6,height=4)

#######################################
# Figure 5c - Mean realtive abundance of bile acids
#######################################
ba.levels <- c('Cholic.acid','Chenodeoxycholic.acid','Glycocholic.acid','Taurocholic.acid','Glycochenodeoxycholic.acid','Taurochenodeoxycholic.acid','Tauro.a.Muricholic.acid','Deoxycholic.acid','Lithocholic.acid','Ursodeoxycholic.acid','Glycodeoxycholic.acid','Taurodeoxycholic.acid','Glycolithocholic.acid','Taurolithocholic.acid','Glycoursodeoxycholic.acid','Tauroursodeoxycholic.acid','Glycohyodeoxycholic.acid','microConj')

# Colors
pal<-c('#8F2D14','#C43E1C',
       '#A69A11','#CBBC15','#E8D721','#ECDE46','#EFE46B',
       '#5B4C6B','#7A668F','#8E7CA2',
       '#2D625C','#3A7E77','#479A91','#57B2A7','#73BFB6','#8FCCC5','#ABD8D3',
       '#090A0B')

df<- df.ratios %>% 
  melt(id.vars = c('Main','Subject','pH','Type','Set','Location', 'total_abundances'), 
       measure.vars = c("Tauro.a.Muricholic.acid", "Tauroursodeoxycholic.acid","Taurocholic.acid",
                        "Glycoursodeoxycholic.acid","Glycohyodeoxycholic.acid",'Glycocholic.acid',
                        'Taurochenodeoxycholic.acid','Taurodeoxycholic.acid','Cholic.acid','Ursodeoxycholic.acid',
                        'Glycochenodeoxycholic.acid','Glycodeoxycholic.acid','Taurolithocholic.acid',
                        'Chenodeoxycholic.acid','Deoxycholic.acid','Glycolithocholic.acid','Lithocholic.acid',
                        'microConj')) %>%
  #'Tyrosocholic.acid','Phenylalanocholic.acid','Leucholic.acid', 'Ala.TriHydroxylated.BA','Arg.Dihydroxylated.BA','Arg.TriHydroxylated.BA','Asn.Dihydroxylated.BA.1','Asn.Dihydroxylated.BA.2','Asn.TriHydroxylated.BA','Cys.Dihydroxylated.BA','Cys.TriHydroxylated.BA','Gln.Dihydroxylated.BA',
  #'Gln.TriHydroxylated.BA','Glu.Dihydroxylated.BA','Glu.TriHydroxylated.BA','His.TriHydroxylated.BA','Lys.Dihydroxylated.BA','Lys.TriHydroxylated.BA','Met.TriHydroxylated.BA','Phe.Dihydroxylated.BA','Ser.Dihydroxylated.BA','Ser.TriHydroxylated.BA','Trp.TriHydroxylated.BA','Tyr.Dihydroxylated.BA')) %>%
  mutate(variable = as.character(variable)) %>%
  mutate(variable = factor(variable, levels = ba.levels)) %>%
  select(Subject, Location, variable, value)

df2plot_c <- df %>% 
  group_by(Subject, Location, variable) %>%
  summarise(ba.total = sum(value)) %>%
  ungroup() %>%
  group_by(Subject, Location) %>%
  mutate(ba.overall.total = sum(ba.total, na.rm=TRUE),
         frequency = ba.total/ba.overall.total,
         Subject = factor(Subject, levels=c('1','2','3','4','5','6','7','8','9','11','12','13','14','10','15')))

(c_pie <- ggplot(df2plot_c, aes(x = factor(1), y = frequency, fill = variable)) +
  geom_bar(stat = "identity", width = 1) +
  facet_wrap(~Location + Subject, nrow=2)+
  coord_polar(theta = "y") +
  scale_fill_manual(values = pal) +
  theme_void() +
  theme(legend.position='none'))

ggsave(paste0(fig_dir_main_subpanels, 'Fig_5c_subpanel_subj_pie_charts.pdf'), c_pie, width=8)


#######################################
# Figure 5d - Percent of liver-conjugated bile acids
#######################################

df2plot_d <- df.ratios %>%
  mutate(total_conjugatedBA = log10(priConjSum+secConjSum),
         totalBA = log10(primarySum+secondarySum+priConjSum+secConjSum+microConj),
         pct_conjugatedBA = total_conjugatedBA/totalBA*100)

my_comparisons <- list(c('Device 1','Device 4'),c('Device 4','Stool'),c('Device 1','Stool'))

stat_test_conjugatedBA <- df2plot_d %>%
  wilcox_test(pct_conjugatedBA ~ Type, p.adjust.method = "bonferroni") %>%
  add_significance() %>%
  add_xy_position(x = "Type") %>%
  filter((group1 == "Device 1" & group2 %in% c("Device 4", "Stool")) | (group1 == "Device 4" & group2 == "Stool"))



## This plot shows the ratio of conjugated/all bile acids
## I think this one more clearly shows that:
## 1) Type 1 Devices have most of the liver-conjugated bile acids
## 2) Most of the deconjugation is occurring between Devices 2-4 
(d <- ggplot(df2plot_d, aes(y=pct_conjugatedBA, x=Type))+
    geom_boxplot(aes(fill=Type), outlier.shape = NA, alpha = 0.9)+
    geom_beeswarm(alpha=0.4)+
    scale_fill_manual(values=CapTypeAndStoolColors) +
    stat_pvalue_manual(stat_test_conjugatedBA, 
                       label = "p.adj.signif",
                       tip.length = 0.01,
                       y.position = c(102, 105, 102), 
                       bracket.shorten = 0.05) +
    labs(y=expression('% of liver-conjguated bile acids', x='')) +
    theme(legend.position='none',
          axis.text = element_text(size = 10),
          axis.title.x=element_blank()))

ggsave(paste0(fig_dir_main_subpanels,'Fig_5d_subpanel_boxplot_pct_conj_bileacids.pdf'), d, width = 4, height = 3)



#######################################
# Figure 5e - Barplots of each sample including antibiotic and concentration metadata
#######################################
highlilght<-c('#8F2D14')

df.relabs<-melt(df.ratios, id.vars = c('Main','Subject','pH','Type','Set','Location','total_abundances'), 
                measure.vars = c('Tauro.a.Muricholic.acid', 'Tauroursodeoxycholic.acid','Taurocholic.acid',
                                 'Glycoursodeoxycholic.acid','Glycohyodeoxycholic.acid','Glycocholic.acid',
                                 'Taurochenodeoxycholic.acid','Taurodeoxycholic.acid','Cholic.acid','Ursodeoxycholic.acid',
                                 'Glycochenodeoxycholic.acid','Glycodeoxycholic.acid','Taurolithocholic.acid',
                                 'Chenodeoxycholic.acid','Deoxycholic.acid','Glycolithocholic.acid','Lithocholic.acid',
                                 'microConj')) %>%
  #'Tyrosocholic.acid','Phenylalanocholic.acid','Leucholic.acid','Ala.TriHydroxylated.BA','Arg.Dihydroxylated.BA','Arg.TriHydroxylated.BA','Asn.Dihydroxylated.BA.1','Asn.Dihydroxylated.BA.2','Asn.TriHydroxylated.BA','Cys.Dihydroxylated.BA','Cys.TriHydroxylated.BA','Gln.Dihydroxylated.BA','Gln.TriHydroxylated.BA','Glu.Dihydroxylated.BA','Glu.TriHydroxylated.BA','His.TriHydroxylated.BA','Lys.Dihydroxylated.BA','Lys.TriHydroxylated.BA','Met.TriHydroxylated.BA','Phe.Dihydroxylated.BA','Ser.Dihydroxylated.BA','Ser.TriHydroxylated.BA','Trp.TriHydroxylated.BA','Tyr.Dihydroxylated.BA')) %>%
  mutate(variable = as.character(variable)) %>%
  mutate(variable = factor(variable, levels=ba.levels),
         rel_abund = value/total_abundances)

# Determine order of samples based on Type & Subj 10 & 15 last
order <- df.ratios %>%
  select(Main, Type, Subject, Location) %>%
  mutate(subj_reorder = factor(Subject, levels=c('1','2','3','4','5','6','7','8','9','11','12','13','14','10','15')))%>%
  arrange(Location,Type,subj_reorder)

order<-as.character(order$Main)

df2plot_e <- df.relabs %>%
  mutate(Main = fct_relevel(Main, order)) %>%
  mutate(abundance_fill = ifelse(total_abundances > 450000, 0,1),
         Antibiotics_Taken = ifelse(Subject %in% c('10','15'), 1, 0))


(p1<-ggplot(df2plot_e, aes(x=Main, y=rel_abund, fill=variable))+
    geom_bar(stat='identity') +
    theme_minimal() +
    scale_fill_manual(labels = paste("<span style='color:",
                                     pal,
                                     "'>",
                                     levels(df2plot_e$variable),
                                     '</span>'), values = pal)+
    scale_x_discrete(limits = rev) +
    coord_flip() +
    theme(legend.position = 'top',
          legend.box = 'vertical',
          legend.margin = margin(),
          legend.spacing.x = unit(1.0, 'pt'),
          legend.text=element_markdown(size=10))+
    guides(color = guide_legend(override.aes = list(shape = NA, fill=NA)))) ## This line removes the point on the legend for the color

# Using the cowplot package
pdf(paste0(fig_dir_main_subpanels,'Fig_5e_subpanel_barplot_relab_legend.pdf'),width=12,height=6)
legend <- cowplot::get_legend(p1)
grid.draw(legend)
dev.off()

#Replot to exclude legend
(p1<-ggplot(df2plot_e, aes(x=Main, y=rel_abund, fill=variable))+
    geom_bar(stat='identity') +
    scale_fill_manual(values = pal) +
    # scale_x_discrete(limits = rev) +
    # coord_flip()+
    theme(legend.position='none',
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    ylab('Relative abundance of bile acid \n'))

(p2<-ggplot(df2plot_e, aes(x=Main, y=total_abundances, fill=abundance_fill)) +
    geom_bar(stat='identity')+
    scale_fill_gradient(labels = paste("<span style='color:",
                                       c(highlilght, 'black'),
                                       "'>",
                                       df2plot_e$abundance_fill,
                                       '</span>'),
                        low=highlilght, high='black')+
    scale_x_discrete(limits = rev) + 
    coord_flip() +
    ylim(0,50000000)+ #Remove outliers for visual
    ylab('Concentration \n')+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = 'none'))


(p3<-ggplot(df2plot_e, aes(x=Main, fill=Antibiotics_Taken)) +
    geom_bar() +
    scale_fill_gradient(low='white', high=highlilght) +
    theme(legend.position='none') +
    scale_x_discrete(limits = rev) +
    coord_flip() +
    ylab('Recent \n antibiotics') +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))


(p4<-ggplot(df2plot_e, aes(x=Main, fill=Type)) +
    geom_bar() +
    scale_fill_manual(values=CapTypeAndStoolColors) +
    theme(legend.position='none')+
    scale_x_discrete(limits = rev) +
    coord_flip() +
    ylab('Sample \n Type') +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()))


# All plots together
(bar.plots<-plot_grid(p4,p1,p2,p3, nrow=1, rel_widths=c(1,6,1,1)))

## The panels for this figure were assembled in Illustrator due to the intricate placement, especially for panel 5e.

