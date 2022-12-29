#######################################
# Supplemental figure: all bar plots with relabs
#######################################

rm(list=ls())
# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

ps.bs.raw <- readRDS(phyloseq_bilesalt) %>% subset_samples(., drop_16s == F)
ps.bs<-subset_samples(ps.bs.raw, !drop_16s & Set %in% c("2", "3", "4", "5", "Stool")) %>%
    filter_taxa(., function(x) sum(x > 3) > (0.05*length(x)), TRUE) %>% # Gets rid of all taxa not found at a count of 3 in at least 5% of samples (that's 14 samples)
    transform_sample_counts(.,function(x) x/sum(x))
                            

                            
df<-psmelt(ps.bs) %>%
    mutate(title=paste0('Subj ',Subject,',\nSet ',Set)) %>%  
    mutate(Type2 = ifelse(Type == 'Capsule 1','D1',
                        ifelse(Type == 'Capsule 2','D2',
                              ifelse(Type == 'Capsule 3',' D3',
                                    ifelse(Type == 'Capsule 4','D4', 'S'))))) %>%
    mutate(name = paste0(Main,'',Type2)) %>%
    mutate(Genus = factor(ifelse(is.na(Genus), ifelse(!is.na(Family), paste0("(Unclassified ", Family, ")"), "(Unclassified family)"), Genus))) %>%
    mutate(Family = factor(ifelse(is.na(Family), ifelse(!is.na(Phylum), paste0("(Unclassified ", Phylum, ")"), "(Unclassified phylum)"), Family)))
              
# Find genera who have max abundance < 5% and categorize as 'Other'
low_abundance <- df %>%
    group_by(Family) %>%
    summarise(max_abund=max(Abundance)) %>%
    filter(max_abund < 0.05)

df.order <- df %>% 
    arrange(Type)

df <- df %>% 
    mutate(Genus = ifelse(Genus %in% low_abundance$Genus, paste0('LowAbund Genera ',Family), as.character(Genus))) %>%
    mutate(Main = factor(Main, levels=unique(df.order$Main))) %>%
    mutate(Location = ifelse(Type %in% 'Stool','Stool','Devices')) %>%
    mutate(Main2 = ifelse(Location %in% 'Devices', paste0(Type2, ', S',Set), Main))

d3 <- df %>% filter(Phylum %in% 'Firmicutes') #Blue
d3.cols<-colorRampPalette(c('#052C39','#D9F3FC'))(length(unique(d3$Family)))

d4 <- df %>% filter(Phylum %in% 'Proteobacteria') #Red
d4.cols<-colorRampPalette(c('#771B18','#F1BDBB'))(length(unique(d4$Family)))

d2 <- df %>% filter(Phylum %in% 'Bacteroidetes') #Purple
d2.cols<-colorRampPalette(c('#371B35','#EBD6E9'))(length(unique(d2$Family)))

d1 <- df %>% filter(Phylum %in% 'Actinobacteria') #Yellow
d1.cols<-colorRampPalette(c('#645002','#FEF2C3'))(length(unique(d1$Family)))

d5 <- df %>% filter(Phylum %in% 'Verrucomicrobia') #Pink
d5.cols<-colorRampPalette(c('#280004','#FFC2C8'))(length(unique(d5$Family)))

main.cols<-data.frame(c(d1.cols,d2.cols,d3.cols,d4.cols,d5.cols))
colnames(main.cols)<-'colors'

genus.order<-df %>% arrange(Phylum, Family) %>% select(Family) %>%
unique() %>% cbind(main.cols) #%>% arrange(Genus)

df2plot <- df %>%
    group_by(Main, Family, Location, Subject) %>%
    dplyr::summarise(Abundance = sum(Abundance)) %>%
    mutate(Family = factor(Family, levels=genus.order$Family))


ggplot(df2plot, aes(x=Main, y=Abundance, fill=Family)) +
    geom_bar(stat='identity', width=1) +
    scale_fill_manual(values=genus.order$colors)+
    facet_wrap(~Location+Subject, scale='free_x', nrow=2) +
    theme_minimal() +
    theme(axis.text.x = element_blank())
ggsave(paste0(fig_dir, 'Figure_S_familyLevel_barplots.pdf'),width=12, height=5)
