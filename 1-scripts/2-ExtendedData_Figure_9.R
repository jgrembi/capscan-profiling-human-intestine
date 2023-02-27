#######################################
# Capsule Study - Gen1 
#
# Extended Data Figure 9
#
# Author: Jake Folz with minor edits by Jess Grembi
#######################################

rm(list=ls())
# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# Tyrosocholic Acid Head To Tail ####

#MSMS spectra 1 (top)
MSMS1_MetaData_t <- read_xlsx(paste0(data_dir, "AA_conjugated_BAs_15S.xlsx"),sheet=2,col_names="meta data", col_types="text",n_max=9)
MSMS1_Spectrum_t <- read_xlsx(paste0(data_dir, "AA_conjugated_BAs_15S.xlsx"),sheet=2,col_names=c("MSMS1_m/z","MSMS1_intensity"), col_types=c("numeric","numeric"),skip=9)
MSMS1_NormIntensity_t <- MSMS1_Spectrum_t$MSMS1_intensity/(max(MSMS1_Spectrum_t$MSMS1_intensity)*0.01)
MSMS1_Spectrum_t <- cbind(MSMS1_Spectrum_t,MSMS1_NormIntensity_t)

#MSMS spectra 2 (bottom)
# MSMS2_MetaData_t <- read_xlsx(paste0(data_dir, "AA_conjugated_BAs_15S.xlsx"),sheet=3,col_names="meta data", col_types="text",n_max=9)
MSMS2_Spectrum_t <- read_xlsx(paste0(data_dir, "AA_conjugated_BAs_15S.xlsx"),sheet=3,col_names=c("MSMS2_m/z","MSMS2_intensity"), col_types=c("numeric","numeric"),skip=9)
MSMS2_NormIntensity_t <- MSMS2_Spectrum_t$MSMS2_intensity/(max(MSMS2_Spectrum_t$MSMS2_intensity)*0.01)
MSMS2_Spectrum_t <- cbind(MSMS2_Spectrum_t,MSMS2_NormIntensity_t)


#Combined MSMS spectra files
MS2_HtT_t <- bind_rows(MSMS1_Spectrum_t,MSMS2_Spectrum_t)

#Head to Tail MS/MS spectra
(p_t <- ggplot(MS2_HtT_t)+
  geom_col(aes(x=`MSMS1_m/z`,y=`MSMS1_NormIntensity_t`),color="black",width=0.2, position = position_dodge()) + 
  geom_text_repel(aes(x=`MSMS1_m/z`,y=`MSMS1_NormIntensity_t`,label = round(`MSMS1_m/z`,4)),data = MS2_HtT_t[MS2_HtT_t$MSMS1_NormIntensity_t>15,]
                  ,nudge_y = 20,nudge_x=-9, size= 3)+
  geom_col(aes(x=`MSMS2_m/z`,y=-`MSMS2_NormIntensity_t`),color="red",width=0.2, position = position_dodge()) + 
  geom_text_repel(aes(x=`MSMS2_m/z`,y=-`MSMS2_NormIntensity_t`,label = round(`MSMS2_m/z`,4)),data = MS2_HtT_t[MS2_HtT_t$MSMS2_NormIntensity_t>15,]
                  ,nudge_y = -20,nudge_x=-9, size= 3)+
  geom_hline(yintercept = 0, color = "black", size=0.5)+
  geom_vline(xintercept = 55, color = "black", size=0.8)+
  ylab("Relative spectral intensity")+
  xlab("m/z")+
  labs(title = (paste("MS/MS spectrum of",gsub("NAME: ","",MSMS1_MetaData_t[1,1]),"&","Tyrosocholic Acid (library)")),
       caption = "")+
  theme_minimal()+
  scale_x_continuous(limits = c(55,max(MSMS1_Spectrum_t$`MSMS1_m/z`)+10), 
                     expand = c(0,0), 
                     breaks = seq(50,max(MSMS1_Spectrum_t$`MSMS1_m/z`)+10 ,
                                  by = round((max(MSMS1_Spectrum_t$`MSMS1_m/z`)+10-50)/15,-1))) +
  scale_y_continuous(limits=c(-100,100), expand= (c(0,2)), labels=c(100,50,0,50,100))+
  theme(plot.title = element_text(size=11, face= "bold",hjust = 0.5,vjust=2),
        axis.title.x=element_text(face="italic")))

ggsave(paste0(fig_dir_ed_subpanels, 'ED_Fig_9a_TyroCAggsave.pdf'), p_t, width=7, height=4.5)


# Phenylalanocholic acid Acid Head To Tail #### 

#MSMS spectra 1 (top)
MSMS1_MetaData_p <- read_xlsx(paste0(data_dir,"AA_conjugated_BAs_15S.xlsx"),sheet=6,col_names="meta data", col_types="text",n_max=9)
MSMS1_Spectrum_p <- read_xlsx(paste0(data_dir,"AA_conjugated_BAs_15S.xlsx"),sheet=6,col_names=c("MSMS1_m/z","MSMS1_intensity"), col_types=c("numeric","numeric"),skip=9)
MSMS1_NormIntensity_p <- MSMS1_Spectrum_p$MSMS1_intensity/(max(MSMS1_Spectrum_p$MSMS1_intensity)*0.01)
MSMS1_Spectrum_p <- cbind(MSMS1_Spectrum_p,MSMS1_NormIntensity_p)

#MSMS spectra 2 (bottom)
# MSMS2_MetaData <- read_xlsx(paste0(data_dir,"AA_conjugated_BAs_15S.xlsx"),sheet=7,col_names="meta data", col_types="text",n_max=9)
MSMS2_Spectrum_p <- read_xlsx(paste0(data_dir,"AA_conjugated_BAs_15S.xlsx"),sheet=7,col_names=c("MSMS2_m/z","MSMS2_intensity"), col_types=c("numeric","numeric"),skip=9)
MSMS2_NormIntensity_p <- MSMS2_Spectrum_p$MSMS2_intensity/(max(MSMS2_Spectrum_p$MSMS2_intensity)*0.01)
MSMS2_Spectrum_p <- cbind(MSMS2_Spectrum_p,MSMS2_NormIntensity_p)

#Combined MSMS spectra files
MS2_HtT_p <- bind_rows(MSMS1_Spectrum_p,MSMS2_Spectrum_p)

#Head to Tail MS/MS spectra
(p_p <- ggplot(MS2_HtT_p)+
    geom_col(aes(x=`MSMS1_m/z`,y=`MSMS1_NormIntensity_p`),color="black",width=0.2, position = position_dodge()) + 
    geom_text_repel(aes(x=`MSMS1_m/z`,y=`MSMS1_NormIntensity_p`,label = round(`MSMS1_m/z`,4)),data = MS2_HtT_p[MS2_HtT_p$MSMS1_NormIntensity_p>15,]
                    ,nudge_y = 20,nudge_x=-9, size= 3)+
    geom_col(aes(x=`MSMS2_m/z`,y=-`MSMS2_NormIntensity_p`),color="red",width=0.2, position = position_dodge()) + 
    geom_text_repel(aes(x=`MSMS2_m/z`,y=-`MSMS2_NormIntensity_p`,label = round(`MSMS2_m/z`,4)),data = MS2_HtT_p[MS2_HtT_p$MSMS2_NormIntensity_p>15,]
                    ,nudge_y = -20,nudge_x=-9, size= 3)+
    geom_hline(yintercept = 0, color = "black", size=0.5)+
    geom_vline(xintercept = 55, color = "black", size=0.8)+
    ylab("Relative spectral intensity")+
    xlab("m/z")+
    labs(title = (paste("MS/MS spectrum of",gsub("NAME: ","",MSMS1_MetaData_p[1,1]),"&","Phenylalanocholic acid (library)")),
         caption = "")+
    theme_minimal()+
    scale_x_continuous(limits = c(55,max(MSMS1_Spectrum_p$`MSMS1_m/z`)+10), 
                       expand = c(0,0), 
                       breaks = seq(50,max(MSMS1_Spectrum_p$`MSMS1_m/z`)+10 ,
                                    by = round((max(MSMS1_Spectrum_p$`MSMS1_m/z`)+10-50)/15,-1))) +
    scale_y_continuous(limits=c(-100,100), expand= (c(0,2)), labels=c(100,50,0,50,100))+
    theme(plot.title = element_text(size=11, face= "bold",hjust = 0.5,vjust=2),
          axis.title.x=element_text(face="italic")))

ggsave(paste0(fig_dir_ed_subpanels, 'ED_Fig_9b_PheCAggsave.pdf'), p_p, width=7, height=4.5)


# Leucocholic Acid Head To Tail #### 

#MSMS spectra 1 (top)
MSMS1_MetaData_l <- read_xlsx(paste0(data_dir, "AA_conjugated_BAs_15S.xlsx"),sheet=4,col_names="meta data", col_types="text",n_max=9)
MSMS1_Spectrum_l <- read_xlsx(paste0(data_dir, "AA_conjugated_BAs_15S.xlsx"),sheet=4,col_names=c("MSMS1_m/z","MSMS1_intensity"), col_types=c("numeric","numeric"),skip=9)
MSMS1_NormIntensity_l <- MSMS1_Spectrum_l$MSMS1_intensity/(max(MSMS1_Spectrum_l$MSMS1_intensity)*0.01)
MSMS1_Spectrum_l <- cbind(MSMS1_Spectrum_l,MSMS1_NormIntensity_l)

#MSMS spectra 2 (bottom)
# MSMS2_MetaData <- read_xlsx(paste0(data_dir, "AA_conjugated_BAs_15S.xlsx"),sheet=5,col_names="meta data", col_types="text",n_max=9)
MSMS2_Spectrum_l <- read_xlsx(paste0(data_dir, "AA_conjugated_BAs_15S.xlsx"),sheet=5,col_names=c("MSMS2_m/z","MSMS2_intensity"), col_types=c("numeric","numeric"),skip=9)
MSMS2_NormIntensity_l <- MSMS2_Spectrum_l$MSMS2_intensity/(max(MSMS2_Spectrum_l$MSMS2_intensity)*0.01)
MSMS2_Spectrum_l <- cbind(MSMS2_Spectrum_l,MSMS2_NormIntensity_l)

#Combined MSMS spectra files
MS2_HtT_l <- bind_rows(MSMS1_Spectrum_l,MSMS2_Spectrum_l)

#Head to Tail MS/MS spectra
(p_l <- ggplot(MS2_HtT_l)+
    geom_col(aes(x=`MSMS1_m/z`,y=`MSMS1_NormIntensity_l`),color="black",width=0.2, position = position_dodge()) + 
    geom_text_repel(aes(x=`MSMS1_m/z`,y=`MSMS1_NormIntensity_l`,label = round(`MSMS1_m/z`,4)),data = MS2_HtT_l[MS2_HtT_l$MSMS1_NormIntensity_l>15,]
                    ,nudge_y = 20,nudge_x=-9, size= 3)+
    geom_col(aes(x=`MSMS2_m/z`,y=-`MSMS2_NormIntensity_l`),color="red",width=0.2, position = position_dodge()) + 
    geom_text_repel(aes(x=`MSMS2_m/z`,y=-`MSMS2_NormIntensity_l`,label = round(`MSMS2_m/z`,4)),data = MS2_HtT_l[MS2_HtT_l$MSMS2_NormIntensity_l>15,]
                    ,nudge_y = -20,nudge_x=-9, size= 3)+
    geom_hline(yintercept = 0, color = "black", size=0.5)+
    geom_vline(xintercept = 55, color = "black", size=0.8)+
    ylab("Relative spectral intensity")+
    xlab("m/z")+
    labs(title = (paste("MS/MS spectrum of",gsub("NAME: ","",MSMS1_MetaData_l[1,1]),"&","Leucholic acid (library)")),
         caption = "")+
    theme_minimal()+
    scale_x_continuous(limits = c(55,max(MSMS1_Spectrum_l$`MSMS1_m/z`)+10), 
                       expand = c(0,0), 
                       breaks = seq(50,max(MSMS1_Spectrum_l$`MSMS1_m/z`)+10 ,
                                    by = round((max(MSMS1_Spectrum_l$`MSMS1_m/z`)+10-50)/15,-1))) +
    scale_y_continuous(limits=c(-100,100), expand= (c(0,2)), labels=c(100,50,0,50,100))+
    theme(plot.title = element_text(size=11, face= "bold",hjust = 0.5,vjust=2),
          axis.title.x=element_text(face="italic")))

ggsave(paste0(fig_dir_ed_subpanels, 'ED_Fig_9c_LeuCAggsave.pdf'), p_l, width=7, height=4.5)

plot_grid(p_t, p_p, p_l, labels = "auto", label_size = 16, nrow = 3)
ggsave(filename = paste0(fig_dir, "ED_Figure_9.pdf"), width = 16, height = 10)
