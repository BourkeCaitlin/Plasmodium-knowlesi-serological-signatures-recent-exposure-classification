library(SeroTrackR)
library(patchwork)
library(tidyverse)


macaque_p1 = SeroTrackR::runPlasmoPipeline(raw_data = here::here("rawData/sabah-macaque/data/plate1/20260302_plate1_Pkclass_CBMLSS.csv"),
                                           platform = "magpix",
                                           panel = here::here("rawData/sabah-macaque/data/panel_macaque.csv"),
                                           plate_layout = here::here("rawData/sabah-macaque/data/plate1/20260302_map-P1_Pkclass_CBMLSS.xlsx"),
                                           std_point = 5,
                                           experiment_name = "macaque-p1", classify =F)

macaque_p2 = runPlasmoPipeline(raw_data = here::here("rawData/sabah-macaque/data/plate2/20230303-PkClass-plate2-MLSS.csv"),
                               platform = "magpix",
                               panel = here::here("rawData/sabah-macaque/data/panel_macaque.csv"),
                               plate_layout = here::here("rawData/sabah-macaque/data/plate2/20260303_P2_Pkclass_MLSS.xlsx"),
                               std_point = 5,
                               experiment_name = "macaque-p2", classify =F)

macaque_p3 = runPlasmoPipeline(raw_data = here::here("rawData/sabah-macaque/data/plate3/20260303-SecondaryDilution_Plate3_MLSS.csv"),
                               platform = "magpix",
                               panel = here::here("rawData/sabah-macaque/data/panel_macaque.csv"),
                               plate_layout = here::here("rawData/sabah-macaque/data/plate3/Plate 3 Secondary dilution 3Mar26.xlsx"),
                               std_point = 5,
                               experiment_name = "macaque-p3", classify =F)


all_mfi_rau = macaque_p1$mfi_outputs$MFI_RAU_long %>% 
  bind_rows(macaque_p2$mfi_outputs$MFI_RAU_long)%>% 
  filter(!str_detect(SampleID, "-400")) %>% filter(!SampleID%in%c("NC", "Blank1")) %>% 
  filter(!str_detect(SampleID, "NC"))

plasma_dil =macaque_p1$mfi_outputs$MFI_RAU_long

secondary_ab = macaque_p3$mfi_outputs$MFI_RAU_long


sample_ids = all_mfi_rau %>% select(SampleID) %>% distinct() 


## sample info
updated_macaque <- readxl::read_excel(here::here("sample-info/Macaque Malaria Samples_master_PCRresults amended Jul26-2.xlsx")) 


new_data <- readxl::read_excel(here::here("sample-info/Pan-Pk-macaque-results-Oct2026.xlsx")) 

macaque <- new_data %>% select(animalid, summary_PCR) %>% 
  left_join(updated_macaque) %>% 
  filter(animalid%in%sample_ids$SampleID) %>% 
  mutate(common_name = tolower(common_name)) %>% 
  dplyr::rename(SampleID = animalid)


# date range

macaque %>% 
  summarise(min = min(collection_date),
            max = max(collection_date))

# how many adults
macaque %>% 
  janitor::tabyl(age)

## macaque species

macaque %>% 
  janitor::tabyl(common_name)%>% 
  mutate(percent = round(percent*100, 1))

# macaque by sex
macaque %>% 
  janitor::tabyl(sex)%>% 
  mutate(percent = round(percent*100, 1))

# macaque by their plasmodium and pk pcr results

macaque %>% 
  janitor::tabyl(summary_PCR)

macaque %>% 
  janitor::tabyl(summary_PCR)


# combine PCR and serology

all_mfi_rau %>% 
  filter(RAU_Method=="Adj_loglog"| (Species=="Pk" & RAU_Method=="loglog")|(Species=="Pf" & RAU_Method=="loglog")|(Antigens=="PvMSP10" & RAU_Method=="loglog")) %>% 
  left_join(macaque) %>%drop_na(common_name) %>% 
  mutate(summary_PCR = factor(summary_PCR, levels = c("Pk", "Plasmodium", "negative"))) %>% 
  select( summary_PCR, collection_date, age, SampleID, common_name) %>% 
  distinct() %>% select(-SampleID) %>% 
  gtsummary::tbl_summary(by = "summary_PCR",
                         digits = list(gtsummary::all_categorical() ~ c(0, 1))) %>% 
  gtsummary::as_flex_table()


# table split by juvenile vs adult macaque

all_mfi_rau %>% 
  filter(RAU_Method=="Adj_loglog"| (Species=="Pk" & RAU_Method=="loglog")|(Species=="Pf" & RAU_Method=="loglog")|(Antigens=="PvMSP10" & RAU_Method=="loglog")) %>% 
  left_join(macaque) %>%drop_na(common_name) %>% 
  mutate(summary_PCR = factor(summary_PCR, levels = c("Pk", "Plasmodium", "negative"))) %>% 
  select( summary_PCR, collection_date, age, SampleID, common_name, sex) %>% 
  distinct() %>% select(-SampleID) %>% 
  gtsummary::tbl_summary(by = "age",
                         digits = list(gtsummary::all_categorical() ~ c(0, 1))) %>% 
  gtsummary::as_flex_table()



# table split by species of macaque

all_mfi_rau %>% 
  filter(RAU_Method=="Adj_loglog"| (Species=="Pk" & RAU_Method=="loglog")|(Species=="Pf" & RAU_Method=="loglog")|(Antigens=="PvMSP10" & RAU_Method=="loglog")) %>% 
  left_join(macaque) %>%drop_na(common_name) %>% 
  select(summary_PCR, SampleID, common_name, sex, age) %>% distinct() %>% select(-SampleID) %>% 
  gtsummary::tbl_summary() %>% 
  gtsummary::as_flex_table()


# classification of macaque data with three modles ------------------------

to_classify = all_mfi_rau %>% 
  filter(RAU_Method=="Adj_loglog"| (Species=="Pk" & RAU_Method=="loglog")) %>% 
  dplyr::rename(mfi = MFI, rau = RAU) %>% 
  select(-c(Species, RAU_Method, mfi)) %>% 
  pivot_wider(names_from = Antigens, values_from = rau) %>% 
  select(sample_id =SampleID, PkMSP10, `PvFAMA-95` = `Pv-fam-a`, `PvMSP1-19`, PvEBP, PvRBP2b, PvMSP8, PvMSP5, PkP12, PkRIPR, `PkPTRAMP-CSS` = `PkTRAMP-CSS`, Pksera3ag1 = PkSERA3Ag2, PkSSP2)

## import rf model from 4-very-recent-pk-modify
rf_ranger = readRDS(here::here("models/rf_ranger_combo2_very-recent.Rds"))
rf_ranger_plasmodium = readRDS(here::here("models/rf_ranger_combo2_plasmodium.Rds"))
rf_model_macaque = readRDS(here::here("models/rf_ranger_combo_macaque_plasmodium.Rds"))


## script for classification with model and data input

unknown_predict <- function(model_input, data_input) {
  predict(model_input[[3]], data_input[,-1], type = "prob") %>%
    bind_cols(data_input) %>%
    mutate(max = if_else(.pred_new>pull(filter(model_input[[2]], threshold=="max"), threshold_value), "new", "old")) %>%
    mutate(sens90 = if_else(.pred_new>pull(filter(model_input[[2]], threshold=="sens90"), threshold_value), "new", "old")) %>%
    mutate(sens85	 = if_else(.pred_new>pull(filter(model_input[[2]], threshold=="sens85"), threshold_value), "new", "old")) %>%
    mutate(spec85 = if_else(.pred_new>pull(filter(model_input[[2]], threshold=="spec85"), threshold_value), "new", "old")) %>%
    mutate(spec90 = if_else(.pred_new>pull(filter(model_input[[2]], threshold=="spec90"), threshold_value), "new", "old")) %>%
    dplyr::select(sample_id, max, .pred_new,  everything()) %>%
    dplyr::mutate(model_name = paste0(unique(model_input[[1]]$model)))
}



all_predictions_macaque_std = unknown_predict(rf_ranger, to_classify) %>% 
  mutate(pred_equal = if_else(.pred_new>0.5, "seropositive", "seronegative"))

all_predictions_macaque_plas = unknown_predict(rf_ranger_plasmodium, to_classify) %>% 
  mutate(pred_equal = if_else(.pred_new>0.5, "seropositive", "seronegative"))

all_predictions_macaque_plas_adapted = unknown_predict(rf_model_macaque, to_classify) %>% 
  mutate(pred_equal = if_else(.pred_new>0.5, "seropositive", "seronegative"))


adapt_data <- function(data) {
  data %>% 
    left_join(macaque %>% 
                select(sample_id = SampleID, common_name, sex, age, summary_PCR, district,captured_site_latitude, captured_site_longitude,collection_date) %>% distinct()) %>% drop_na(common_name)
}


all_predictions_macaque_std <- all_predictions_macaque_std %>% adapt_data
all_predictions_macaque_plas <- all_predictions_macaque_plas %>% adapt_data
all_predictions_macaque_plas_adapted <- all_predictions_macaque_plas_adapted %>% adapt_data


# standard pk classification ----------------------------------------------

all_predictions_macaque_std %>% 
  mutate(age = factor(age, levels  = c("Juvenile", "Adult"))) %>% 
  mutate(pred_equal = fct_relevel(pred_equal, "seropositive")) %>% 
  
  select(pred_equal) %>% 
  
  gtsummary::tbl_summary(digits = list(gtsummary::all_categorical() ~ c(0, 1))) %>% 
  gtsummary::add_ci() %>% 
  gtsummary::as_flex_table()



# plasmodium standard -----------------------------------------------------

all_predictions_macaque_plas %>% 
  mutate(age = factor(age, levels  = c("Juvenile", "Adult"))) %>% 
  mutate(pred_equal = fct_relevel(pred_equal, "seropositive")) %>% 
  
  select(pred_equal) %>% 
  
  gtsummary::tbl_summary(digits = list(gtsummary::all_categorical() ~ c(0, 1))) %>% 
  gtsummary::add_ci() %>% 
  gtsummary::as_flex_table()



# model adapted for macaques ----------------------------------------------

all_predictions_macaque_plas_adapted %>% 
  mutate(age = factor(age, levels  = c("Juvenile", "Adult"))) %>% 
  mutate(pred_equal = fct_relevel(pred_equal, "seropositive")) %>% 
  
  select(pred_equal) %>% 
  gtsummary::tbl_summary(digits = list(gtsummary::all_categorical() ~ c(0, 1))) %>% 
  
  gtsummary::add_ci() %>% 
  gtsummary::as_flex_table()


# supplementary figure 5: showing the pred new value distributions --------

boxplots_pred_macaque_std <- all_predictions_macaque_std %>% 
  mutate(summary_PCR = factor(summary_PCR, levels = c("Pk", "Plasmodium", "negative"))) %>%
  mutate(age = factor(age, levels  = c("Juvenile", "Adult"))) %>% 
  ggplot(aes(x = age, y = .pred_new))+
  geom_boxplot(width =2, outlier.shape = NA)+
  facet_wrap(~summary_PCR, labeller = label_wrap_gen(width = 20))+
  ggpubr::stat_compare_means()+
  geom_jitter(alpha = 0.6, width = 0.2, 
              mapping = aes(shape = common_name) )+
  theme_minimal()+
  scale_y_continuous(limits = c(0,1.1))+
  labs(y = "prediction values", colour = "", x = "",
       subtitle = "Pk recent: human optimised")+
  theme(panel.border = element_rect(colour = "black",fill = NA),
        axis.text.x = element_text(angle = 45, hjust= 1),
        legend.position = "bottom",
        strip.text = element_text(size = 10))



boxplots_pred_macaque_adapted = all_predictions_macaque_plas_adapted %>% 
  mutate(summary_PCR = factor(summary_PCR, levels = c("Pk", "Plasmodium", "negative"))) %>%
  mutate(age = factor(age, levels  = c("Juvenile", "Adult"))) %>% 
  ggplot(aes(x = age, y = .pred_new))+
  geom_boxplot(width =2, outlier.shape = NA)+
  facet_wrap(~summary_PCR, labeller = label_wrap_gen(width = 20))+
  ggpubr::stat_compare_means()+
  geom_jitter(alpha = 0.6, width = 0.2 , 
              mapping = aes(shape = common_name))+
  theme_minimal()+
  scale_y_continuous(limits = c(0,1.1))+
  labs(y = "prediction values", colour = "", x = "",
       subtitle = "Plasmodium recent: macaque optimised")+
  theme(panel.border = element_rect(colour = "black",fill = NA),
        axis.text.x = element_text(angle = 45, hjust= 1),
        legend.position = "bottom",
        strip.text = element_text(size = 10))

svg(here::here("figures/figure5-supp.svg"),width = 12, height = 6)
print(boxplots_pred_macaque_std +boxplots_pred_macaque_adapted+
        plot_annotation(tag_levels   = "A")+
        plot_layout(guides = "collect")& 
        
        theme(axis.text = element_text(size = 12), 
              axis.title = element_text(size =12),
              plot.tag = element_text(size = 15),
              legend.position = "bottom")) & labs(shape = "")

dev.off()

# main figure 5 the macaque abs -------------------------------------------


source(here::here("scripts-analysis/0_setup.R"))
source(here::here("scripts-analysis/2_data-tidy.R"))

pk_median <- median_wk1_pk %>% 
  dplyr::rename(Antigens = protein_label_clean) %>% 
  filter(Antigens%in%c("Pk8", "PkMSP10", "PkP12", "PkRIPR", "Pksera3ag1", "Pksera3ag2", "PkSSP2", "PkPTRAMP-CSS", "PvFAMA-95", "PvMSP10", "PvEBP", "PvMSP1-19", "PvMSP5", "PvMSP8", "PvRBP2b")) %>% 
  mutate(sample_cat = ("P. knowlesi"))


pv_median <- median_wk1_pv %>% 
  dplyr::rename(Antigens = protein_label_clean) %>% 
  filter(Antigens%in%c("Pk8", "PkMSP10", "PkP12", "PkRIPR", "Pksera3ag1", "Pksera3ag2", "PkSSP2", "PkPTRAMP-CSS", "PvFAMA-95", "PvMSP10", "PvEBP", "PvMSP1-19", "PvMSP5", "PvMSP8", "PvRBP2b")) %>% 
  mutate(sample_cat = ("P. vivax"))


neg_median = cut_off %>% 
  left_join(protein_names %>% dplyr::select(protein_label, protein_label_clean ) %>% distinct()) %>% 
  dplyr::rename(Antigens = protein_label_clean) %>% 
  filter(Antigens%in%c("Pk8", "PkMSP10", "PkP12", "PkRIPR", "Pksera3ag1", "Pksera3ag2", "PkSSP2", "PkPTRAMP-CSS", "PvFAMA-95", "PvMSP10", "PvEBP", "PvMSP1-19", "PvMSP5", "PvMSP8", "PvRBP2b")) %>% 
  mutate(sample_cat = "Thai/Aus malaria naive")

symnum.args <- list(
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
  symbols   = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "ns"))

macaque_rau <- all_mfi_rau %>% 
  filter(RAU_Method=="Adj_loglog"| (Species=="Pk" & RAU_Method=="loglog")|(Species=="Pf" & RAU_Method=="loglog")|(Antigens=="PvMSP10" & RAU_Method=="loglog")) %>% 
  left_join(macaque) %>%drop_na(common_name) %>% 
  mutate(summary_PCR = factor(summary_PCR, levels = c("Pk", "Plasmodium", "negative"))) %>%
  filter(Antigens!="PPF-Pk10") %>% 
  filter(Antigens!="PvPTEX150") %>% 
  filter(Antigens!="PvCSS") %>% 
  filter(Species=="Pk"|Species=="Pv") %>% 
  mutate(Antigens = Antigens %>% 
           replace_values(
             "PkTRAMP-CSS"~ "PkPTRAMP-CSS",
             "PkSERA3Ag2" ~ "Pksera3ag2",
             "PkSERA3ag1" ~ "Pksera3ag1",
             "Pv-fam-a" ~ "PvFAMA-95"
           )) %>% 
  mutate(Antigens = factor(Antigens, levels = c("PkMSP10", "PkP12", "PkRIPR", "PkPTRAMP-CSS", "Pksera3ag1", "Pksera3ag2", "PkSSP2","Pk8", "PvMSP1-19", "PvMSP8", "PvFAMA-95", "PvMSP5", "PvEBP", "PvRBP2b", "PvMSP10"))) %>% 
  ggplot(aes(x = summary_PCR, y = RAU))+
  geom_hline(yintercept = 0.02, linetype = 2)+
  geom_hline(yintercept = 1.95e-05, linetype = 2)+
  geom_hline(data = pk_median, aes(yintercept = median_week1, colour = sample_cat ), linetype = 6, size = 1, alpha = 0.7)+
  geom_hline(data = pv_median, aes(yintercept = median_week1, colour = sample_cat), linetype = 5, size = 1, alpha = 0.85)+
  geom_hline(data = neg_median, aes(yintercept = median_neg, colour = sample_cat), linetype = 6, size = 1, alpha = 0.7)+
  geom_boxplot(outlier.shape = NA, colour = "darkgrey", alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 0.8, colour = "darkgrey") +
  scale_y_log10(labels = scales::comma, limits = c(1.95e-05, 0.035))+
  facet_wrap(~as.factor(Antigens), nrow = 3,
             labeller = as_labeller(\(x) dplyr::case_match(x,
                                                           "PkP12" ~ "PkP12*",
                                                           "PkRIPR" ~"PkRIPR*",
                                                           "PkPTRAMP-CSS" ~ "PkPTRAMP-CSS*",
                                                           "Pksera3ag1" ~ "Pksera3ag1*",
                                                           
                                                           .default=x)))+
  theme_minimal()+
  theme(panel.border = element_rect(colour = "black",fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 12),
        legend.position = "bottom",
        axis.title = element_text(size = 12.5),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 12.5, hjust = 0))+
  labs(x = "PCR group", colour = "reference RAU median")+
  ggpubr:: stat_compare_means(label = "p.signif", vjust = 0.58, symnum.args = symnum.args)+
  scale_color_manual(
    values = c("P. knowlesi" = "#DB4F48", "P. vivax" = "#C49A1C", "Thai/Aus malaria naive" = "#5E9ABF"),
    labels = c(expression(italic("P. knowlesi")), expression(italic("P. vivax")), "Thai/Aus malaria naive")
  )



svg(here::here("figures/figure5.svg"),width = 10, height = 8.5)
print(macaque_rau)
dev.off()

# comparison of any infection vs neg for each antigen ---------------------
all_mfi_rau %>% 
  filter(RAU_Method=="Adj_loglog"| (Species=="Pk" & RAU_Method=="loglog")|(Species=="Pf" & RAU_Method=="loglog")|(Antigens=="PvMSP10" & RAU_Method=="loglog")) %>% 
  left_join(macaque) %>%drop_na(common_name) %>% 
  mutate(summary_PCR = factor(summary_PCR, levels = c("Pk", "Plasmodium", "negative"))) %>%
  filter(Antigens!="PPF-Pk10") %>% 
  filter(Antigens!="PvPTEX150") %>% 
  filter(Antigens!="PvCSS") %>% 
  filter(Species=="Pk"|Species=="Pv") %>% 
  mutate(Antigens = Antigens %>% 
           replace_values(
             "PkTRAMP-CSS"~ "PkPTRAMP-CSS",
             "PkSERA3Ag2" ~ "Pksera3ag2",
             "PkSERA3ag1" ~ "Pksera3ag1",
             "Pv-fam-a" ~ "PvFAMA-95"
           )) %>% 

  select( SampleID,RAU, summary_PCR, Antigens) %>% 
  mutate(RAU = log(RAU)) %>% 
  group_by(Antigens) %>% 
  nest() %>% 
  mutate(
    test = map(data, ~ t.test(pull(filter(.x, summary_PCR=="negative"), "RAU"), pull(filter(.x, summary_PCR!="negative"), "RAU"), paired = FALSE, exact = FALSE)), 
    tidied = map(test, broom::tidy)
  ) %>% 
  unnest(cols = tidied) %>% 
  select(-data, -test) %>% 
  filter(p.value<0.05) %>% arrange(p.value) %>% 
  mutate(adj.p = p.adjust(p.value, method = "bonferroni")) %>% 
  select(adj.p, everything())



# export rau for supplementary file ---------------------------------------


supp_rau <- all_mfi_rau %>% 
  filter(RAU_Method=="Adj_loglog"| (Species=="Pk" & RAU_Method=="loglog")|(Species=="Pf" & RAU_Method=="loglog")|(Antigens=="PvMSP10" & RAU_Method=="loglog")) %>% 
  dplyr::rename(mfi = MFI, rau = RAU) %>% 
  select(-c(Species, RAU_Method, mfi)) %>% 
  pivot_wider(names_from = Antigens, values_from = rau) %>% 
  select(SampleID, PkMSP10, `PvFAMA-95` = `Pv-fam-a`, `PvMSP1-19`, PvEBP, PvRBP2b, PvMSP8, PvMSP5, PkP12, PkRIPR, `PkPTRAMP-CSS` = `PkTRAMP-CSS`, Pksera3ag1 = PkSERA3Ag2, PkSSP2, PvMSP10)

pred_values <- all_predictions_macaque_plas_adapted %>% 
  select(sample_id, model_macaque_pred_new = .pred_new) %>% 
  left_join(all_predictions_macaque_std %>% select(sample_id, pkrecent_pred_new = .pred_new)) %>% 
  left_join(all_predictions_macaque_plas %>% select(sample_id, eitherplasmodium_pred_new = .pred_new)) %>% 
  dplyr::rename(SampleID = sample_id)

#write copy with all information
macaque %>% 
  select(SampleID, summary_PCR, sex, age, collection_date, common_name) %>% 
  left_join(supp_rau) %>% left_join(pred_values) %>% 
  arrange(collection_date) %>% 
  mutate(macaque_id = paste0("macaque-", row_number())) %>% 
  select(SampleID, macaque_id, everything()) %>% 
  write_csv(here::here("processed-data/macaque_supplementary_rau.csv"))


#write copy where there is a different manuscript id 'macaque1'
macaque %>% 
  select(SampleID, summary_PCR, sex, age, collection_date, common_name) %>% 
  left_join(supp_rau) %>% left_join(pred_values) %>%
  arrange(collection_date) %>% 
  mutate(macaque_id = paste0("macaque-", row_number())) %>% 
  select(-c(SampleID, collection_date)) %>% 
  select(macaque_id, everything())%>%   
  write_csv(here::here("processed-data/supplementary_files/macaque_rau.csv"))



         
         

  

