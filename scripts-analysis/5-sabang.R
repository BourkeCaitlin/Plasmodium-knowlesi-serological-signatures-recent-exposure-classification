library(tidyverse)
library(tidymodels)
library(patchwork)

## maps of aceh

#map of sabang/pulau weh + n sumatra
aceh_boundaries_adm1 <- rgeoboundaries::geoboundaries(country = "indonesia",
                                                      adm_lvl =  1)
aceh_boundaries_adm2 <- rgeoboundaries::geoboundaries(country = "indonesia",
                                                      adm_lvl =  2)
# create map of island--------- 
# here only using rgeoboundaries admin two levels
# use rgeospatial to add scale bar and N direction
sabang_map <- aceh_boundaries_adm2 %>% 
  ggplot()+
  geom_sf(data = aceh_boundaries_adm2, colour = "darkgrey")+
  # annotate(x = 95.32446179607246,
  #          y = 5.892733850193278,
  #          geom = "point",
  #          size = 1)+
  
  # annotate(geom = "text",
  #          x = 95.32446179607246,
  #          y = 5.9,
  #          label = "Sabang",
  #          size = 4.6)+
  
  # annotate(geom = "text",
  #          y= 5.872982814192166, x = 95.25,
  #          label = "Iboih",
  #          fontface = 'italic',
  #          alpha = 0.6,
  #          size = 3.6)+
  # 
  # annotate(geom = "text",
  #          x = 95.297,
  #          y = 5.797,
  #          label = "Sukajaya",
  #          fontface = 'italic',
  #          alpha = 0.6,
  #          size = 3.6)+
  # 
  # annotate(geom = "text",
  #          x = 95.315,
  #          y = 5.835,
  #          label = "Sukakarya",
  #          alpha = 0.6,
  #          fontface = 'italic',
  #          size = 3.6)+
  annotate(geom = "text",
           x = 95.315,
           y = 5.825,
           label = "SABANG",
           alpha = 0.9,
           fontface = 'italic',
           size = 8)+
  
  ggspatial::annotation_scale(location = "tr", text_cex = 1.2 ) +
  ggspatial::annotation_north_arrow(location = "br", which_north = "true")+
  scale_y_continuous(limits = c(5.6, 5.92))+
  scale_x_continuous(limits = c(95.1, 95.52))+
  theme(panel.background =  element_rect( "#85BAEA"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())



# create map of northern sumatra to inset to show where sabang actually is  -----
sumatra <-  ggplot()+
  geom_sf(data = aceh_boundaries_adm1, colour = "darkgrey", fill = "lightgrey", alpha = 0.8)+
  ggspatial::annotation_scale(location = "bl", text_cex = 0.6) +
  scale_y_continuous(limits = c(2.5, 7.5))+
  scale_x_continuous(limits = c(93, 100))+
  annotate("point", x = 95.3, y =5.85 , size = 4, pch = 1, color = "darkblue", stroke = 1)+
  annotate(geom = "text",
           x = 96.5,
           y = 7.2,
           label = "SUMATRA",
           fontface = 'italic',
           size = 6)+
  
  annotate(geom = "text",
           x = 97,
           y = 4.5,
           label = "Aceh",
           fontface = 'italic',
           alpha = 0.6,
           size = 2.5)+  
  annotate(geom = "text",
           x = 99,
           y = 3,
           alpha = 0.6,
           label = "North\nSumatra",
           fontface = 'italic',
           size = 2.5)+  
  theme_void()+
  theme(
    plot.background = element_rect(fill = "#6493bd", color = "black", linetype = "dotted", linewidth = 1),
    panel.background = element_rect(fill = "transparent"))


# svg(here::here("figures/figure4-sabangmap.svg"),width = 8, height = 8)
# print(sabang_map + patchwork::inset_element(sumatra, left = 0.01, bottom = 0.01, right = 0.5, top = 0.4))
# dev.off()

## import data from 5-process-zoomagic -------
all_zoomagic = read_csv(here::here("processed-data/zoomagic_ehi.csv"))

## read in meta data shared by matt -------
zoomagic_meta = read_rds(here::here("sample-info/ZOOMAGIC_combined_141025.rds"))

zoomagic_meta %>% 
  select(studysite, code, contains("occup")) %>% filter(str_detect(studysite, "Sabang")) %>% 
  filter(occup3=="")

## extract the RAU appropriate for the right species

all = all_zoomagic %>% 
  filter(RAU_Method=="Adj_loglog"| (Species=="Pk" & RAU_Method=="loglog")) %>% 
  dplyr::rename(mfi = MFI, rau = RAU) %>% 
  select(-c(Species, RAU_Method, mfi)) %>% 
  pivot_wider(names_from = Antigens, values_from = rau)


## import rf model from 4-very-recent-pk-modify & either plasmodium
rf_ranger = readRDS(here::here("models/rf_ranger_combo2_very-recent.Rds"))
rf_ranger_plasmodium = readRDS(here::here("models/rf_ranger_combo2_plasmodium.Rds"))



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


all = all %>%
  dplyr::rename(`PvFAMA-95` = matches("Pv-fam-a")) 

to_predict = all %>% dplyr::rename(sample_id = SampleID)

#create merged db with new Pk and old bead re-run Pv
split = to_predict %>% 
  group_split(beads)

split[[4]] <- split[[1]] %>% dplyr::select(sample_id, Plate, PkMSP10, PkP12) %>% 
  left_join(split[[2]] %>% 
              dplyr::select(-c(PkMSP10, PkP12))) %>% mutate(beads = "combined2026")

to_predict <- bind_rows(split)

all_predictions = unknown_predict(rf_ranger, to_predict) %>% 
  mutate(pred_equal = if_else(.pred_new>0.5, "seropositive", "seronegative"))

all_predictions_plasmodium = unknown_predict(rf_ranger_plasmodium, to_predict) %>% 
  mutate(pred_equal = if_else(.pred_new>0.5, "seropositive", "seronegative"))

## for now use pred_equal which means 50% of trees voted for recent infection

all_predictions %>% 
  filter(str_detect(sample_id, "S")) %>% 
  group_by(beads) %>% janitor::tabyl(beads,pred_equal) %>% 
  janitor::adorn_percentages() %>% 
  janitor::adorn_pct_formatting(digits= 1) %>% 
  janitor::adorn_ns()

all_predictions_plasmodium %>% 
  filter(str_detect(sample_id, "S")) %>% 
  group_by(beads) %>% janitor::tabyl(beads,pred_equal) %>% 
  janitor::adorn_percentages() %>% 
  janitor::adorn_pct_formatting(digits= 1) %>% 
  janitor::adorn_ns()

## correlation plot of pred pk vs pred plasmodium----------------
pred_pk_pred_plasmodium = 
  all_predictions %>% 
  filter(str_detect(sample_id, "S")) %>% 
  mutate(sample_id= case_when(
    str_detect(sample_id, "SBK-130") ~ "SBK-103",
    T ~ sample_id
  )) %>% 
  filter(!sample_id%in%c("SBS-006a", "SBS-007a", "SBS-038a", "SBK-042B")) %>% 
  select(sample_id, pkrecent = .pred_new, beads ) %>% filter(beads=="newbatch") %>% 
  left_join(all_predictions_plasmodium %>% 
              mutate(sample_id= case_when(
                str_detect(sample_id, "SBK-130") ~ "SBK-103",
                T ~ sample_id
              )) %>% 
              # filter(!sample_id%in%c("SBS-006a", "SBS-007a", "SBS-038a", "SBK-042B")) %>% 
              select(sample_id, plasmodium = .pred_new, beads) %>% 
              filter(beads=="newbatch") %>% 
              distinct()) %>% 
  mutate(pcr_current = if_else(sample_id%in%c("SBK-006", "SBI-007", "SBI-048", "SBS-099"),"PCR Pk","neg PCR")) %>%
  arrange(desc(pcr_current)) %>% 
  ggplot(aes(x = plasmodium, y = pkrecent))+
  geom_point(aes(colour = pcr_current,shape = pcr_current),size = 2.5, alpha = 0.4)+
  ggpubr::stat_cor(p.accuracy = 0.01, r.accuracy = 0.01,method = "spearman", 
                   colour = "black",label.y = 1, inherit.aes = TRUE, size = 4.5,cor.coef.name = "rho")+
  theme_minimal()+
  labs(y = "prediction value 'Pk recent'", x = "prediction value 'either Plasmodium'")+
  coord_equal()+
  scale_colour_manual(values = c("black", "red"))+
  labs(colour = "", shape = "")+
  theme(panel.background = element_rect(colour = "black"),
        legend.position = "bottom")

# svg(here::here("figures/figure4-corrpreds.svg"),width = 3, height = 3)
# print(pred_pk_pred_plasmodium)
# dev.off()

############# 
# create dataset for analysis-------------
class_prot = c('PkMSP10', 'PvMSP5','PvFAMA-95', 'PvMSP1-19', 'PvRBP2b', 'PvMSP8', 'PvEBP')
epi_var = c("studysite", "age", "sex", "village", "agegroup", "plat", "monk", "pcr_eimb", "thrombocyt", "forestnh", "clearfornh", "insrepel", "eavesgaps", "occup3", "bednet", "sleepout", "forestin")

zoomagic_analysis = all_predictions %>% 
  filter(str_detect(sample_id, "S")) %>% 
  mutate(sample_id= case_when(
    str_detect(sample_id, "SBK-130") ~ "SBK-103",
    T ~ sample_id
  )) %>% 
  filter(!sample_id%in%c("SBS-006a", "SBS-007a", "SBS-038a", "SBK-042B")) %>% 
  left_join(zoomagic_meta %>% dplyr::rename(sample_id = code)) %>% distinct() %>% filter(beads=="newbatch") %>% 
  select(sample_id, pred_equal, sens85, sens90, .pred_new, all_of(class_prot), all_of(epi_var)) %>% 
  mutate(across(epi_var[!epi_var%in%c("age", "plat")], as.factor)) %>% 
  mutate(across(c("pred_equal", "sens85", "sens90"), as.factor)) %>% 
  mutate(pred_equal = fct_relevel(pred_equal, "seropositive")) %>% 
  mutate(sens85 = fct_relevel(sens85, "new")) %>% 
  mutate(sens90 = fct_relevel(sens90, "new")) 

zoomagic_analysis_plas = all_predictions_plasmodium %>% 
  filter(str_detect(sample_id, "S")) %>% 
  mutate(sample_id= case_when(
    str_detect(sample_id, "SBK-130") ~ "SBK-103",
    T ~ sample_id
  )) %>% 
  filter(!sample_id%in%c("SBS-006a", "SBS-007a", "SBS-038a", "SBK-042B")) %>% 
  left_join(zoomagic_meta %>% dplyr::rename(sample_id = code)) %>% distinct() %>% filter(beads=="newbatch") %>% 
  select(sample_id, pred_equal, sens85, sens90, .pred_new, all_of(class_prot), all_of(epi_var)) %>% 
  mutate(across(epi_var[!epi_var%in%c("age", "plat")], as.factor)) %>% 
  mutate(across(c("pred_equal", "sens85", "sens90"), as.factor)) %>% 
  mutate(pred_equal = fct_relevel(pred_equal, "seropositive")) %>% 
  mutate(sens85 = fct_relevel(sens85, "new")) %>% 
  mutate(sens90 = fct_relevel(sens90, "new")) 



zoomagic_analysis_long = zoomagic_analysis %>% 
  select(sample_id, all_of(class_prot),pred_equal, .pred_new) %>% 
  pivot_longer(-c(sample_id, pred_equal, .pred_new)) %>% 
  
  left_join(zoomagic_meta %>% dplyr::rename(sample_id = code)) %>% 
  mutate(pcr_eimb = if_else(pcr_eimb=="P.genus", "Neg", pcr_eimb)) %>% 
  mutate(plot_class = case_when(
    pcr_eimb=="Pk" & pred_equal=="seropositive" ~ "seropositive and PCR positive",
    pcr_eimb=="Neg" & pred_equal=="seronegative" ~"seronegative and PCR negative",
    pcr_eimb=="Pk" & pred_equal=="seronegative" ~"seronegative and PCR positive",
    pcr_eimb=="Neg" & pred_equal=="seropositive" ~"seropositive and PCR negative",
  ))%>% 
  select(sample_id, pred_equal,  .pred_new, name, value, all_of(epi_var), plot_class)

zoomagic_analysis_long_plas = zoomagic_analysis_plas %>% 
  select(sample_id, all_of(class_prot),pred_equal, .pred_new) %>% 
  pivot_longer(-c(sample_id, pred_equal, .pred_new)) %>% 
  left_join(zoomagic_meta %>% dplyr::rename(sample_id = code)) %>% 
  mutate(pcr_eimb = if_else(pcr_eimb=="P.genus", "Neg", pcr_eimb)) %>% 
  mutate(plot_class = case_when(
    pcr_eimb=="Pk" & pred_equal=="seropositive" ~ "seropositive and PCR positive",
    pcr_eimb=="Neg" & pred_equal=="seronegative" ~"seronegative and PCR negative",
    pcr_eimb=="Pk" & pred_equal=="seronegative" ~"seronegative and PCR positive",
    pcr_eimb=="Neg" & pred_equal=="seropositive" ~"seropositive and PCR negative",
  ))%>% 
  select(sample_id, pred_equal,  .pred_new, name, value,plot_class, all_of(epi_var))

### here try to figure out how to combine both models...maybe just do RAU plasmodium recent in supplementary??
combined_long_class <- zoomagic_analysis_long %>% 
  left_join(zoomagic_analysis_long_plas %>% select(sample_id, plasmodium = plot_class) %>% distinct()) %>% 
  filter(plot_class!=plasmodium) %>% 
  select(sample_id,pred_equal, plot_class, plasmodium) 


## calculate study stats-------

zoomagic_analysis %>% 
  summarise(median(age))

zoomagic_analysis %>% 
  summarise(median = median(age),
            lower = quantile(age, 1/4),
            upper = quantile(age, 3/4))


zoomagic_analysis %>% janitor::tabyl(sex)

zoomagic_analysis %>% janitor::tabyl(studysite)


## sex stratifications --------------
zoomagic_analysis %>% 
  select(pred_equal) %>% 
  gtsummary::tbl_summary() %>% 
  gtsummary::add_ci()

zoomagic_analysis_plas %>% 
  select(pred_equal, sex) %>% 
  gtsummary::tbl_summary(by = pred_equal,percent = "row") %>% 
  gtsummary::add_ci() %>% 
  gtsummary::add_p()

zoomagic_analysis %>% 
  select(pred_equal, sex) %>% 
  gtsummary::tbl_summary(by = pred_equal,percent = "row") %>% 
  gtsummary::add_ci() %>% 
  gtsummary::add_p()

zoomagic_analysis %>% 
  select(pred_equal, sex, agegroup) %>% 
  mutate(cat = paste0(sex, "-", agegroup)) %>% 
  gtsummary::tbl_summary(by = pred_equal) %>% 
  gtsummary::add_ci() %>% 
  gtsummary::add_p()

zoomagic_analysis %>% 
  filter(agegroup=="<15") %>% 
  select(pred_equal, sex) %>% 
  gtsummary::tbl_summary(by = sex) %>% 
  gtsummary::add_ci() %>% 
  gtsummary::add_p()

zoomagic_analysis %>% 
  filter(agegroup=="15-49") %>% 
  select(pred_equal, sex) %>% 
  gtsummary::tbl_summary(by = sex) %>% 
  gtsummary::add_ci() %>% 
  gtsummary::add_p()

zoomagic_analysis %>% 
  filter(agegroup=="≥50") %>% 
  select(pred_equal, sex) %>% 
  
  gtsummary::tbl_summary(by = sex) %>% 
  gtsummary::add_ci() %>% 
  gtsummary::add_p()

zoomagic_analysis_plas %>% 
  filter(agegroup=="≥50") %>% 
  select(pred_equal, sex) %>% 
  
  gtsummary::tbl_summary(by = sex) %>% 
  gtsummary::add_ci() %>% 
  gtsummary::add_p()



# figure 4B - comparison of PCR positive and seropositive
plot_rau_zoomagic = zoomagic_analysis_long %>% 
  mutate(sex2 = if_else(sex=="Laki-laki", "Male", "Female")) 


RAU_zoomagic <- plot_rau_zoomagic %>% 
  filter(sex=="Perempuan") %>% 
  ggplot(aes(x = agegroup, y = value, colour = plot_class))+
  geom_hline(yintercept = 1.95e-05, linetype = 2)+
  geom_hline(yintercept = 0.02, linetype = 2)+
  geom_jitter(data = filter(plot_rau_zoomagic, plot_class=="seronegative and PCR negative"), alpha = 0.6, size = 0.9, height = 0, width = 0.2)+
  geom_jitter(data = filter(plot_rau_zoomagic, plot_class=="seronegative and PCR positive"), alpha = 0.7, size = 2, , height = 0, width = 0.2, shape = 17)+
  geom_jitter(data = filter(plot_rau_zoomagic, plot_class=="seropositive and PCR positive"), alpha = 0.8, size = 2, width = 0.2, height = 0, shape = 8)+
  geom_jitter(data = filter(plot_rau_zoomagic, plot_class=="seropositive and PCR negative"), alpha = 1, size = 1, width = 0.2, height = 0)+
  scale_colour_manual(values = c("darkgrey", "#006ECF", "#D92B2B", "#92063C"))+
  scale_size_manual(values = c(1, 2))+
  scale_y_log10(labels = scales::comma)+
  labs(y = "RAU", colour = "", shape = "")+
  facet_grid(sex2~name)+
  theme_minimal()+
  guides(color = guide_legend(nrow = 2))+
  theme(panel.background = element_rect(colour = "black"),
        legend.position = "bottom",
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12))


# svg(here::here("figures/figure4-RAU.svg"),width = 10, height = 5)
# print(RAU_zoomagic)
# dev.off()
# 
### create patchwork of figure-------

(patchwork::wrap_plots(sabang_map + patchwork::inset_element(sumatra, left = 00.001, bottom = 0.00, right = 0.37, top = 0.37)) + pred_pk_pred_plasmodium)/RAU_zoomagic + plot_annotation(tag_levels  = "A") + plot_layout(widths = c(5,2), heights = c(3,2.5))


svg(here::here("figures/figure4-all.svg"),width = 10, height = 10)
print((wrap_plots(sabang_map + labs(tag = "A")+
                    patchwork::inset_element(sumatra, left = 0.01, bottom = 0.01, right = 0.5, top = 0.5)) +
         pred_pk_pred_plasmodium+labs(tag = "B"))/
        RAU_zoomagic + labs(tag = "C") + 
        # plot_annotation(tag_levels  = "A") + 
        plot_layout(widths = c(5,2), heights = c(3,2.5))) & 
  theme(plot.tag = element_text(size = 22))
dev.off()



zoomagic_analysis_reg = zoomagic_analysis %>% 
  mutate(ocu_cat = case_when(
    occup3=="Agriculture" ~ "Agriculture",
    occup3=="" ~ NA,
    occup3==NA~NA,
    T ~ "other"
  )) %>%  
  mutate(occup3 = as.factor(occup3)) %>% 
  mutate(pred_equal = if_else(pred_equal=="seropositive", 1,0)) %>% 
  mutate(sens85 = if_else(sens85=="new", 1,0)) %>% 
  mutate(sens90 = if_else(sens90=="new", 1,0)) %>% 
  
  mutate(sex = fct_relevel(sex, "Perempuan")) %>% 
  mutate(insrepel = fct_relevel(insrepel, "Y")) %>% 
  mutate(monk = if_else(monk=="Y", "Y", "N")) %>% 
  mutate(clearfornh = if_else(clearfornh=="Y", "Y", "N")) %>% 
  mutate(forestnh = if_else(forestnh=="Y", "Y", "N")) %>% 
  
  mutate(monk = fct_relevel(monk, "N")) %>% 
  mutate(ocu_cat = fct_relevel(ocu_cat, "other")) %>% 
  mutate(bednet = if_else(bednet=="Y", "Y", "N")) %>% 
  
  mutate(bednet = fct_relevel(bednet, "N")) %>% 
  mutate(sleepout = if_else(sleepout=="Y", "Y", "N")) %>% 
  
  mutate(sleepout = fct_relevel(sleepout, "N")) %>% 
  mutate(forestin = if_else(forestin=="Y", "Y", "N")) %>% 
  
  mutate(forestin = fct_relevel(forestin, "N"))

### demographics split

zoomagic_analysis_reg %>% 
  gtsummary::tbl_summary(include = c("age", "sex", "monk", "thrombocyt", "ocu_cat", "clearfornh", "forestnh","bednet",  "sleepout", "forestin"),
                         digits = gtsummary::all_categorical() ~ 1,
                         missing_text = "Unknown",
                         missing_stat = "{N_miss}")





#####
zoomagic_analysis_reg

uv = zoomagic_analysis_reg %>% 
  select(pred_equal, age, sex, monk, thrombocyt, ocu_cat, forestnh, clearfornh, forestin, bednet, sleepout) %>% 
  gtsummary::tbl_uvregression(                         
    method = glm,                           
    y = pred_equal,                            
    method.args = list(family = binomial),  pvalue_fun = gtsummary::label_style_pvalue(digits = 2),
    exponentiate = TRUE                     
  ) %>% 
  gtsummary::add_nevent() 



prop_sero <- zoomagic_analysis_reg %>% 
  gtsummary::tbl_summary(include = c("age", "sex", "monk", "thrombocyt", "ocu_cat", "clearfornh", "forestnh","bednet",  "sleepout", "forestin"),
                         by = pred_equal,digits = gtsummary::all_categorical() ~ 1,
                         missing_text = "Unknown",
                         missing_stat = "{N_miss}") %>% 
  gtsummary::as_flex_table()

# flextable::save_as_docx(
#   "prop_sero" = prop_sero, 
#   path = here::here("figures/prop_sero-prop.docx"))



#####



# export data for supplementary data 2 -----------------------------------


sabang_export <- zoomagic_analysis %>% 
  select(sample_id, recentpk_pred_new=.pred_new, PkMSP10, PvMSP5, `PvFAMA-95`, `PvMSP1-19`, PvRBP2b, PvMSP8, PvEBP, sex, agegroup) %>% 
  left_join(zoomagic_analysis_plas %>% 
              select(sample_id, eitherplasmodium_pred_new=.pred_new)) %>% 
  arrange(agegroup, sex) %>% 
  mutate(sabang_id = (paste0('sabang-', row_number()))) %>% 
  select(sabang_id,sample_id,  sex, agegroup,contains('pred'), everything()) 

sabang_export %>% 
  write_csv(here::here("processed-data/sabang_supplementary_all.csv"))

sabang_export %>% 
  select(-sample_id) %>% 
  write_csv(here::here("processed-data/supplementary_files/sabang_rau.csv"))

