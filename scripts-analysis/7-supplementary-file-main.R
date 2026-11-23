source(here::here("scripts-analysis/0_setup.R"))
source(here::here("scripts-analysis/2_data-tidy.R"))

library(tidyverse)


# main data for model and describing ab -----------------------------------


cats_for_plot %>% 
  select(sample_cat, id) %>% 
  distinct() %>% 
  group_by(sample_cat) %>% 
  mutate(id2 = paste0(sample_cat,"-", row_number())) %>% write_csv(here("sample-info/key_destudy_id.csv"))


key <- read_csv(here("sample-info/key_destudy_id.csv")) %>% 
  select(id, id2)

cats_for_plot %>% 
  select(protein_label_clean, sample_cat, id, day_cat, rau, day) %>% 
  filter(protein_label_clean%in%proteins_paper_clean) %>% 
  pivot_wider(names_from = protein_label_clean, values_from = rau) %>% 
  left_join(key) %>% 
  select(id2, day_cat, day, everything()) %>% select(-c("id", "sample_cat")) %>% 
  dplyr::rename(id_ref = id2) %>% 
  write_csv(here::here("processed-data/supplementary_files/", "main_rau.csv"))

