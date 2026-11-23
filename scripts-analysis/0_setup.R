## setup
library(tidyverse)
## have the key to change between the different protein name options
proteins <- read_csv(here::here("protein_info/protein_naming.csv"))

proteins_paper_clean <- proteins$protein_label_clean
proteins_paper <- proteins$protein_label

## ggplot correlations
order_combined <- proteins %>% 
  dplyr::select(order_1, order_2, protein_label_clean) %>% distinct() %>% 
  mutate(prot2 = protein_label_clean) %>% 
  dplyr::rename(prot1 = protein_label_clean)

order_1 <- order_combined %>% dplyr::select(order_1, prot1)
order_2 <- order_combined %>% dplyr::select(order_2, prot2)


proteins_paper_clean <- unique(proteins_paper_clean)


### FUNCTIONS

classification_fn <- function(data, status_choice) {
  data %>% 
    dplyr::left_join(proteins %>% 
                       dplyr::select(protein_label, protein_label_clean, order_1, order_2) %>% 
                       distinct())  %>% 
    filter(protein_label%in%proteins_paper) %>% 
    dplyr::select(-c(protein_label, order_1, order_2)) %>% 
    dplyr::rename(status={{status_choice}}) %>%
    dplyr::select(sample_id, status, protein_label_clean, rau) %>% 
    pivot_wider(names_from = "protein_label_clean", values_from = "rau") %>% 
    mutate(status = as.factor(status))
}

classification_fn_mfi <- function(data, status_choice) {
  data %>% 
    dplyr::left_join(proteins %>% 
                       dplyr::select(protein_label, protein_label_clean, order_1, order_2) %>% 
                       distinct())  %>% 
    filter(protein_label%in%proteins_paper) %>% 
    dplyr::select(-c(protein_label, order_1, order_2)) %>% 
    dplyr::rename(status={{status_choice}}) %>%
    dplyr::select(sample_id, status, protein_label_clean, mfi) %>% 
    pivot_wider(names_from = "protein_label_clean", values_from = "mfi") %>% 
    mutate(status = as.factor(status))
}


dash_rename <- function(data) {
  
  data %>% 
    mutate(combo = str_replace(combo, "MSP1-19", "MSP1_19")) %>% 
    mutate(combo = str_replace(combo, "FAMA-95", "FAMA_95"))
  
}

dash_revert<- function(data) {
  data %>% 
    pivot_longer(-mean) %>% 
    mutate(value = str_replace(value, "MSP1\\_19", "MSP1-19")) %>% 
    mutate(value = str_replace(value, "FAMA_95", "FAMA-95")) %>% 
    pivot_wider(names_from = name, values_from = value)
  
}


