library(tidyverse)
library(tidymodels)
#--------------- JUST TIDYING UP SOME OF PROTEIN NAMES ETC
proteins <- read_csv(here::here("protein_info/protein_naming.csv"))

proteins_paper_clean <- proteins$protein_label_clean
proteins_paper <- proteins$protein_label
proteins_paper_clean <- unique(proteins_paper_clean)


#############################################################################

# functions

##############################################################################

# ----- FUNCTION TO PLOT ROC CURVES

roc_plot <- function(data) {
  data %>% 
    ggplot(aes(x = 1-specificity, y = sensitivity, group = model, colour = model))+
    geom_path(size = 1.5)+
    geom_abline(linetype = 2)+
    coord_equal()+
    scale_x_continuous(limits = c(0,1), breaks = c(0.2, 0.4, 0.6, 0.8,1), expand = c(0,0))+
    scale_y_continuous(limits = c(0,1), breaks = c(0,0.2, 0.4, 0.6, 0.8,1), expand = c(0,0))+
    theme_linedraw()+
    theme(panel.grid = element_line(colour = "grey"))
}

#############################################################################

# ----- FUNCTION THAT PLOTS THE CONVERGENCE OF SENS/SPEC AND .THRESHOLD

thresholds_plot <- function(data) {
  data %>% 
    dplyr::select(.threshold, specificity, sensitivity, model) %>% 
    pivot_longer(c(specificity, sensitivity)) %>% 
    ggplot(aes( x = .threshold, y = value, colour = name))+
    geom_path(size = 1.5)+
    scale_colour_manual(values = c("#0A7373", "#C43302"))+
    coord_equal()+
    scale_x_continuous(limits = c(0,1), breaks = c(0.2, 0.4, 0.6, 0.8,1), expand = c(0,0))+
    scale_y_continuous(limits = c(0,1), breaks = c(0,0.2, 0.4, 0.6, 0.8,1), expand = c(0,0))+
    theme_linedraw()+
    theme(panel.grid = element_line(colour = "grey"))+
    labs(y = "sensitivity/specificity")
}  

###############################################################################

# ----- FUNCTION TO GET TABLE WITH .threshold VALUES FOR THE FIVE CUTOFFS


get_thresholds <- function(data) {
  data %>% 
    dplyr::select(.threshold, specificity, sensitivity) %>% 
    dplyr::mutate(across(where(is.numeric), round, 2)) %>% 
    distinct() %>% 
    dplyr::mutate(max = if_else(specificity==sensitivity, .threshold, NA)) %>% 
    dplyr::mutate(spec85 = if_else(specificity==0.850, .threshold, NA)) %>% 
    dplyr::mutate(spec90 = if_else(specificity==0.90, .threshold, NA)) %>% 
    dplyr::mutate(sens85 = if_else(sensitivity==0.85, .threshold, NA)) %>%  
    dplyr::mutate(sens90 = if_else(sensitivity==0.90, .threshold, NA)) %>%  
    dplyr::mutate(sens90 = if_else(sensitivity==0.90, .threshold, NA)) %>% 
    pivot_longer(c(max, spec85, spec90, sens85, sens90)) %>% 
    drop_na(value) %>% 
    dplyr::select(name, value) %>% 
    dplyr::group_by(name) %>% 
    summarise(value = mean(value)) %>% 
    dplyr::mutate(value = round(value, 3)) %>% 
    dplyr::rename(threshold = name, threshold_value = value)
  
}

###############################################################################

# ---- FUNCTION WORKING ON LIST OBJECT FROM THE 'model-objects' FOLDER
#--------model-objects are in order (metrics table, threshold values, model fit)

test_predict <- function(model_input) {
  predict(model_input[[3]], test[,-1], type = "prob") %>%  #remove col 1 because this is known status
    mutate(sample_id = test_id) %>% #put sampleIDs back
    bind_cols(test) %>% #merge RAU data
    mutate(max = if_else(.pred_new>pull(filter(model_input[[2]], threshold=="max"), threshold_value), "new", "old")) %>% 
    mutate(sens90 = if_else(.pred_new>pull(filter(model_input[[2]], threshold=="sens90"), threshold_value), "new", "old")) %>% 
    mutate(sens85	 = if_else(.pred_new>pull(filter(model_input[[2]], threshold=="sens85"), threshold_value), "new", "old")) %>% 
    mutate(spec85 = if_else(.pred_new>pull(filter(model_input[[2]], threshold=="spec85"), threshold_value), "new", "old")) %>% 
    mutate(spec90 = if_else(.pred_new>pull(filter(model_input[[2]], threshold=="spec90"), threshold_value), "new", "old")) %>% 
    mutate(assess = case_when(
      max=="new" & status=="new" ~ "true positive",
      max=="new" & status=="old" ~ "false positive",
      max=="old" & status=="new" ~ "false negative",
      max=="old" & status=="old" ~ "true negative",
    )) %>% 
    dplyr::select(sample_id, assess, status, max, .pred_new, protein_subset, everything())
}

################################################################################

# ------ FUNCTION SPECIFIC FOR TEST DATA TO PLOT TP/TN/FP/FN

test_ab_dist <- function(data) {
  data %>% 
    pivot_longer(protein_subset) %>% 
    ggplot(aes(x = assess, y = value, fill = assess))+
    geom_boxplot()+
    facet_wrap(~name)+ theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_log10()+
    scale_fill_manual(values = c("#A60530", "#D9851E", "#84B1D9", "#3B668C" ))+
    labs(y = "RAU")
}

# --- FUNCTION SPECIFIC FOR TEST DATA TO PLOT CONT MATRIX

cont_matrix <- function(data) {
  data %>% dplyr::select(status, max) %>% 
    count(status, max) %>% 
    ggplot(aes(x = max, y = fct_rev(status)))+
    geom_tile(colour = "white", size = 2, fill = "lightgrey")+
    geom_text(aes(label = n), size = 5)+
    theme_minimal()+
    labs(y = "predict", x = "truth")+
    scale_x_discrete(position = "top")+
    coord_equal()
}


################################################################################

# -----FUNCTION TO PREDICT ON UNKNOWN DATA (NOT THE TEST DATA)
# ----------data_input assumes sample_id is in first col

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




