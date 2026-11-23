###############################################################################
###############################################################################
library(tidyverse)
library(tidymodels)
library(patchwork)
#--------------- ALL TIDIED DATA
classification = readRDS(here::here("processed-data/classification.rds"))

#--------------- JUST TIDYING UP SOME OF PROTEIN NAMES ETC
proteins <- read_csv(here::here("protein_info/protein_naming.csv"))

proteins_paper_clean <- proteins$protein_label_clean
proteins_paper <- proteins$protein_label
proteins_paper_clean <- unique(proteins_paper_clean)

protein_subset = unique(classification$protein_label_clean)

#############################################################################

#------FUNCTION TO CHOOSE WHAT IS STATUS
#------------this will either be (classify_month_pk, classify_month_plasmodium, or classify_any_pk) from classification object

classification_fn <- function(data, status_choice) {
  data %>% 
    dplyr::left_join(proteins %>% 
                       dplyr::select(protein_label, protein_label_clean, order_1, order_2) %>% 
                       distinct())  %>% 
    # filter(protein_label%in%proteins_paper) %>% 
    dplyr::select(-c(protein_label, order_1, order_2)) %>% 
    dplyr::rename(status={{status_choice}}) %>%
    dplyr::select(sample_id, status, protein_label_clean, rau) %>% 
    pivot_wider(names_from = "protein_label_clean", values_from = "rau") %>% 
    mutate(status = as.factor(status))
}
###############################################################################

### here we are defining the training object
### classification_fn() is defined above - can differ the second part to choose the status (the two different scenarios in the paper)

data = classification_fn(classification, classify_month_plasmodium) %>% 
  dplyr::select(sample_id, status, proteins_paper_clean)

###----------REMOVE ACUTE INFECTIONS
data = data %>%
  filter(!str_detect(sample_id, "-0")) %>%
  filter(!str_detect(sample_id, "w0")) %>% 
  mutate(status = factor(status, levels = c("new", "old"))) 

data_id = data %>% pull(sample_id)

training  = data %>% dplyr::select(-sample_id)




###-------------------- classification models with very-recently --------------
source(here::here("scripts-analysis/0_setup-classification.R"))



##-----------RECIPE
## Recipe is using factor STATUS 


training_macaque = training %>% dplyr::select(status,
                                              PkP12,
                                              PkRIPR,
                                              `PkPTRAMP-CSS`,
                                              Pksera3ag1)

recipe_macaque = recipe(status~., data = training_macaque) %>% 
  step_log(all_numeric())


training_list = list(training_macaque)
recipe_list = list(recipe_macaque)

## - model evaluation with workflowsets

# take forward all 6 models

##########################################

lr = logistic_reg(mode = "classification", engine = "glm")

lr_workflow = NULL
lr_resampled = NULL
for (i in 1:length(recipe_list)) {
  lr_workflow[[i]] = workflow() %>% add_model(lr) %>% add_recipe(recipe_list[[i]])
  
  lr_resampled[[i]] = lr_workflow[[i]] %>% fit_resamples(resamples = vfold_cv(training_list[[i]], v = 10, repeats = 5, strata = status),
                                                         control = control_resamples(save_pred = T,
                                                                                     save_workflow = T),
                                                         metrics = metric_set(roc_auc, pr_auc, accuracy, sensitivity, specificity, brier_class))
  
}

roc_lr  = NULL
metrics_lr = NULL
for (i in 1:length(lr_resampled)) {
  roc_lr[[i]] = lr_resampled[[i]] %>% 
    collect_predictions(summarize = F) %>% 
    roc_curve(status, .pred_new, event_level = "first") %>% 
    mutate(combo = paste0("logisticregression_combo", i))
  
  metrics_lr[[i]] = lr_resampled[[i]] %>% 
    collect_metrics() %>% 
    mutate(combo = paste0("logisticregression_combo", i))
}


##########################################
library(discrim)
lda = discrim_linear() %>% 
  set_engine("MASS",
             seed = 1234) %>% #sets seed
  set_mode("classification")

lda_workflow = NULL
lda_resampled = NULL
for (i in 1:length(recipe_list)) {
  lda_workflow[[i]] = workflow() %>% add_model(lda) %>% add_recipe(recipe_list[[i]])
  
  lda_resampled[[i]] = lda_workflow[[i]] %>% fit_resamples(resamples = vfold_cv(training_list[[i]], v = 10, repeats = 5, strata = status),
                                                           control = control_resamples(save_pred = T,
                                                                                       save_workflow = T),
                                                           metrics = metric_set(roc_auc, pr_auc, accuracy, sensitivity, specificity, brier_class))
  
}

roc_lda  = NULL
metrics_lda = NULL
for (i in 1:length(lda_resampled)) {
  roc_lda[[i]] = lda_resampled[[i]] %>% 
    collect_predictions(summarize = F) %>% 
    roc_curve(status, .pred_new, event_level = "first") %>% 
    mutate(combo = paste0("lda_combo", i))
  
  metrics_lda[[i]] = lda_resampled[[i]] %>% 
    collect_metrics() %>% 
    mutate(combo = paste0("lda_combo", i))
}
########################################
rf = rand_forest(mode = "classification", 
                 # mtry = tune(),
                 # min_n = tune(),
                 trees = 5000) %>% 
  set_engine("ranger", importance = "permutation")



rf_workflow = NULL
rf_resampled = NULL
for (i in 1:length(recipe_list)) {
  rf_workflow[[i]] = workflow() %>% add_model(rf) %>% add_recipe(recipe_list[[i]])
  
  rf_resampled[[i]] = rf_workflow[[i]] %>% fit_resamples(resamples = vfold_cv(training_list[[i]], v = 10, repeats = 5, strata = status),
                                                         control = control_resamples(save_pred = T,
                                                                                     save_workflow = T),
                                                         metrics = metric_set(roc_auc, pr_auc, accuracy, sensitivity, specificity, brier_class))
  
}

roc_rf  = NULL
metrics_rf = NULL
for (i in 1:length(rf_resampled)) {
  roc_rf[[i]] = rf_resampled[[i]] %>% 
    collect_predictions(summarize = F) %>% 
    roc_curve(status, .pred_new, event_level = "first") %>% 
    mutate(combo = paste0("rf_combo", i))
  
  metrics_rf[[i]] = rf_resampled[[i]] %>% 
    collect_metrics() %>% 
    mutate(combo = paste0("rf_combo", i))
}

######################################

bt = boost_tree(
  # tree_depth = tune(),
  # learn_rate = tune(),
  # loss_reduction = tune(),
  # min_n = tune(),
  trees = 5000) %>%
  set_engine("xgboost", event_level = "second") %>% 
  set_mode("classification") %>% 
  translate()

bt_workflow = NULL
bt_resampled = NULL
for (i in 1:length(recipe_list)) {
  bt_workflow[[i]] = workflow() %>% add_model(bt) %>% add_recipe(recipe_list[[i]])
  
  bt_resampled[[i]] = bt_workflow[[i]] %>% fit_resamples(resamples = vfold_cv(training_list[[i]], v = 10, repeats = 5, strata = status),
                                                         control = control_resamples(save_pred = T,
                                                                                     save_workflow = T),
                                                         metrics = metric_set(roc_auc, pr_auc, accuracy, sensitivity, specificity, brier_class))
  
}

roc_bt  = NULL
metrics_bt = NULL
for (i in 1:length(bt_resampled)) {
  roc_bt[[i]] = bt_resampled[[i]] %>% 
    collect_predictions(summarize = F) %>% 
    roc_curve(status, .pred_new, event_level = "first") %>% 
    mutate(combo = paste0("bt_combo", i))
  
  metrics_bt[[i]] = bt_resampled[[i]] %>% 
    collect_metrics() %>% 
    mutate(combo = paste0("bt_combo", i))
}

##########################################

bind_rows(roc_lr) %>% 
  bind_rows(roc_lda) %>% 
  bind_rows(roc_rf) %>% 
  bind_rows(roc_bt) %>% 
  separate(combo, into = c("model", "combo")) %>% 
  ggplot(aes(x =1-specificity, y = sensitivity, group = interaction(combo, model), colour = model ))+
  geom_path()+facet_wrap(~combo)+
  geom_abline(linetype = 2)+
  coord_fixed()

combinations_metrics = bind_rows(metrics_lr) %>% 
  bind_rows(metrics_lda) %>% 
  bind_rows(metrics_rf) %>% 
  bind_rows(metrics_bt) %>% 
  separate(combo, into = c("model", "combo")) %>% 
  mutate(combo = str_replace(combo, "combo", "")) %>% 
  mutate(.metric = factor(.metric, levels = c("roc_auc", "sensitivity", "specificity", "pr_auc", "accuracy", "brier_class"))) %>% 
  ggplot(aes(x = combo, y = mean, group = model, colour = model))+
  geom_point(size = 3)+
  facet_wrap(~.metric, scales = "free")+
  scale_colour_manual(
    labels = c("boosted tree", "linear discriminant \nanalysis", "logistic regression", "random forest"),                            
    values = c("#ABA0F2", "#F29057", "#732514", "#F25A38"))+
  theme_linedraw()+
  theme(strip.background = element_rect(fill = NA, colour = NA),
        strip.text = element_text(colour = "black", size = 13, hjust = 0),
        axis.text.x = element_text(size = 13, colour = "white"),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.position = "right",
        plot.margin=unit(c(0,0,0,0), "pt"))+
  # guides(color = guide_legend(nrow = 2, byrow = TRUE))+
  labs(y = "cross-validated mean", x = "", colour = "")
# svg(here::here("figures/combo_metric_comparison-plasmodium-macaquemodel.svg"), width = 8, height = 5.5)
# combinations_metrics
# dev.off()


single_roc_recent = read_csv(here::here("figures/single_roc_plasmodium.csv")) %>% filter(protein=="PkMSP10") %>% 
  select(sensitivity, specificity, combo = protein)

roc_rf_recent = bind_rows(roc_rf)%>% 
  select(-.threshold) %>% 
  bind_rows(single_roc_recent) %>% 
  mutate(line = case_when(
    combo=="PkMSP10" ~ "single antigen",
    T~ "multi-antigen\nrandom forest"
  )) %>% arrange(desc(combo)) %>% 
  mutate(combo = if_else(combo=="rf_combo1", "macaque_adapted", combo) ) %>% 
  mutate(combo = factor(combo, levels = c("macaque_adapted", "PkMSP10"))) %>% 
  ggplot(aes(x = 1-specificity, y = sensitivity, colour = fct_rev(combo), group = combo, linetype = line))+
  geom_path(size = 1.5, alpha = 0.8)+
  geom_abline(linetype = 2)+
  scale_x_continuous(limits = c(0,1), breaks = c(0.2, 0.4, 0.6, 0.8,1), expand = c(0,0))+
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.2, 0.4, 0.6, 0.8,1), expand = c(0,0))+
  scale_colour_manual(values = c("darkgrey", "#45858C" ),
                      labels = c("PkMSP10",expression(italic("macaque_adapted"))))+
  theme_linedraw()+
  labs(linetype = "")+
  theme(
    panel.border = element_rect(colour = "grey", fill = NA),
    legend.position = "inside",
    legend.position.inside = c(0.75, 0.26),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 12),
    legend.text = element_text( size = 12),
    plot.tag = element_text(size = 22))+
  coord_equal()+labs(colour = "")


#################################################################################################
# extract model 1-for macaque

rf_resampled[[1]]

##----------get the predictions from the cross-validation
class = rf_resampled[[1]] %>% 
  collect_predictions(summarize = T) %>% arrange(.row)

##----------get the metrics for AUC and se values
auc = rf_resampled[[1]] %>% 
  collect_metrics() %>% 
  filter(.metric=="roc_auc")

##----------make a table of threshold values, sens, spec and add se auc and a model name to use later on
roc = roc_rf[[1]] %>% 
  mutate(se = auc$std_err,
         auc = auc$mean,
         model = paste0("rf_combo_macaque"))

thresholds = roc %>% 
  get_thresholds() %>% 
  arrange(threshold_value)

model_macaque = fit(rf_workflow[[1]], training_macaque)

rf_model_macaque = list(roc, thresholds, model_macaque)

saveRDS(rf_model_macaque, here::here("models/rf_ranger_combo_macaque_plasmodium.Rds"))

rf_model_macaque = readRDS(here::here("models/rf_ranger_combo_macaque_plasmodium.Rds"))

####################################
source(here::here("scripts-analysis/2_data-tidy.R"))
set.seed("123")
cross_validation_pred = class %>% 
  mutate(sample_id = data_id) %>% 
  left_join(cats_for_plot %>% 
              dplyr::select(sample_id, day_cat, sample_cat) %>% distinct()) %>% 
  mutate(sample_cat = factor(sample_cat, levels = c("knowlesi", "vivax", "ARC", "VBDR", "TRC", "SNC"))) %>% 
  ggplot(aes(x = interaction(day_cat), y = .pred_new, colour = sample_cat))+
  geom_jitter(width = 0.35)+
  geom_boxplot(position = position_dodge2(preserve = "single"),outlier.shape = NA, alpha =0.6)+
  scale_colour_manual(values = c('#F24141',  "#F2B705", "#5E9ABF", "#768FA6", "#B3D5F2", "#88E8F2"),
                      labels = c(expression(italic("P. knowlesi")), expression(italic("P. vivax")), "ARC", 
                                 "VBDR", "TRC", "SNC"),
                      guide = guide_legend(override.aes = list(fill = c('#F24141',  "#F2B705", "#5E9ABF", "#768FA6", "#B3D5F2", "#88E8F2"),
                                                               linetype = c(0,0,0,0,0,0),
                                                               alpha = c(1,1,1,1,1,1))))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  theme_minimal()+
  labs(y = "cross-validation prediction value", fill = "", x = "sample follow-up time", colour = "",
       title = c(expression(paste("Either ", italic("Plasmodium"), " random forest: key macaque antigens")))
  )+
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.tag = element_text(size = 22),
        plot.title = element_text(size = 14),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        panel.border = element_rect(colour = "black",fill = NA),
        strip.text = element_text(size = 14, hjust = 0))




### predict on the data not included in model - the acute infections---------


predict_current = classification %>% filter(day==0) %>%
  dplyr::select(sample_id, protein_label_clean, rau, sample_cat) %>%
  pivot_wider(names_from = protein_label_clean, values_from = rau)

predict_current %>% count(sample_cat)

# 
current_classified = unknown_predict(rf_model_macaque, predict_current)
# 
set.seed(123)
acute_preds = current_classified %>%
  mutate(prefix = str_extract(sample_id, "[A-Z]+") ) %>%
  ggplot(aes(x = sample_cat, y = .pred_new, colour=sample_cat))+
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA)+
  geom_jitter(width = 0.1, alpha =0.6, size = 0.7)+
  theme_minimal()+
  scale_colour_manual(values = c('#F24141',  "#F2B705"),
                      labels = c(expression(italic("P. knowlesi")), expression(italic("P. vivax"))), 
                      guide = guide_legend(override.aes = list(fill = c('#F24141',  "#F2B705"),
                                                               linetype = c(0,0),
                                                               alpha = c(1,1))))+
  theme(axis.text.x = element_text(size = 12, colour = NA),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.tag = element_text(size = 22),
        plot.title = element_text(size = 12),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        panel.border = element_rect(colour = "black",fill = NA),
        strip.text = element_text(size = 14, hjust = 0))+
  labs(y = "prediction value for\nacute infections", x = "", colour = "acute\ninfection",
       title = c(expression(paste("Either ", italic("Plasmodium"), ": macaque model"))),
       subtitle = "predictions for acute infections")



### plot supp figure 12

design <- "
AAAABBBB
AAAABBBB
CCCCCDDD
"


svg(here::here("figures/figure_s12.svg"),width = 14, height = 10)
print(combinations_metrics + roc_rf_recent+
        cross_validation_pred + acute_preds+
        patchwork::plot_annotation(tag_levels = "A")+
        patchwork::plot_layout(design = design )& 
        
        theme(plot.tag = element_text(size = 22))
)
dev.off()

