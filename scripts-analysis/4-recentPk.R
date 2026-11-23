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

data = classification_fn(classification, classify_month_Pk) %>% 
  dplyr::select(sample_id, status, proteins_paper_clean)

###----------REMOVE ACUTE INFECTIONS
data = data %>%
  filter(!str_detect(sample_id, "-0")) %>%
  filter(!str_detect(sample_id, "w0")) %>% 
  mutate(status = factor(status, levels = c("new", "old"))) 

data_id = data %>% pull(sample_id)

training  = data %>% dplyr::select(-sample_id)



################################################################################
##       single antigen ROC curves - ALL DATA 
################################################################################

roc_single_split  = data %>% 
  pivot_longer(-c(sample_id, status)) %>% 
  group_split(name)


roc_single = roc_single_split %>% 
  map(~
        pROC::roc(data = .x,
                  response = "status", 
                  predictor = "value", 
                  levels = c("old", "new"),
                  auc = T,
                  ci = T)
  ) 

all_data = NULL
for (i in 1:length(roc_single)) {
  all_data[[i]]  = cbind(sensitivity = roc_single[[i]]$sensitivities,
                         specificity = roc_single[[i]]$specificities,
                         protein = unique(roc_single_split[[i]]$name),
                         auc = unique(roc_single[[i]]$auc),
                         min_c =min(roc_single[[i]]$ci),
                         max_c =max(roc_single[[i]]$ci)) %>% as_tibble()
  
  all_data[[i]] = all_data[[i]] %>% 
    dplyr::mutate(across(c(sensitivity, specificity, auc, min_c, max_c), as.numeric))
  
}

roc_single_very_recent_pk = bind_rows(all_data) %>% 
  mutate(species = case_when(
    str_detect(protein, "Pk") ~ "Pk",
    str_detect(protein, "Pv") ~ "Pv"
  )) %>% 
  ggplot(aes(x = 1-specificity, y = sensitivity, colour = species))+
  geom_path(aes(group = protein), alpha = 0.6, size = 1.1)+
  geom_abline(linetype = 2, colour = "darkgrey")+
  scale_color_manual(values = c("#DB4F48","#E8B840"))+
  scale_x_continuous(limits = c(0,1), breaks = c(0.2, 0.4, 0.6, 0.8,1), expand = c(0,0))+
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.2, 0.4, 0.6, 0.8,1), expand = c(0,0))+
  coord_equal()+
  theme_minimal()+
  theme(legend.position = "inside",
        legend.position.inside = c(0.75, 0.25),
        panel.border = element_rect(colour = "grey", fill = NA),
        strip.text = element_text(hjust = 0, size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.text = element_text(face = "italic", size = 14),
        legend.title = element_text(size = 16),
        legend.margin=margin(-10, 0, 0, 0),
        plot.tag = element_text(size = 22))+
  labs(colour = "species of \nrecombinant\nprotein", tag = "A")



bind_rows(all_data) %>% 
  mutate(species = case_when(
    str_detect(protein, "Pk") ~ "Pk",
    str_detect(protein, "Pv") ~ "Pv"
  )) %>% 
  mutate(scenario = "pk-recent") %>% 
  write_csv(here::here("figures/single_roc_pk-recent.csv"))



### ---two by two LDA
vec = colnames(training)[-1]
two = gtools::combinations(n = 21, r = 2, v = vec)

two_list = NULL
for (i in 1:nrow(two)) {
  two_list[[i]] = training %>% 
    dplyr::select(status, two[i,1], two[i,2]) %>% 
    dplyr::mutate(combo = paste0(two[i,1]," & ", two[i,2])) %>% 
    dplyr::mutate(no_combo = "2")
}

library(discrim)
##----------DEFINE MODEL
set.seed(1234)
lda = discrim_linear() %>%
  set_engine("MASS",
             seed = 1234) %>% #sets seed
  set_mode("classification")

rec_list = NULL
roc = NULL
class = NULL
metrics = NULL

for (i in 1:length(two_list)) {
  recipe = recipe(status~., data = two_list[[i]]) %>%
    update_role(combo, no_combo, new_role = "id variable") %>%
    step_log(all_numeric_predictors(), base = 10)
  
  workflow = workflow() %>%
    add_model(lda) %>%
    add_recipe(recipe)
  
  rec_list[[i]] = workflow %>%
    fit_resamples(resamples = vfold_cv(two_list[[i]], v = 10, repeats = 5, strata = status),
                  control = control_resamples(save_pred = T,
                                              save_workflow = T),
                  metrics = metric_set(roc_auc, pr_auc, accuracy, sensitivity, specificity, brier_class))
  
  roc[[i]]  = rec_list[[i]]%>%
    collect_predictions(summarize = F) %>%
    roc_curve(status, .pred_new, event_level = "first") %>%
    mutate(combo = paste0(unique(two_list[[i]]$combo)),
           no_combo = paste0(unique(two_list[[i]]$no_combo)))
  
  class[[i]]  = rec_list[[i]]%>%
    collect_predictions(summarize = T)
  
  metrics[[i]] = rec_list[[i]] %>%
    collect_metrics()  %>%
    dplyr::mutate(combo = paste0(unique(two_list[[i]]$combo)),
                  no_combo = paste0(unique(two_list[[i]]$no_combo)))
}


proteins = proteins %>% mutate(prot_label2 = janitor::make_clean_names(protein_label_clean))

order_1 = proteins %>% dplyr::select(prot1 = protein_label_clean , order_1) 
order_2 = proteins %>% dplyr::select(prot2 = protein_label_clean , order_2)


two_by_two_plot = bind_rows(metrics) %>% 
  filter(no_combo=="2") %>% 
  filter(.metric=="roc_auc") %>% 
  dplyr::select(mean, combo) %>% 
  distinct() %>% 
  separate(combo, into = c("protein1", "protein2"), sep = " & ") %>% 
  mutate(c1 = paste0(protein1, "&", protein2),
         c2 = paste0(protein2, "&", protein1)) %>% 
  dplyr::select(mean, c1, c2) %>% 
  pivot_longer(-mean) %>% 
  separate(value, into = c("prot1", "prot2"), sep = "&") %>% 
  left_join(order_1) %>% 
  left_join(order_2) %>% distinct() %>% 
  filter(order_2>=order_1) %>% 
  ggplot(aes(x = fct_reorder(prot1, order_1), y = fct_reorder(prot2, -order_2), fill = mean))+
  geom_tile()+
  scale_fill_viridis_c(option = "plasma", limits = c(0.55, 1))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 16),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.85),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 20),
        plot.tag = element_text(size = 22),
        plot.title = element_text(size = 16),
        plot.margin=unit(c(0,0,0,0), "pt"))+
  coord_equal()+
  labs(fill = "resampled\nmean\nAUC",
       x = "", y = "", tag = "B",
       title = "")

# pdf(here::here("figures/fig3b.pdf"), width = 6, height = 6)
# two_by_two_plot
# dev.off()
# 
# jpeg(here::here("figures/fig3b.jpeg"), width = 6, height = 6,units = "in", res = 400)
# two_by_two_plot
# dev.off()


## check range of AUC values for combos with PkMSP10
bind_rows(metrics) %>% 
  filter(no_combo=="2") %>% 
  filter(.metric=="roc_auc") %>% arrange(desc(mean)) %>% filter(str_detect(combo, "PkMSP10")) %>% print(n = 240)


###-------------------- classification models with very-recently --------------
source(here::here("scripts-analysis/0_setup-classification.R"))



##-----------RECIPE
## Recipe is using factor STATUS 

training_1 = training
training_2 = training %>% dplyr::select(status,
                                        PkMSP10,
                                        PvEBP,
                                        PvRBP2b,
                                        PvMSP8,
                                        `PvMSP1-19`,
                                        `PvFAMA-95`,
                                        PvMSP5
)

training_3 = training %>% dplyr::select(-PkMSP10)

training_4 = training %>% dplyr::select(status,
                                        PkP12,
                                        PvEBP,
                                        PvRBP2b,
                                        PvMSP8,
                                        `PvMSP1-19`,
                                        `PvFAMA-95`,
                                        PvMSP5)

training_5 = training %>% dplyr::select(status,
                                        PvEBP,
                                        PvRBP2b,
                                        PvMSP8,
                                        `PvMSP1-19`,
                                        `PvFAMA-95`,
                                        PvMSP5)



recipe_model1 = recipe(status~., data = training) %>% 
  step_log(all_numeric())

recipe_model2 = recipe(status~., data = training_2) %>% 
  step_log(all_numeric())

recipe_model3 = recipe(status~., data = training_3) %>% 
  step_log(all_numeric())

recipe_model4 = recipe(status~., data = training_4) %>% 
  step_log(all_numeric())

recipe_model5 = recipe(status~., data = training_5) %>% 
  step_log(all_numeric())

training_list = list(training_1, training_2, training_3,  training_4,training_5)
recipe_list = list(recipe_model1, recipe_model2, recipe_model3,  recipe_model4, recipe_model5)

## - model evaluation with workflowsets

# take forward all 5 models

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
  ggplot(aes(x = combo, y = mean, group = model, colour = model, shape = model, size = model))+
  geom_point()+
  facet_wrap(~.metric, scales = "free")+
  scale_colour_manual(
    labels = c("boosted tree", "linear discriminant \nanalysis", "logistic regression", "random forest"),                            
    values = c("#ABA0F2", "#F29057", "#732514", "#F25A38"))+
  scale_shape_manual(values = c(18, 16, 15, 17),
                     labels = c("boosted tree", "linear discriminant \nanalysis", "logistic regression", "random forest"))+
  scale_size_manual(values = c(5,4,5,5),
                    labels = c("boosted tree", "linear discriminant \nanalysis", "logistic regression", "random forest"))+
  theme_linedraw()+
  theme(strip.background = element_rect(fill = NA, colour = NA),
        strip.text = element_text(colour = "black", size = 12, hjust = 0),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 11))+
  labs(y = "cross-validated mean", x = "antigen combination")

svg(here::here("figures/figure_s4_recent.svg"), width = 12, height = 5.5)
combinations_metrics
dev.off()


single_roc_recent = read_csv(here::here("figures/single_roc_pk-recent.csv")) %>% filter(protein=="PkMSP10") %>% 
  select(sensitivity, specificity, combo = protein)

roc_rf_recent = bind_rows(roc_rf)%>% 
  select(-.threshold) %>% 
  bind_rows(single_roc_recent) %>% 
  mutate(line = case_when(
    combo=="PkMSP10" ~ "single antigen",
    T~ "multi-antigen\nrandom forest"
  )) %>% arrange(desc(combo)) %>% 
  mutate(combo = factor(combo, levels = c("rf_combo6", "rf_combo5", "rf_combo4", "rf_combo3", "rf_combo2", "rf_combo1", "PkMSP10"))) %>% 
  ggplot(aes(x = 1-specificity, y = sensitivity, colour = fct_rev(combo), group = combo, linetype = line))+
  geom_path(size = 1.5, alpha = 0.8)+
  geom_abline(linetype = 2)+
  scale_x_continuous(limits = c(0,1), breaks = c(0.2, 0.4, 0.6, 0.8,1), expand = c(0,0))+
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.2, 0.4, 0.6, 0.8,1), expand = c(0,0))+
  scale_colour_manual(values = c("darkgrey", "#B0D1D9","#45858C",  "#467302",  "#A62626", "#F2B705"))+
  theme_linedraw()+
  theme(
    panel.border = element_rect(colour = "grey", fill = NA),
    axis.title = element_text(size = 16),
    legend.position = "inside",
    legend.position.inside = c(0.82, 0.35),
    axis.text = element_text(size = 14),
    legend.text = element_text( size = 14),
    plot.tag = element_text(size = 22))+
  coord_equal()+labs(colour = "", linetype ="", tag = "C")

# svg(here::here("figures/roc_rf_recent.svg"), width = 8, height = 6)
# roc_rf_recent
# dev.off


svg(here::here("figures/fig3-pkrecent.svg"), width = 13, height = 13)
patchwork::wrap_plots(roc_single_very_recent_pk, two_by_two_plot, roc_rf_recent
) + patchwork::plot_layout(nrow = 2)&theme(plot.tag = element_text(size = 25))
dev.off()



jpeg(here::here("figures/fig3-pkrecent.jpeg"), width = 13, height = 13, units = "in", res = 500)
patchwork::wrap_plots(roc_single_very_recent_pk, two_by_two_plot, roc_rf_recent
) + patchwork::plot_layout(nrow = 2) &theme(plot.tag = element_text(size = 25))
dev.off()


#################################################################################################
# extract model 2

rf_resampled[[2]]

##----------get the predictions from the cross-validation
class = rf_resampled[[2]] %>% 
  collect_predictions(summarize = T) %>% arrange(.row)

##----------get the metrics for AUC and se values
auc = rf_resampled[[2]] %>% 
  collect_metrics() %>% 
  filter(.metric=="roc_auc")

##----------make a table of threshold values, sens, spec and add se auc and a model name to use later on
roc = roc_rf[[2]] %>% 
  mutate(se = auc$std_err,
         auc = auc$mean,
         model = paste0("rf_combo2"))

thresholds = roc %>% 
  get_thresholds() %>% 
  arrange(threshold_value)

set.seed(1234)
model_2 = fit(rf_workflow[[2]], training_2)

rf_model2 = list(roc, thresholds, model_2)

saveRDS(rf_model2, here::here("models/rf_ranger_combo2_very-recent.Rds"))

rf_ranger = readRDS(here::here("models/rf_ranger_combo2_very-recent.Rds"))

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
       title = c(expression(paste("Recent ", italic("P. knowlesi"), " random forest")))
  )+
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.tag = element_text(size = 22),
        plot.title = element_text(size = 14),
        legend.position = "right",
        legend.text = element_text(size = 12),
        panel.border = element_rect(colour = "black",fill = NA),
        strip.text = element_text(size = 14, hjust = 0))


# svg(here::here("figures/supp-cross-validation-pred-new-pk.svg"),width = 9, height = 4)
# print(cross_validation_pred)
# dev.off()

## also save ggplot as rds so can merge for combining as svg

saveRDS(cross_validation_pred, file = here::here("figures/rds/cross_validation_pred-recentpk.rds"))


### predict on the data not included in model - the acute infections---------

predict_current = classification %>% filter(day==0) %>%
  dplyr::select(sample_id, protein_label_clean, rau, sample_cat) %>%
  pivot_wider(names_from = protein_label_clean, values_from = rau)

predict_current %>% count(sample_cat)

# 
current_classified = unknown_predict(rf_ranger, predict_current)
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
        legend.position = "right",
        legend.text = element_text(size = 12),
        panel.border = element_rect(colour = "black",fill = NA),
        strip.text = element_text(size = 14, hjust = 0))+
  labs(y = "prediction value for acute infections", x = "", colour = "acute\ninfection",
       title = c(expression(paste("Recent ", italic("P. knowlesi")))),
       subtitle = "predictions for acute infections")


# svg(here::here("figures/acute-pred-recent.svg"),width = 4, height = 4)
# print(acute_preds)
# dev.off()

saveRDS(acute_preds, file = here::here("figures/rds/acute_preds-recentpk.rds"))
###################################################################################


##------------------------------- setup data where a test set and training set have NO ids in common. 
partition = classification %>% 
  filter(day!=0|is.na(day)) %>%
  dplyr::select(sample_id, protein_label_clean, rau, sample_cat, id, classify_month_Pk) %>%
  pivot_wider(names_from = protein_label_clean, values_from = rau) %>% mutate(status = as.factor(classify_month_Pk))


repeat_ids = partition %>% pull(id) %>% unique()

choose_holdout_n = round(length(repeat_ids)*0.2)
set.seed(123)
holdout_ids = sample(repeat_ids,size = choose_holdout_n)


training_partition = partition %>% filter(!id%in%holdout_ids) %>% dplyr::select(sample_id, status, PkMSP10, `PvFAMA-95`, PvRBP2b, `PvMSP1-19`, PvEBP, PvMSP5, PvMSP8)
holdout_partition = partition %>% filter(id%in%holdout_ids) %>% dplyr::select(sample_id,status, PkMSP10, `PvFAMA-95`, PvRBP2b, `PvMSP1-19`, PvEBP, PvMSP5, PvMSP8)


### model
rf = rand_forest(mode = "classification", 
                 # mtry = tune(),
                 # min_n = tune(),
                 trees = 5000) %>% 
  set_engine("ranger", importance = "permutation")

training_partition_ids = training_partition %>% pull(sample_id)
training_partition = training_partition %>% dplyr::select(-sample_id)

recipe_model_holdout = recipe(status~., data = training_partition) %>% 
  # update_role(sample_id, new_role = "id") %>% 
  step_log(all_numeric())

rf_workflow_holdout = workflow() %>% add_model(rf) %>% add_recipe(recipe_model_holdout)

rf_resampled_holdout = rf_workflow_holdout %>% fit_resamples(resamples = vfold_cv(training_partition, v = 10, repeats = 5, strata = status),
                                                             control = control_resamples(save_pred = T,
                                                                                         save_workflow = T),
                                                             metrics = metric_set(roc_auc, pr_auc, accuracy, sensitivity, specificity, brier_class))


roc_rf_holdout= rf_resampled_holdout %>% 
  collect_predictions(summarize = F) %>% 
  roc_curve(status, .pred_new, event_level = "first") %>% 
  mutate(combo = paste0("rf_combo2_holdout"))

metrics_rf_holdout = rf_resampled_holdout%>% 
  collect_metrics() %>% 
  mutate(combo = paste0("rf_combo2_holdout"))


class_holdout = rf_resampled_holdout %>% 
  collect_predictions(summarize = T) %>% arrange(.row)

##----------get the metrics for AUC and se values
auc = rf_resampled_holdout %>% 
  collect_metrics() %>% 
  filter(.metric=="roc_auc")

##----------make a table of threshold values, sens, spec and add se auc and a model name to use later on
roc = roc_rf_holdout%>% 
  mutate(se = auc$std_err,
         auc = auc$mean,
         model = paste0("rf_combo2_holdout"))

thresholds = roc %>% 
  get_thresholds() %>% 
  arrange(threshold_value)

model_2_holdout = fit(rf_workflow_holdout, training_partition)

rf_model2_holdout = list(roc, thresholds, model_2_holdout)


### now predict on the data that doesn't have any of the same ids...

predict_holdout = unknown_predict(rf_model2_holdout, holdout_partition)

set.seed(123)
holdout_plot = predict_holdout %>% 
  left_join(classification) %>% 
  mutate(sample_cat = factor(sample_cat, levels = c("knowlesi", "vivax", "ARC", "VBDR", "TRC", "SNC"))) %>% 
  dplyr::select(day_cat,sample_id, .pred_new, status, sample_cat ) %>% distinct() %>% 
  ggplot(aes(x = day_cat, y = .pred_new, colour = sample_cat))+
  geom_jitter(width = 0.3)+
  geom_boxplot(position = position_dodge2(preserve = "single"),outlier.shape = NA, alpha =0.6) +
  scale_colour_manual(values = c('#F24141',  "#F2B705", "#5E9ABF", "#768FA6", "#B3D5F2", "#88E8F2"),
                      labels = c(expression(italic("P. knowlesi")), expression(italic("P. vivax")), "ARC", 
                                 "VBDR", "TRC", "SNC"),
                      guide = guide_legend(override.aes = list(fill = c('#F24141',  "#F2B705", "#5E9ABF", "#768FA6", "#B3D5F2", "#88E8F2"),
                                                               linetype = c(0,0,0,0,0,0),
                                                               alpha = c(1,1,1,1,1,1))))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  theme_minimal()+
  labs(y = "prediction value", fill = "", x = "sample follow-up time", colour = "",
       title = c(expression(paste("Recent ", italic("P. knowlesi"), " random forest"))),
       subtitle = "Holdout Unique IDs"
  )+
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.tag = element_text(size = 22),
        plot.title = element_text(size = 14),
        legend.position = "right",
        legend.text = element_text(size = 12),
        panel.border = element_rect(colour = "black",fill = NA),
        strip.text = element_text(size = 14, hjust = 0))

# svg(here::here("figures/holdout-pred-new-pk.svg"),width = 9, height = 4)
# print(holdout_plot)
# dev.off()

### save rds
saveRDS(holdout_plot, file = here::here("figures/rds/holdout_plot-recentpk.rds"))

### correlation of holdout model predicted vs non-holdout normal

correlation_pred_new = predict_holdout %>% dplyr::select(sample_id, .pred_new_holdout = .pred_new) %>% 
  left_join(class %>% 
              mutate(sample_id = data_id) %>% 
              left_join(cats_for_plot %>% 
                          dplyr::select(sample_id, day_cat, sample_cat) %>% distinct()) ) %>%
  mutate(sample_cat = factor(sample_cat, levels = c("knowlesi", "vivax", "ARC", "VBDR", "TRC", "SNC"))) %>% 
  ggplot(aes(x = .pred_new_holdout, y = .pred_new, colour = sample_cat))+
  geom_abline(linetype = 2, colour = "darkgrey")+
  geom_point(size = 2)+
  scale_colour_manual(values = c('#F24141',  "#F2B705", "#5E9ABF", "#768FA6", "#B3D5F2", "#88E8F2"),
                      labels = c(expression(italic("P. knowlesi")), expression(italic("P. vivax")), "ARC", 
                                 "VBDR", "TRC", "SNC"),
                      guide = guide_legend(override.aes = list(fill = c('#F24141',  "#F2B705", "#5E9ABF", "#768FA6", "#B3D5F2", "#88E8F2"),
                                                               linetype = c(0,0,0,0,0,0),
                                                               alpha = c(1,1,1,1,1,1))))+
  theme_minimal()+
  scale_x_continuous(breaks = c(0.25, 0.5, 0.75, 1))+
  labs(colour = "sample category")+
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1.2),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        plot.tag = element_text(size = 22),
        plot.title = element_text(size = 14),
        legend.position = "right",
        legend.text = element_text(size = 12),
        panel.border = element_rect(colour = "black",fill = NA),
        strip.text = element_text(size = 12, hjust = 0)) +facet_wrap(~day_cat, nrow = 2)+
  coord_equal()+
  geom_abline(linetype = 2, colour = "darkgrey")+
  labs(y = "cross-validation prediction values", x = "prediction values on holdout")

# svg(here::here("figures/corr-holdout-pred-new-pk.svg"),width = 12, height = 5)
# print(correlation_pred_new)
# dev.off()


## save rds

saveRDS(correlation_pred_new, file = here::here("figures/rds/correlation_pred_new-recentpk.rds"))


### use 'holdout model' to predict on acute infections and restrict to sample_ids not in training

predict_holdout_acute = unknown_predict(rf_model2_holdout, predict_current)

### is there the same pattern correlating predictions on the acute infections 

current_all_comparison = predict_holdout_acute %>% left_join(classification ) %>% 
  select(sample_id, id, .pred_new, day_cat, sample_cat) %>% 
  filter(!id%in%holdout_ids | str_detect(id, "QDM")) %>% 
  dplyr::rename(.pred_new_holdout = .pred_new) %>% 
  left_join(current_classified) %>% dplyr::rename(.pred_new_original = .pred_new) %>% 
  ggplot(aes(y=.pred_new_original, x =  .pred_new_holdout, colour = sample_cat))+
  geom_point(size = 3, alpha = 0.7)+
  geom_abline(linetype = 2)+
  facet_wrap(~day_cat)+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1.2),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        plot.tag = element_text(size = 22),
        plot.title = element_text(size = 14),
        legend.position = "right",
        legend.text = element_text(size = 12),
        panel.border = element_rect(colour = "black",fill = NA),
        strip.text = element_text(size = 12, hjust = 0)) +facet_wrap(~day_cat, nrow = 2)+
  coord_equal()+
  scale_colour_manual(values = c('#F24141',  "#F2B705"),
                      labels = c(expression(italic("P. knowlesi")), expression(italic("P. vivax"))),
                      guide = guide_legend(override.aes = list(fill = c('#F24141',  "#F2B705"),
                                                               linetype = c(0,0),
                                                               alpha = c(1,1))))


# svg(here::here("figures/corr-holdout-pred-new-pk_current.svg"),width = 4, height = 4)
# print(current_all_comparison)
# dev.off()




########################################################

### ---------- repeat trial with one sample per individual

set.seed(123)
samples_no_repeat <- classification %>% 
  filter(day!=0|is.na(day)) %>%
  dplyr::select(sample_id, protein_label_clean, rau, sample_cat, id, classify_month_Pk, day) %>%
  pivot_wider(names_from = protein_label_clean, values_from = rau) %>% mutate(status = as.factor(classify_month_Pk))

ids_no_repeat = samples_no_repeat %>% 
  dplyr::select(id, day) %>% 
  distinct() %>% 
  group_by(id) %>% 
  sample_n(1) %>% 
  mutate(sample_id = case_when(
    is.na(day) ~ id,
    TRUE ~ paste0(id, "-", day)))

classification_downselected = samples_no_repeat %>% 
  mutate(sample_id_filter = case_when(
    is.na(day) ~ id,
    TRUE ~ paste0(id, "-", day))) %>% 
  filter(sample_id_filter%in%c(ids_no_repeat$sample_id)) %>% 
  dplyr::select(sample_id, status, PkMSP10, `PvFAMA-95`, PvRBP2b, `PvMSP1-19`, PvEBP, PvMSP5, PvMSP8)

ids_downselected = classification_downselected %>% dplyr::select(sample_id)
training_downselected = classification_downselected %>% dplyr::select(-sample_id)

### new recipe with classification_downselected

### model
rf = rand_forest(mode = "classification", 
                 # mtry = tune(),
                 # min_n = tune(),
                 trees = 5000) %>% 
  set_engine("ranger", importance = "permutation")

recipe_model_downselect = recipe(status~., data = training_downselected) %>% 
  # update_role(sample_id, new_role = "id") %>% 
  step_log(all_numeric())

rf_workflow_downselect = workflow() %>% add_model(rf) %>% add_recipe(recipe_model_downselect)

rf_resampled_downselect = rf_workflow_downselect %>% fit_resamples(resamples = vfold_cv(training_downselected, v = 10, repeats = 5, strata = status),
                                                                   control = control_resamples(save_pred = T,
                                                                                               save_workflow = T),
                                                                   metrics = metric_set(roc_auc, pr_auc, accuracy, sensitivity, specificity, brier_class))


roc_rf_downselect= rf_resampled_downselect %>% 
  collect_predictions(summarize = F) %>% 
  roc_curve(status, .pred_new, event_level = "first") %>% 
  mutate(combo = paste0("rf_combo2_downselect"))

metrics_rf_downselect = rf_resampled_downselect%>% 
  collect_metrics() %>% 
  mutate(combo = paste0("rf_combo2_downselect"))


class_downselect = rf_resampled_downselect %>% 
  collect_predictions(summarize = T) %>% arrange(.row)

##----------get the metrics for AUC and se values
auc = rf_resampled_downselect %>% 
  collect_metrics() %>% 
  filter(.metric=="roc_auc")

##----------make a table of threshold values, sens, spec and add se auc and a model name to use later on
roc = roc_rf_downselect%>% 
  mutate(se = auc$std_err,
         auc = auc$mean,
         model = paste0("rf_combo2_downselect"))

thresholds = roc %>% 
  get_thresholds() %>% 
  arrange(threshold_value)

model_2_downselect = fit(rf_workflow_downselect, training_downselected)

rf_model2_downselect = list(roc, thresholds, model_2_downselect)


##########

training_reduced_ind = class_downselect %>% 
  dplyr::mutate(ids_downselected) %>% 
  left_join(classification) %>% 
  mutate(sample_cat = factor(sample_cat, levels = c("knowlesi", "vivax", "ARC", "VBDR", "TRC", "SNC"))) %>% 
  dplyr::select(day_cat,sample_id, .pred_new, status, sample_cat ) %>% distinct() %>% 
  ggplot(aes(x = day_cat, y = .pred_new, colour = sample_cat))+
  geom_jitter(width = 0.3)+
  geom_boxplot(position = position_dodge2(preserve = "single"),outlier.shape = NA, alpha =0.6) +
  scale_colour_manual(values = c('#F24141',  "#F2B705", "#5E9ABF", "#768FA6", "#B3D5F2", "#88E8F2"),
                      labels = c(expression(italic("P. knowlesi")), expression(italic("P. vivax")), "ARC", 
                                 "VBDR", "TRC", "SNC"),
                      guide = guide_legend(override.aes = list(fill = c('#F24141',  "#F2B705", "#5E9ABF", "#768FA6", "#B3D5F2", "#88E8F2"),
                                                               linetype = c(0,0,0,0,0,0),
                                                               alpha = c(1,1,1,1,1,1))))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  theme_minimal()+
  labs(y = "prediction value", fill = "", x = "sample follow-up time", colour = "",
       title = c(expression(paste("Recent ", italic("P. knowlesi"), " random forest"))),
       subtitle = ""
  )+
  theme(axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.tag = element_text(size = 22),
        plot.title = element_text(size = 14),
        legend.position = "right",
        legend.text = element_text(size = 12),
        panel.border = element_rect(colour = "black",fill = NA),
        strip.text = element_text(size = 14, hjust = 0))


# svg(here::here("figures/training_reduced_individuals_pred_new.svg"),width = 12, height = 5)
# print(training_reduced_ind)
# dev.off()


## also save ggplot as rds so can merge for combining as svg

saveRDS(training_reduced_ind, file = here::here("figures/rds/training_reduced_ind_pk.rds"))
