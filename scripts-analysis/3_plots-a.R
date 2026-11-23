library(tidymodels)
source(here::here("scripts-analysis/0_setup.R"))
source(here::here("scripts-analysis/2_data-tidy.R"))


dir.create(here::here("figures"))
#############################
############ RAU PLOTS
##############################

plot_lists_rau <- cats_for_plot_list %>% 
  map(~
        ggplot( data = .x, aes (x = day_cat_short, y = rau, fill = sample_cat))+
        geom_boxplot(position = position_dodge(preserve = "single"), outlier.size = 0.3, size = 0.5, width = 1.2)+
        geom_hline(yintercept = 0.02, linetype=2, size = 0.2)+
        geom_hline(yintercept = 1.95e-05, linetype=2, size = 0.2)+
        facet_wrap(~protein_label_clean, scales = "free")+
        scale_y_log10(labels = scales::comma)+
        scale_fill_manual(values = c('#F24141',  "#F2B705", "#5E9ABF", "#768FA6", "#B3D5F2", "#88E8F2"),
                          labels = c(expression(italic("P. knowlesi")), expression(italic("P. vivax")), "ARC", 
                                     "VBDR", "TRC", "SNC"),
                          guide = guide_legend(override.aes = list(linetype = c(0,0,0,0,0,0))))+
        theme_minimal()+
        scale_x_discrete(labels = c("current", "1 wk", "2 wk", "1 mo", "3 mo", "6 mo", "9 mo+", "controls"))+
        labs(y = "RAU", fill = "")+
        theme(axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 14),
              plot.tag = element_text(size = 22),
              legend.position = "right",
              panel.border = element_rect(colour = "black",fill = NA),
              strip.text = element_text(size = 14, hjust = 0)))


RAU = (plot_lists_rau$PkMSP10 +labs(tag = "A", y = "RAU")) + 
  (plot_lists_rau$PvMSP10 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())) + 
  (plot_lists_rau$PvEBP + theme(axis.text.y = element_blank(), axis.title.y = element_blank()))+ 
  (plot_lists_rau$PvRBP2b + theme(axis.text.y = element_blank(), axis.title.y = element_blank()))+ 
  (plot_lists_rau$PvEXP1 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())) +
  plot_spacer() + 
  
  (plot_lists_rau$`PkPTRAMP-CSS` + labs(tag = "B", y = "RAU")) +
  (plot_lists_rau$PkRIPR+ theme(axis.text.y = element_blank(), axis.title.y = element_blank())) +
  (plot_lists_rau$PkP12 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())) +   
  (plot_lists_rau$Pksera3ag1+ theme(axis.text.y = element_blank(), axis.title.y = element_blank())) +
  (plot_lists_rau$Pksera3ag2+ theme(axis.text.y = element_blank(), axis.title.y = element_blank())) +
  (plot_lists_rau$PkSSP2+ theme(axis.text.y = element_blank(), axis.title.y = element_blank())) +
  
  (plot_lists_rau$PkP38+labs( y = "RAU"))+
  (plot_lists_rau$PvMSP8+ theme(axis.text.y = element_blank(), axis.title.y = element_blank()))  +
  
  (plot_lists_rau$`PvMSP1-19` + theme(axis.text.y = element_blank(), axis.title.y = element_blank()))  + 
  (plot_lists_rau$`PvFAMA-95` + theme(axis.text.y = element_blank(), axis.title.y = element_blank()))  + 
  (plot_lists_rau$PvP12 + theme(axis.text.y = element_blank(), axis.title.y = element_blank()))  + 
  (plot_lists_rau$PvMSP5 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())) + 
  
  
  (plot_lists_rau$Pk_031930ag1 +labs(tag = "C", y = "RAU"))+
  (plot_lists_rau$Pk_031930ag2+ theme(axis.text.y = element_blank(), axis.title.y = element_blank()))  +
  (plot_lists_rau$Pk_0800300+ theme(axis.text.y = element_blank(), axis.title.y = element_blank())) +
  (plot_lists_rau$Pk8+ theme(axis.text.y = element_blank(), axis.title.y = element_blank())) +
  plot_spacer()+ guide_area()+
  
  plot_layout(guides = "collect", ncol = 6) & 
  theme(axis.title.x =  element_blank(), 
        legend.text.align = 0,
        legend.text = element_text(size = 18)) 

svg(here::here("figures/fig1_rau.svg"), width = 12, height = 16)
print(RAU)
dev.off()

#############################
############ MFI PLOTS
##############################

plot_lists_mfi <- cats_for_plot_list %>% 
  map(~ ggplot(data = .x, aes (x = day_cat_short, y = mfi, fill = sample_cat))+
        geom_boxplot(position = position_dodge(preserve = "single"), outlier.size = 0.3, size = 0.5, width = 1.2)+
        facet_wrap(~protein_label_clean)+
        scale_y_log10(labels = scales::comma)+
        scale_fill_manual(values = c('#F24141',  "#F2B705", "#5E9ABF", "#768FA6", "#B3D5F2", "#88E8F2"),
                          labels = c(expression(italic("P. knowlesi")), expression(italic("P. vivax")), "ARC", 
                                     "VBDR", "TRC", "SNC"),
                          guide = guide_legend(override.aes = list(linetype = c(0,0,0,0,0,0))))+
        theme_minimal()+
        scale_x_discrete(labels = c("current", "1 wk", "2 wk", "1 mo", "3 mo", "6 mo", "9 mo+", "controls"))+
        labs(y = "MFI", fill = "")+
        theme(axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 14),
              plot.tag = element_text(size = 22),
              legend.position = "right",
              panel.border = element_rect(colour = "black",fill = NA),
              strip.text = element_text(size = 14, hjust = 0)))


MFI = (plot_lists_mfi$PkMSP10 +labs(tag = "A", y = "MFI")) + 
  (plot_lists_mfi$PvMSP10 + theme( axis.title.y = element_blank())) + 
  (plot_lists_mfi$PvEBP + theme( axis.title.y = element_blank()))+ 
  (plot_lists_mfi$PvRBP2b + theme( axis.title.y = element_blank()))+ 
  (plot_lists_mfi$PvEXP1 + theme(axis.title.y = element_blank())) +
  plot_spacer() + 
  
  
  (plot_lists_mfi$`PkPTRAMP-CSS` + labs(tag = "B",y = "MFI")) +
  (plot_lists_mfi$PkRIPR+ theme( axis.title.y = element_blank())) +
  (plot_lists_mfi$PkP12 + theme( axis.title.y = element_blank())) +   
  (plot_lists_mfi$Pksera3ag1+ theme(axis.title.y = element_blank())) +
  (plot_lists_mfi$Pksera3ag2+ theme( axis.title.y = element_blank())) +
  (plot_lists_mfi$PkSSP2+ theme( axis.title.y = element_blank())) +
  
  
  (plot_lists_mfi$PkP38 +labs( y = "MFI"))+
  (plot_lists_mfi$PvMSP8+ theme( axis.title.y = element_blank()))  +
  (plot_lists_mfi$`PvMSP1-19` + theme( axis.title.y = element_blank()))  + 
  (plot_lists_mfi$`PvFAMA-95` + theme( axis.title.y = element_blank()))  + 
  (plot_lists_mfi$PvP12 + theme( axis.title.y = element_blank()))  + 
  (plot_lists_mfi$PvMSP5 + theme(axis.title.y = element_blank())) + 
  
  
  (plot_lists_mfi$Pk_031930ag1 +labs(tag = "C", y = "MFI"))+
  (plot_lists_mfi$Pk_031930ag2+ theme( axis.title.y = element_blank()))  +
  (plot_lists_mfi$Pk_0800300+ theme( axis.title.y = element_blank())) +
  (plot_lists_mfi$Pk8+ theme( axis.title.y = element_blank())) +
  plot_spacer()+ guide_area()+
  
  plot_layout(guides = "collect", ncol = 6) & 
  theme(axis.title.x =  element_blank(), 
        legend.text.align = 0,
        legend.text = element_text(size = 18)) 

svg(here::here("figures/supp_1_MFI.svg"), width = 14, height = 15)
print(MFI)
dev.off()

#############################
############ FOLD CHANGE PLOTS
##############################

fold_change_plot <- fold_change %>% 
  ggplot(aes(x = fct_reorder(protein_label_clean, order), y = fold_change))+
  geom_line(aes(group = protein_label, colour = species), size = 2, alpha = 0.6)+
  geom_point(aes(colour = sample_cat), size = 5, alpha = 0.9)+
  scale_y_log10()+
  scale_color_manual(values = c("#9EB3BE", "#DB4F48","#E8B840" ))+
  theme_linedraw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y = "fold change relative to \nweek one cut-off",
       x  = "")


fold_change_plot_subset <- fold_change %>%
  filter(protein_label%in%proteins_paper) %>% 
  mutate(species_facet = case_when(
    str_detect(protein_label_clean, "Pk") ~ "P. knowlesi",
    str_detect(protein_label_clean, "Pv") ~ "P. vivax"
  )) %>% 
  mutate(sample_cat = paste0("P. ", sample_cat)) %>% 
  mutate(species_facet = factor(species_facet, levels = c("P. vivax", "P. knowlesi"))) %>% 
  ggplot(aes(x = fct_reorder(protein_label_clean, order), y = fold_change))+
  geom_hline(yintercept = 1, linetype = 2, colour = "black")+
  geom_line(aes(group = protein_label, colour = species_facet), size = 1.5, alpha = 0.6)+
  geom_point(aes(colour = sample_cat), size = 4, alpha = 0.8)+
  scale_y_log10(limits = c(0.1, 170))+
  facet_grid(~species_facet, scales = "free_x", space = "free", labeller = labeller(species_facet = c("P. vivax" = "P. vivax proteins", "P. knowlesi" = "P. knowlesi proteins") ) )+
  scale_color_manual(values = c("#E8B840" ,"#DB4F48"),
                     guide = guide_legend(override.aes = list(linetype = c(0,0))))+
  theme_linedraw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 17),
        strip.background = element_blank(),
        
        strip.text = element_text(colour = "black", face = "italic",hjust = 0, size  = 18),
        legend.text = element_text(face = "italic", size = 18),
        legend.position = "bottom",
        panel.grid = element_line(colour = "black"))+
  labs(y = "fold-change one week\npost-infection compared to\nsero-positivity cut-off",
       x  = "recombinant protein",
       colour = "species cohort")


##########################################################
############ FOLD CHANGE and PROTEIN SIMILARITY PLOTS
###########################################################

identity <- readxl::read_excel(here::here("protein_info/SupplementaryFile1-ProteinInfo.xlsx"), sheet = "identity") %>% 
  dplyr::rename(protein_label = 1,
                identity = 2) %>% 
  mutate(identity = as.numeric(identity))

fold_change_orthologue <-fold_change %>% 
  filter(protein_label%in%proteins_paper) %>% 
  left_join(identity) %>% 
  mutate(species_facet = case_when(
    str_detect(protein_label_clean, "Pk") ~ "P. knowlesi",
    str_detect(protein_label_clean, "Pv") ~ "P. vivax"
  )) %>% 
  mutate(species_facet = factor(species_facet, levels = c("P. vivax", "P. knowlesi"))) %>% 
  mutate(sample_cat = paste0("P. ", sample_cat)) %>% 
  mutate(sample_cat = factor(sample_cat, levels = c("P. vivax", "P. knowlesi"))) %>% 
  mutate(species_facet = factor(species_facet, levels = c("P. vivax", "P. knowlesi"))) %>% 
  filter(!species_facet==sample_cat) %>% 
  ggplot(aes(y = fold_change, x = identity, colour = species_facet))+
  geom_point(size = 3)+
  scale_color_manual(values = c("#DB4F48", "#E8B840" ),
                     guide = guide_legend(override.aes = list(linetype = c(0,0))))+
  geom_hline(yintercept = 1, linetype = 2, colour = "grey")+
  ggnewscale::new_scale_color()+
  ggrepel::geom_label_repel(aes(label = protein_label_clean, colour = species_facet), size = 4.3,force = 3,force_pull = 3,label.padding = 0.1,alpha = 0.9,)+
  scale_color_manual(values = c("#DB4F48", "#BF6B04" ),
                     guide = guide_legend(override.aes = list(linetype = c(0,0))))+
  scale_y_log10(limits = c(0.1, 170))+
  ggpubr::stat_cor(p.accuracy = 0.01, r.accuracy = 0.01,method = "spearman", colour = "black",label.y = 2.1, inherit.aes = TRUE, size = 4.5,cor.coef.name = "rho")+
  facet_wrap(~sample_cat, labeller = labeller(sample_cat = c("P. vivax" = "P. vivax cohort", "P. knowlesi" = "P. knowlesi cohort")))+
  theme_linedraw()+
  labs(y = "fold change compared to\nsero-positivity cut-off",
       x  = "recombinant fragment protein identity (%) \nwith Pv and Pk orthologues ")+
  theme(axis.text.x = element_text(angle = 0, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 17),
        panel.grid = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", face = "italic",hjust = 0, size  = 18),
        legend.position = "none")

fold_change_orthologue

### same plot but with both species in both panels to help merge with previous figure...add as supplementary
includes_both_spec = fold_change %>% 
  filter(protein_label%in%proteins_paper) %>% 
  left_join(identity) %>% 
  mutate(species_facet = case_when(
    str_detect(protein_label_clean, "Pk") ~ "P. knowlesi",
    str_detect(protein_label_clean, "Pv") ~ "P. vivax"
  )) %>% 
  mutate(species_facet = factor(species_facet, levels = c("P. vivax", "P. knowlesi"))) %>% 
  mutate(sample_cat = paste0("P. ", sample_cat)) %>% 
  mutate(sample_cat = factor(sample_cat, levels = c("P. vivax", "P. knowlesi"))) %>% 
  mutate(species_facet = factor(species_facet, levels = c("P. vivax", "P. knowlesi"))) %>% 
  # filter(!species_facet==sample_cat) %>% 
  ggplot(aes(y = fold_change, x = identity, colour = species_facet))+
  geom_point(size = 3)+
  ggrepel::geom_label_repel(aes(label = protein_label_clean), size = 4,force = 3,force_pull = 3,label.padding = 0.1,alpha = 0.8)+
  scale_color_manual(values = c("#DB4F48", "#E8B840" ),
                     guide = guide_legend(override.aes = list(linetype = c(0,0))))+
  scale_y_log10(limits = c(0.1, 170))+
  ggpubr::stat_cor(p.accuracy = 0.01, r.accuracy = 0.01,method = "spearman", colour = "black")+
  facet_wrap(~sample_cat, labeller = labeller(sample_cat = c("P. vivax" = "P. vivax cohort", "P. knowlesi" = "P. knowlesi cohort")))+
  theme_linedraw()+
  labs(y = "fold change relative to \nweek one cut-off",
       x  = "recombinant fragment protein identity (%) \nwith Pv and Pk orthologues ")+
  theme(axis.text.x = element_text(angle = 0, size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 17),
        panel.grid = element_line(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(colour = "darkgrey", face = "italic",hjust = 0, size  = 18),
        legend.position = "none")


for_correlation <- fold_change %>% 
  filter(protein_label%in%proteins_paper) %>% 
  left_join(identity) %>% 
  mutate(species_facet = case_when(
    str_detect(protein_label_clean, "Pk") ~ "P. knowlesi",
    str_detect(protein_label_clean, "Pv") ~ "P. vivax"
  )) %>% 
  mutate(species_facet = factor(species_facet, levels = c("P. vivax", "P. knowlesi"))) %>% 
  mutate(sample_cat = paste0("P. ", sample_cat)) %>% 
  mutate(sample_cat = factor(sample_cat, levels = c("P. vivax", "P. knowlesi"))) %>% 
  mutate(species_facet = factor(species_facet, levels = c("P. vivax", "P. knowlesi"))) %>% 
  filter(!species_facet==sample_cat) %>% 
  group_by(species)
###########################################################################################


design <- "
AAABBB
"

figure2 = (fold_change_plot_subset + fold_change_orthologue)  +plot_annotation(tag_levels = "A")+ plot_layout(design = design) &
  theme(plot.tag = element_text(size = 22)) 

svg(here::here("figures/fig2_crossreactivity.svg"), width = 12.8, height = 5.5)
print(figure2)
dev.off()



###########################################################################################

#### double check correlations calculated manually
correlation_fold <- for_correlation %>% 
  drop_na(identity) %>% 
  dplyr::select(protein_label_clean, species, fold_change, identity) %>% 
  group_by(species) %>% 
  nest() %>% 
  mutate(
    test = map(data, ~ cor.test(.x$fold_change, .x$identity, method = "spearman")), 
    tidied = map(test, tidy)
  ) %>% 
  unnest(cols = tidied) %>% 
  select(-data, -test)
######----------Plot sample to day

time_individuals <- cats_for_plot %>% 
  dplyr::select(id, day, day_cat, sample_cat) %>% 
  distinct() %>% 
  dplyr::group_by(id) %>% add_count() %>% filter(n>1) %>% 
  dplyr::mutate(max_day = max(day)) %>% 
  ungroup() %>% 
  
  ggplot(aes(x = day, y = fct_reorder(id,max_day), colour = sample_cat, shape =day_cat))+
  geom_line(aes(group = id))+
  geom_point(size = 1.5, alpha = 0.9)+
  theme_void()+
  theme(axis.text.x = element_text(colour = "black"),
        axis.title = element_text(colour = "black"),
        axis.title.y = element_text(angle = 90),
        panel.grid.major.x = element_line(colour = "grey"),
        panel.grid.minor.x = element_line(colour = "lightgrey"),
        legend.position = "bottom")+
  scale_shape_manual(values = c(16,2, 4, 13, 8,9,15))+
  labs(x = "day", y = "individual",
       colour = "Infection species",
       shape = "time grouping")+
  scale_colour_manual(values = c('#F24141',  "#F2B705"),
                      labels = c(expression(italic("P. knowlesi")), expression(italic("P. vivax"))))

svg(here::here("figures/figs14_time.svg"), width = 10, height = 12)
print(time_individuals)
dev.off()

cats_for_plot %>% 
  select(id, day, day_cat, sample_cat) %>% 
  distinct() %>% 
  group_by(id) %>% add_count() %>% filter(n>1) %>% ungroup() %>% group_by(id, sample_cat)

cats_for_plot %>% filter(day!=0) %>% 
  select(id, day, day_cat, sample_cat) %>% 
  distinct() %>% count( sample_cat)

cats_for_plot %>% 
  filter(day!=0) %>% 
  select(id, day, day_cat, sample_cat) %>% 
  distinct() %>% select(id, sample_cat) %>% distinct() %>% count(sample_cat)

cats_for_plot %>% 
  select(id, day, day_cat, sample_cat) %>% 
  distinct() %>% select(id, sample_cat) %>% distinct() %>% count(sample_cat)

