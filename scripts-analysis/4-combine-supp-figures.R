# a script to combine the supplementary figure files which many are from both recent PK and either Plasmodium

library(patchwork)

rds_figures <- list.files(here::here("figures/rds"),full.names = T)

figures = mapply(read_rds, rds_figures)


figure_s6a <- figures$`/Users/cb/Library/CloudStorage/OneDrive-CharlesDarwinUniversity/projects/1-Serology/ClassificationProject/PkSeroClassification/figures/rds/cross_validation_pred-recentpk.rds`
figure_s6b <- figures$`/Users/cb/Library/CloudStorage/OneDrive-CharlesDarwinUniversity/projects/1-Serology/ClassificationProject/PkSeroClassification/figures/rds/cross_validation_pred-plasmodium.rds`


figure_S6 <- figure_s6a/figure_s6b + plot_layout(guides = "collect") + plot_annotation(tag_level = "A")


svg(here::here("figures/figure_S6.svg"),width = 12, height = 9)
print(figure_S6)
dev.off()

figure_s7a <- figures$`/Users/cb/Library/CloudStorage/OneDrive-CharlesDarwinUniversity/projects/1-Serology/ClassificationProject/PkSeroClassification/figures/rds/acute_preds-recentpk.rds`
figure_s7b <- figures$`/Users/cb/Library/CloudStorage/OneDrive-CharlesDarwinUniversity/projects/1-Serology/ClassificationProject/PkSeroClassification/figures/rds/acute_preds_plasmodium.rds`

figure_S7 <- figure_s7a+figure_s7b + plot_layout(guides = "collect") + plot_annotation(tag_level = "A")

svg(here::here("figures/figure_S7.svg"),width = 8, height = 5)
print(figure_S7)
dev.off()

figure_s8a <- figures$`/Users/cb/Library/CloudStorage/OneDrive-CharlesDarwinUniversity/projects/1-Serology/ClassificationProject/PkSeroClassification/figures/rds/holdout_plot-recentpk.rds`
figure_s8b <- figures$`/Users/cb/Library/CloudStorage/OneDrive-CharlesDarwinUniversity/projects/1-Serology/ClassificationProject/PkSeroClassification/figures/rds/holdout_plot_plasmodium.rds`

figure_S8 <- figure_s8a/figure_s8b + plot_layout(guides = "collect") + plot_annotation(tag_level = "A")

svg(here::here("figures/figure_S8.svg"),width = 12, height = 9)
print(figure_S8)
dev.off()

figure_s9a <- figures$`/Users/cb/Library/CloudStorage/OneDrive-CharlesDarwinUniversity/projects/1-Serology/ClassificationProject/PkSeroClassification/figures/rds/correlation_pred_new-recentpk.rds`
figure_s9b <- figures$`/Users/cb/Library/CloudStorage/OneDrive-CharlesDarwinUniversity/projects/1-Serology/ClassificationProject/PkSeroClassification/figures/rds/correlation_pred_new_plasmodium.rds`

figure_S9 <- figure_s9a/figure_s9b + plot_layout(guides = "collect") + plot_annotation(tag_level = "A")


svg(here::here("figures/figure_S9.svg"),width = 12, height = 10)
print(figure_S9)
dev.off()

figure_s10a <- figures$`/Users/cb/Library/CloudStorage/OneDrive-CharlesDarwinUniversity/projects/1-Serology/ClassificationProject/PkSeroClassification/figures/rds/training_reduced_ind_pk.rds`
figure_s10b <- figures$`/Users/cb/Library/CloudStorage/OneDrive-CharlesDarwinUniversity/projects/1-Serology/ClassificationProject/PkSeroClassification/figures/rds/training_reduced_ind_plasmodium.rds`

figure_S10 <- figure_s10a/figure_s10b + plot_layout(guides = "collect") + plot_annotation(tag_level = "A")


svg(here::here("figures/figure_S10.svg"),width = 12, height = 10)
print(figure_S10)
dev.off()

