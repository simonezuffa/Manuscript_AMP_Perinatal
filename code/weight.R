setwd("~/Desktop/Manuscript_AMP_Perinatal")

library(tidyverse)
library(readxl)
library(ggpubr)
library(effectsize)
library(emmeans)

metadata_microbiome <- read_csv("data_weight/meta_microbiome.csv") %>%
  dplyr::filter(!str_detect(pattern = "Zymo|Blank", SampleID))
metadata_metabolomics <- read_csv("data_weight/meta_metabolomics.csv") %>%
  dplyr::filter(Batch == 1)

metadata_integration <- read_csv("data_weight/meta_integration.csv")

unique_microbiome <- metadata_microbiome %>% distinct(Cage, Hypothesis, Batch, .keep_all = TRUE)
unique_metabolomics <- metadata_metabolomics %>% distinct(Cage, Hypothesis, Batch, .keep_all = TRUE)

################################
# Calculate weight differences #
################################
hypo1_pup_weight <- read_excel("data_weight/weight.xlsx", sheet = 2)
hypo2_pup_weight <- read_excel("data_weight/weight.xlsx", sheet = 3)

# Antepartum
weight_plot_hypo1 <- hypo1_pup_weight %>%
  ggboxplot(x = "AMP", y = "Weight", add = "jitter", 
            add.params = list(color = "AMP", alpha = 0.8, size = 1),
            palette = c("#287DAB", "#E5BF86"), title = "Pups weight at weaning",
            xlab = FALSE, ylab = "Weight (g)", legend = "none") +  ylim(0, 25) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

model1 <- hypo1_pup_weight %>%
  lm(formula = Weight ~ AMP + Sex)
summary(model1)

# Postpartum
weight_plot_hypo2 <- hypo2_pup_weight %>%
  ggboxplot(x = "Treatment", y = "Weight", add = "jitter", 
            add.params = list(color = "Treatment", alpha = 0.8, size = 1),
            palette = c("#287DAB", "#E5BF86"), title = "Pups weight at weaning",
            xlab = FALSE, ylab = "Weight (g)", legend = "none") +  ylim(0, 25) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

model2 <- hypo2_pup_weight %>% 
  lm(formula = Weight ~ Treatment + Sex)
summary(model2)

weight_plot <- ggarrange(weight_plot_hypo1, weight_plot_hypo2, nrow = 1)

#ggsave(plot = weight_plot, filename = "weight_pups.svg", device = "svg", dpi = "retina", width = 3, height = 1.85)


# Combined data as suggested by Reviewer and use estimated marginal means (emmeans)
hypo1_final <- hypo1_pup_weight %>%
  dplyr::select(1,5,3,4)
colnames(hypo1_final) <- c("hypo", "group", "sex", "weight")

hypo2_final <- hypo2_pup_weight %>%
  dplyr::select(1,5,3,4)
colnames(hypo2_final) <- c("hypo", "group", "sex", "weight")

hypo_combined <- rbind(hypo1_final, hypo2_final) %>%
  dplyr::mutate(group = gsub("TRUE", "AMP", group)) %>%
  dplyr::mutate(group = gsub("FALSE", "PBS", group))

model3 <- hypo_combined %>% 
  lm(formula = weight ~ group * hypo + sex)
summary(model3)

em_results <- emmeans(model3, ~ group | hypo)
group_comparisons <- pairs(em_results)
summary(group_comparisons)
pairs(group_comparisons, by = NULL)
