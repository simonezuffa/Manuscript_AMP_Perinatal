setwd("~/Desktop/Manuscript_AMP_Perinatal")

library(tidyverse)
library(vegan)
library(caret)

# Import
microbiome_ante <- read_csv("data_integration/ante_microbiome.csv")
microbiome_post <- read_csv("data_integration/post_microbiome.csv")
metadata_microbiome <- read_csv("data_integration/metadata_microbiome.csv") %>% 
  dplyr::filter(SampleID %in% c(microbiome_ante$SampleID, microbiome_post$SampleID)) %>%
  dplyr::mutate(animalID = gsub("O_", "E_", SampleID)) %>%
  dplyr::mutate(animalID = gsub("_Hypo1_1", "_1", animalID)) %>%
  dplyr::mutate(animalID = gsub("_Hypo2_1", "_2", animalID)) %>%
  dplyr::mutate(animalID = gsub("ME_9", "ME_C9", animalID))

metabolome_ante <- read_csv("data_integration/ante_metabolome.csv")
metabolome_post <- read_csv("data_integration/post_metabolome.csv")
metadata_metabolome <- read_csv("data_integration/metadata_metabolome.csv") %>% 
  dplyr::filter(Sample %in% c(metabolome_ante$SampleID, metabolome_post$SampleID)) %>%
  dplyr::mutate(animalID = gsub("-", "_", Sample)) %>%
  dplyr::mutate(animalID = paste(animalID, Hypothesis, sep = "_"))

metadata_combined <- metadata_microbiome %>% 
  full_join(metadata_metabolome, by = "animalID") %>%
  dplyr::filter(!is.na(Mother_Infant.x)) %>%
  dplyr::filter(!is.na(Mother_Infant.y)) %>%
  dplyr::filter(animalID != "ME_C4_B21_2")

all(metadata_combined$Cage.x == metadata_combined$Cage.y)
all(metadata_combined$Hypothesis.x == metadata_combined$Hypothesis.y)
all(metadata_combined$Batch.x == metadata_combined$Batch.y)

all(metadata_combined$Amp == metadata_combined$AMP) # ME_C9_2 is Amp TRUE

metadata_final <- metadata_combined %>% 
  dplyr::select(animalID, SampleID, Sample, AMP, Days_Post_Birth, Cage.x, Hypothesis.x)

# Fix sample names
microbiome_ante_final <- microbiome_ante %>% 
  dplyr::filter(SampleID %in% metadata_final$SampleID) %>%
  left_join(metadata_final %>% dplyr::select(SampleID, animalID)) %>%
  dplyr::select(-SampleID) %>% arrange(animalID) %>% 
  column_to_rownames("animalID") %>% 
  decostand(method = "rclr") %>% 
  select_at(vars(-one_of(caret::nearZeroVar(., names = TRUE))))

microbiome_post_final <- microbiome_post %>% 
  dplyr::filter(SampleID %in% metadata_final$SampleID) %>%
  left_join(metadata_final %>% dplyr::select(SampleID, animalID)) %>%
  dplyr::select(-SampleID) %>% arrange(animalID) %>% 
  column_to_rownames("animalID") %>% 
  decostand(method = "rclr") %>% 
  select_at(vars(-one_of(caret::nearZeroVar(., names = TRUE))))

metabolome_ante_final <- metabolome_ante %>% 
  dplyr::filter(SampleID %in% metadata_final$Sample) %>%
  left_join(metadata_final %>% dplyr::select(Sample, animalID), by = c("SampleID" = "Sample")) %>%
  dplyr::select(-SampleID) %>% arrange(animalID) %>% 
  column_to_rownames("animalID") %>% 
  decostand(method = "rclr") %>% 
  select_at(vars(-one_of(caret::nearZeroVar(., names = TRUE))))

metabolome_post_final <- metabolome_post %>% 
  dplyr::filter(SampleID %in% metadata_final$Sample) %>%
  left_join(metadata_final %>% dplyr::select(Sample, animalID), by = c("SampleID" = "Sample")) %>%
  dplyr::select(-SampleID) %>% arrange(animalID) %>% 
  column_to_rownames("animalID") %>% 
  decostand(method = "rclr") %>% 
  select_at(vars(-one_of(caret::nearZeroVar(., names = TRUE))))

metabolomics_overlap <- read_csv("data_integration/metabolomics_overlap.csv") %>%
  dplyr::mutate(dataset.x = gsub("MSV000089558_", "", dataset.x)) %>%
  dplyr::mutate(dataset.y = gsub("MSV000092652_", "", dataset.y)) %>%
  dplyr::filter(!(dataset.x %in% c(6366, 20524, 23473, 16188))) %>% 
  dplyr::filter(str_detect(pattern = "bile acid|Glucopyranosyl|QUINOLINECARBOXYLIC|BETAINE|
                           CARNITINE|carnitine|putrescine|Spermidine|Ketodeoxycholic", Compound_Name) | 
                  str_detect(pattern ="Carnitine", syn_lib))

# I have 50 animals from the Antepartum cohort and 53 from the Postpartum cohort 
# with matching microbiome and metabolomics samples

meta_ante <- metadata_final %>% dplyr::filter(Hypothesis.x == 1) %>% arrange(animalID)
meta_post <- metadata_final %>% dplyr::filter(Hypothesis.x == 2) %>% arrange(animalID)

all(rownames(microbiome_ante_final) == rownames(metabolome_ante_final))
all(rownames(microbiome_ante_final) == meta_ante$animalID)

all(rownames(microbiome_post_final) == rownames(metabolome_post_final))
all(rownames(microbiome_post_final) == meta_post$animalID)

metabolome_ante_overlap <- metabolome_ante_final %>% 
  dplyr::select(metabolomics_overlap$dataset.x)

metabolome_post_overlap <- metabolome_post_final %>% 
  dplyr::select(metabolomics_overlap$dataset.y)

# PLS between the two omics datasets
pls_ante <- mixOmics::pls(metabolome_ante_final, microbiome_ante_final, ncomp = 2, scale = TRUE)
pls_ante_plot <- mixOmics::plotIndiv(pls_ante, comp = 1:2, rep.space= 'XY-variate', 
                                     group = meta_ante$AMP, ind.names = FALSE, legend = TRUE, 
                                     title = "PLS comp 1 - 2, XY-space", pch = 20, 
                                     style = "lattice", centroid = TRUE, ellipse = TRUE)
cor(pls_ante$variates$X, pls_ante$variates$Y) %>% diag()

pls_post <- mixOmics::pls(metabolome_post_final, microbiome_post_final, ncomp = 2, scale = TRUE)
pls_post_plot <- mixOmics::plotIndiv(pls_post, comp = 1:2, rep.space= 'XY-variate', 
                                     group = meta_post$AMP, ind.names = FALSE, legend = TRUE, 
                                     title = "PLS comp 1 - 2, XY-space", pch = 20, 
                                     style = "lattice", centroid = TRUE, ellipse = TRUE)
# there is an outlier (ME_C4_B21_2) --> remove from the beginning 
cor(pls_post$variates$X, pls_post$variates$Y) %>% diag()


pls_ante_overlap <- mixOmics::pls(metabolome_ante_overlap, microbiome_ante_final, ncomp = 2, scale = TRUE)
pls_ante_overlap_plot <- mixOmics::plotIndiv(pls_ante, comp = 1:2, rep.space= 'XY-variate', 
                                     group = meta_ante$AMP, ind.names = FALSE, legend = TRUE, 
                                     title = "PLS comp 1 - 2, XY-space", pch = 20, 
                                     style = "lattice", centroid = TRUE, ellipse = TRUE)
cor(pls_ante_overlap$variates$X, pls_ante_overlap$variates$Y) %>% diag()

pls_post_overlap <- mixOmics::pls(metabolome_post_overlap, microbiome_post_final, ncomp = 2, scale = TRUE)
pls_post_overlap_plot <- mixOmics::plotIndiv(pls_post, comp = 1:2, rep.space= 'XY-variate', 
                                             group = meta_post$AMP, ind.names = FALSE, legend = TRUE, 
                                             title = "PLS comp 1 - 2, XY-space", pch = 20, 
                                             style = "lattice", centroid = TRUE, ellipse = TRUE)
cor(pls_post_overlap$variates$X, pls_post_overlap$variates$Y) %>% diag()


# Prepare data for integration
data_ante <- list(Metabo = as.matrix(metabolome_ante_final), 
                  Microb = as.matrix(microbiome_ante_final))
response_ante <- factor(meta_ante$AMP)
design_ante <- matrix(0.9, ncol = length(data_ante), nrow = length(data_ante), 
                      dimnames = list(names(data_ante), names(data_ante)))
diag(design_ante) = 0

data_post <- list(Metabo = as.matrix(metabolome_post_final), 
                  Microb = as.matrix(microbiome_post_final))
response_post <- factor(meta_post$AMP)
design_post <- matrix(0.9, ncol = length(data_post), nrow = length(data_post), 
                      dimnames = list(names(data_post), names(data_post)))
diag(design_post) = 0

data_ante_ov <- list(Metabo = as.matrix(metabolome_ante_overlap), 
                  Microb = as.matrix(microbiome_ante_final))
response_ante <- factor(meta_ante$AMP)
design_ante_ov <- matrix(0.9, ncol = length(data_ante_ov), nrow = length(data_ante_ov), 
                      dimnames = list(names(data_ante_ov), names(data_ante_ov)))
diag(design_ante_ov) = 0

data_post_ov <- list(Metabo = as.matrix(metabolome_post_overlap), 
                     Microb = as.matrix(microbiome_post_final))
response_post <- factor(meta_post$AMP)
design_post_ov <- matrix(0.9, ncol = length(data_post_ov), nrow = length(data_post_ov), 
                         dimnames = list(names(data_post_ov), names(data_post_ov)))
diag(design_post_ov) = 0


# DIABLO 
list.keepX = list(Metabo = c(20, 20), Microb = c(20, 20)) # keep top 20

DIABLO_ante <- mixOmics::block.splsda(X = data_ante, Y = response_ante, 
                                      ncomp = 2, design = design_ante,
                                      keepX = list.keepX)
perf_DIABLO_ante <- mixOmics::perf(DIABLO_ante, validation = 'loo')
plot(perf_DIABLO_ante, legend = FALSE)

DIABLO_post <- mixOmics::block.splsda(X = data_post, Y = response_post, 
                                      ncomp = 2, design = design_post,
                                      keepX = list.keepX)
perf_DIABLO_post <- mixOmics::perf(DIABLO_post, validation = 'loo')
plot(perf_DIABLO_post, legend = FALSE)

DIABLO_ante_ov <- mixOmics::block.splsda(X = data_ante_ov, Y = response_ante, 
                                      ncomp = 2, design = design_ante_ov,
                                      keepX = list.keepX)
perf_DIABLO_ante_ov <- mixOmics::perf(DIABLO_ante_ov, validation = 'loo')
plot(perf_DIABLO_ante_ov, legend = FALSE)

DIABLO_post_ov <- mixOmics::block.splsda(X = data_post_ov, Y = response_post, 
                                         ncomp = 2, design = design_post_ov,
                                         keepX = list.keepX)
perf_DIABLO_post_ov <- mixOmics::perf(DIABLO_post_ov, validation = 'loo')
plot(perf_DIABLO_post_ov, legend = FALSE)

# Score plots 
mixOmics::plotIndiv(DIABLO_ante, legend = TRUE,
                    title = 'DIABLO', col.per.group = c("#287DAB", "#E5BF86"),
                    style = "lattice", ellipse = TRUE, centroid = TRUE, pch = 20, 
                    X.label = "Component 1", Y.label = "Component 2")

mixOmics::plotIndiv(DIABLO_ante_ov, legend = TRUE,
                    title = 'DIABLO', col.per.group = c("#287DAB", "#E5BF86"),
                    style = "lattice", ellipse = TRUE, centroid = TRUE, pch = 20, 
                    X.label = "Component 1", Y.label = "Component 2")

mixOmics::plotIndiv(DIABLO_post, legend = TRUE,
                    title = 'DIABLO', col.per.group = c("#287DAB", "#E5BF86"),
                    style = "lattice", ellipse = TRUE, centroid = TRUE, pch = 20, 
                    X.label = "Component 1", Y.label = "Component 2")

mixOmics::plotIndiv(DIABLO_post_ov, legend = TRUE,
                    title = 'DIABLO', col.per.group = c("#287DAB", "#E5BF86"),
                    style = "lattice", ellipse = TRUE, centroid = TRUE, pch = 20, 
                    X.label = "Component 1", Y.label = "Component 2")

# Circos plots 
#pdf("circosPlot_DIABLO_ante_overlap.pdf", width = 4, height = 4)
mixOmics::circosPlot(DIABLO_ante_ov, cutoff = 0.7, line = FALSE, size.labels = 0.5, comp = 1,
                     color.blocks = c("#287DAB", "#E5BF86"), showIntraLinks = FALSE,
                     size.variables = 0.5, size.legend = 0.5)
#dev.off()

#pdf("circosPlot_DIABLO_post_overlap.pdf", width = 10, height = 10)
mixOmics::circosPlot(DIABLO_post_ov, cutoff = 0.8, line = FALSE, size.labels = 0.5, comp = 1,
                     color.blocks = c("#287DAB", "#E5BF86"), showIntraLinks = FALSE,
                     size.variables = 0.5, size.legend = 0.5)
#dev.off()
