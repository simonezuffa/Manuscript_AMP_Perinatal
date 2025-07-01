setwd("~/Desktop/Manuscript_AMP_Perinatal")

library(tidyverse)
library(ggpubr)
library(vegan)
library(caret)
library(patchwork)
library(UpSetR)
library(Spectra)
library(MsBackendMgf)

##########
# Import #
##########
data_weight <- read_csv("data_metabolomics/weight.csv")
metadata <- read_tsv("data_metabolomics/metadata_metabolomics.tsv")
metadata_filter <- metadata %>%
  dplyr::filter(Mother_Infant == "Mother") %>%
  dplyr::filter(Hypothesis != 3) %>%
  dplyr::filter(Type == "Fecal") %>%
  dplyr::filter(Days_Post_Birth != 28)

rev_cosine <- read_delim("data_metabolomics/rev_cosine.tsv") %>% 
  dplyr::select(1,2) %>% dplyr::rename(rev_cosine = COMPOUND_NAME)

syn_lib <- read_csv("data_metabolomics/syn_lib.csv") %>%
  dplyr::select(1, 48) %>% group_by(SpectrumID) %>%
  summarise(ba_lib = paste(unique(synlib_compound_name), collapse = ";"), .groups = "drop")

ba_massql_pre <- read_tsv("data_metabolomics/ba_filter_antepartum.tsv")
ba_massql_post <- read_tsv("data_metabolomics/ba_filter_postpartum.tsv")

# Metadata contains info on all acquired samples. Initial study design contained also
# hypo3, looking at direct AMP injection in infants. Due to severe batch problems and
# lack of samples for infants, analysis was performed only on maternal fecal samples
# from hypo1 and hypo2. Cage identify the animals, as they were single housed.

# Feature tables acquired at different timepoints so they will be analyze separately
feature_MSV000089558 <- read_csv("data_metabolomics/gnps_quant_MSV000089558.csv")
feature_MSV000092652 <- read_csv("data_metabolomics/gnps_quant_MSV000092652.csv")

# Annotations from FBMN
annotations_MSV000089558 <- read.delim("data_metabolomics/annotations_MSV000089558.tsv") %>%
  dplyr::filter(!str_detect(pattern = "REFRAME", LibraryName)) # remove drug library
annotations_MSV000089558$X.Scan. <- as.character(annotations_MSV000089558$X.Scan.)

annotations_MSV000092652 <- read.delim("data_metabolomics/annotations_MSV000092652.tsv") %>%
  dplyr::filter(!str_detect(pattern = "REFRAME", LibraryName)) # remove drug library
annotations_MSV000092652$X.Scan. <- as.character(annotations_MSV000092652$X.Scan.)

# Annotations from synthesis library
annotations_MSV000089558_syn <- read.delim("data_metabolomics/annotations_MSV000089558_syn.tsv") %>%
  dplyr::filter(str_detect(pattern = "Synthesis_library_filter.mgf|Carnitines_library_beta.mgf", LibraryName)) %>%
  dplyr::select(2, 15)
annotations_MSV000089558_syn$X.Scan. <- as.character(annotations_MSV000089558_syn$X.Scan.)

annotations_MSV000092652_syn <- read.delim("data_metabolomics/annotations_MSV000092652_syn.tsv") %>%
  dplyr::filter(str_detect(pattern = "Synthesis_library_filter.mgf|Carnitines_library_beta.mgf", LibraryName)) %>%
  dplyr::select(2, 15)
annotations_MSV000092652_syn$X.Scan. <- as.character(annotations_MSV000092652_syn$X.Scan.)

info_feature_MSV000089558 <- feature_MSV000089558 %>% dplyr::select(1:3,7)
colnames(info_feature_MSV000089558) <- c("Feature", "mz", "RT", "Corr_ID")
info_feature_MSV000089558$Feature <- as.character(info_feature_MSV000089558$Feature)
info_feature_MSV000089558_complete <- info_feature_MSV000089558 %>% 
  left_join(annotations_MSV000089558, by = c("Feature" = "X.Scan.")) %>% 
  dplyr::select(1:5,18,24) %>%
  left_join(rev_cosine, by = c("SpectrumID" = "query_id")) %>%
  left_join(syn_lib, by = c("SpectrumID" = "SpectrumID")) %>%
  left_join(annotations_MSV000089558_syn, by = c("Feature" = "X.Scan.")) %>%
  dplyr::rename(Compound_Name = `Compound_Name.x`) %>%
  dplyr::rename(syn_lib = `Compound_Name.y`)

info_feature_MSV000092652 <- feature_MSV000092652 %>% dplyr::select(1:3,7)
colnames(info_feature_MSV000092652) <- c("Feature", "mz", "RT", "Corr_ID")
info_feature_MSV000092652$Feature <- as.character(info_feature_MSV000092652$Feature)
info_feature_MSV000092652_complete <- info_feature_MSV000092652 %>% 
  left_join(annotations_MSV000092652, by = c("Feature" = "X.Scan.")) %>% 
  dplyr::select(1:5,18,24) %>%
  left_join(rev_cosine, by = c("SpectrumID" = "query_id")) %>%
  left_join(syn_lib, by = c("SpectrumID" = "SpectrumID")) %>%
  left_join(annotations_MSV000092652_syn, by = c("Feature" = "X.Scan.")) %>%
  dplyr::rename(Compound_Name = `Compound_Name.x`) %>%
  dplyr::rename(syn_lib = `Compound_Name.y`)


# Transpose tables and keep only maternal samples of interest
data_MSV000089558 <- feature_MSV000089558 %>%
  column_to_rownames("row ID") %>% dplyr::select(contains("Peak")) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("SampleID") %>% 
  arrange(SampleID) %>% distinct(SampleID, .keep_all = TRUE)
data_MSV000089558$SampleID <- gsub(".mzML Peak area", "", data_MSV000089558$SampleID)

data_MSV000089558_filter <- data_MSV000089558 %>% 
  dplyr::filter(SampleID %in% metadata_filter$Sample | 
                  str_detect(pattern = regex("QC|Blank", ignore_case = TRUE), SampleID))

data_MSV000092652 <- feature_MSV000092652 %>%
  column_to_rownames("row ID") %>% dplyr::select(contains("Peak")) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("SampleID") %>% 
  arrange(SampleID) %>% distinct(SampleID, .keep_all = TRUE)
data_MSV000092652$SampleID <- gsub(".mzML Peak area", "", data_MSV000092652$SampleID)

data_MSV000092652_filter <- data_MSV000092652 %>% 
  dplyr::filter(SampleID %in% metadata_filter$Sample | 
                  str_detect(pattern = regex("QC|Blank", ignore_case = TRUE), SampleID))


##################
# PCA - Raw Data #
##################

# MSV000089558
data_sample_clr <- decostand(data_MSV000089558_filter %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sample_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("Sample") %>% left_join(metadata_filter)

PCA_MSV000089558_plots <- list()

for (i in c("Hypothesis", "Batch", "Days_Post_Birth", "Type", "AMP", "Cage", "Run")) {
  
  PCA_plot <- PCA_whole_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_MSV000089558_plots[[i]] <- PCA_plot
  
}

PCA_MSV000089558_plots_final <- wrap_plots(PCA_MSV000089558_plots, nrow = 3)
# There are some samples from hypo2 --> remove them (ME_C10_B1, ME_C5_B14, ME_C5_B21, ME_C5_B7)
# Good clustering of QCs


# MSV000092652
data_sample_clr <- decostand(data_MSV000092652_filter %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sample_clr %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("Sample") %>% left_join(metadata_filter)

PCA_MSV000092652_plots <- list()

for (i in c("Hypothesis", "Batch", "Days_Post_Birth", "Type", "AMP", "Cage", "Run")) {
  
  PCA_plot <- PCA_whole_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_MSV000092652_plots[[i]] <- PCA_plot
  
}

PCA_MSV000092652_plots_final <- wrap_plots(PCA_MSV000092652_plots, nrow = 3)
# Good clustering of QCs


# Fecal samples
data_MSV000089558_feces <- data_MSV000089558_filter %>%
  dplyr::filter(!(SampleID %in% c("ME_C10_B1", "ME_C5_B14", "ME_C5_B21", "ME_C5_B7"))) %>% # removed because from hypo2
  dplyr::filter(!(SampleID %in% c("ME_C8-B1", "ME_C10-B14", "ME_C8-B14"))) %>% # removed because IS was missing or twice as much (problem sample preparation)
  dplyr::filter(!(SampleID %in% c("Blank_01_no_inj", "Blank_02_no_inj")))

data_MSV000092652_feces <- data_MSV000092652_filter %>%
  dplyr::filter(SampleID != "Blank_001_no_inj")


# Investigate total peak area
data_TIC_MSV000089558 <- data.frame(TIC = rowSums(data_MSV000089558_feces %>% 
                                                    column_to_rownames("SampleID"))) %>% rownames_to_column("SampleID")
data_TIC_MSV000092652 <- data.frame(TIC = rowSums(data_MSV000092652_feces %>% 
                                                    column_to_rownames("SampleID"))) %>% rownames_to_column("SampleID")

# Check TIC Sample, QCmix and Blank

# MSV000089558
data_TIC_MSV000089558 %>% dplyr::filter(!(str_detect(pattern = regex("QC|Blank", ignore_case = TRUE), SampleID))) %>% 
  dplyr::mutate(Order = seq_len(n())) %>% # fake order because it was not provided
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 4e9) +
  stat_cor()

data_TIC_MSV000089558 %>% dplyr::filter(str_detect(pattern = regex("QC", ignore_case = TRUE), SampleID)) %>% 
  dplyr::mutate(Order = seq_len(n())) %>% # this order is correct according to SampleID
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 1.5e9) +
  stat_cor()

data_TIC_MSV000089558 %>% dplyr::filter(str_detect(pattern = regex("Blank", ignore_case = TRUE), SampleID)) %>% 
  dplyr::mutate(Order = seq_len(n())) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 1e9) +
  stat_cor()

#MSV000092652
data_TIC_MSV000092652 %>% dplyr::filter(!(str_detect(pattern = regex("QC|Blank", ignore_case = TRUE), SampleID))) %>% 
  dplyr::mutate(Order = seq_len(n())) %>%
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 2.5e9) +
  stat_cor()

data_TIC_MSV000092652 %>% dplyr::filter(str_detect(pattern = regex("QC", ignore_case = TRUE), SampleID)) %>% 
  dplyr::filter(!(str_detect(pattern = "_00", SampleID))) %>% # samples with similar ID are the problem (remove)
  dplyr::mutate(Order = seq_len(n())) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 3e9) +
  stat_cor()

data_TIC_MSV000092652 %>% dplyr::filter(str_detect(pattern = regex("Blank", ignore_case = TRUE), SampleID)) %>% 
  dplyr::mutate(Order = seq_len(n())) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 7e8) +
  stat_cor()


# Check internal standard 

# MSV000089558
fbmn_IS_MSV000089558 <- annotations_MSV000089558 %>% 
  dplyr::filter(str_detect(Compound_Name, regex("sulf", ignore_case = TRUE))) %>% 
  distinct(X.Scan., .keep_all = TRUE) %>% dplyr::filter(Organism != "BILELIB19")

# Extract IS features for the table
table_IS_MSV000089558 <- data_MSV000089558_feces %>% column_to_rownames("SampleID") %>% 
  t() %>% as.data.frame() %>% rownames_to_column("ID") %>% dplyr::filter(ID %in% fbmn_IS_MSV000089558$X.Scan.) %>% 
  column_to_rownames("ID") %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% 
  dplyr::filter(!(str_detect(Sample, pattern = regex("QC|Blank", ignore_case = TRUE)))) %>%
  dplyr::select(Sample, `17396`) 
# Removed 2 samples with high IS and one with no (ME_C8-B1, ME_C10-B14, ME_C8-B14)

colnames(table_IS_MSV000089558)[2] <- "Sulfadimethoxine"

table_IS_MSV000089558 %>% 
  dplyr::mutate(Order = seq_len(n())) %>% # fake order
  ggscatter(x = "Order", y = "Sulfadimethoxine", add = c("reg.line")) +
  ylim(0, 5e6) + stat_cor()

table_IS_MSV000089558 %>% 
  dplyr::mutate(Order = seq_len(n())) %>%
  ggbarplot(x = "Order", y = "Sulfadimethoxine", xlab = "Run Order", 
            ylab = "Peak Area Sulfadimethoxine", title = "Internal Standard Acquisition") +
  geom_hline(yintercept = mean(table_IS_MSV000089558$Sulfadimethoxine, na.rm = TRUE), linetype = "dashed", color = "blue")

cv_is <- sd(table_IS_MSV000089558$Sulfadimethoxine)/mean(table_IS_MSV000089558$Sulfadimethoxine)

# MSV000092652
fbmn_IS_MSV000092652 <- annotations_MSV000092652 %>% 
  dplyr::filter(str_detect(Compound_Name, regex("sulf", ignore_case = TRUE))) %>% 
  distinct(X.Scan., .keep_all = TRUE) %>% dplyr::filter(Organism != "BILELIB19")

# Extract IS features for the table
table_IS_MSV000092652 <- data_MSV000092652_feces %>% column_to_rownames("SampleID") %>% 
  t() %>% as.data.frame() %>% rownames_to_column("ID") %>% dplyr::filter(ID %in% fbmn_IS_MSV000092652$X.Scan.) %>% 
  column_to_rownames("ID") %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% 
  dplyr::filter(!(str_detect(Sample, pattern = regex("QC|Blank", ignore_case = TRUE)))) %>%
  dplyr::select(Sample, `17160`) 

colnames(table_IS_MSV000092652)[2] <- "Sulfadimethoxine"

table_IS_MSV000092652 %>% 
  dplyr::mutate(Order = seq_len(n())) %>% # fake order
  ggscatter(x = "Order", y = "Sulfadimethoxine", add = c("reg.line")) +
  ylim(0, 1.5e7) + stat_cor()

table_IS_MSV000092652 %>% 
  dplyr::mutate(Order = seq_len(n())) %>% # fake order 
  ggbarplot(x = "Order", y = "Sulfadimethoxine", xlab = "Run Order", 
            ylab = "Peak Area Sulfadimethoxine", title = "Internal Standard Acquisition") +
  geom_hline(yintercept = mean(table_IS_MSV000092652$Sulfadimethoxine, na.rm = TRUE), linetype = "dashed", color = "blue")

cv_is <- sd(table_IS_MSV000092652$Sulfadimethoxine)/mean(table_IS_MSV000092652$Sulfadimethoxine)


# Check features quality per dataset

# MSV000089558
data_blank_MSV000089558 <- data_MSV000089558_feces %>% dplyr::filter(str_detect(pattern = "Blank", SampleID))
data_sixmix_MSV000089558 <- data_MSV000089558_feces %>% dplyr::filter(str_detect(pattern = "Mix", SampleID))
data_sample_MSV000089558 <- data_MSV000089558_feces %>% dplyr::filter(!(str_detect(pattern = "Blank|Mix", SampleID)))

# Blank
blanks_feature_info_MSV000089558 <- data.frame(Feature = colnames(data_blank_MSV000089558)[-1],
                                               Mean_blank = data_blank_MSV000089558 %>% column_to_rownames("SampleID") %>% colMeans(), 
                                               SD_blank =  data_blank_MSV000089558 %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_blank = SD_blank/Mean_blank) %>% left_join(annotations_MSV000089558, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature_MSV000089558) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                                         Precursor_MZ, Mean_blank, SD_blank, CV_blank) %>% 
  dplyr::filter(Mean_blank > 0) %>% arrange(desc(Mean_blank))

# QCmix
sixmix_feature_info_MSV000089558 <- data.frame(Feature = colnames(data_sixmix_MSV000089558)[-1],
                                               Mean_sixmix = data_sixmix_MSV000089558 %>% column_to_rownames("SampleID") %>% colMeans(), 
                                               SD_sixmix = data_sixmix_MSV000089558 %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_sixmix = SD_sixmix/Mean_sixmix) %>% left_join(annotations_MSV000089558, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature_MSV000089558) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                                         Precursor_MZ, Mean_sixmix, SD_sixmix, CV_sixmix) %>% 
  dplyr::filter(Mean_sixmix > 0) %>% arrange(desc(Mean_sixmix))

# Sample
sample_feature_info_MSV000089558 <- data.frame(Feature = colnames(data_sample_MSV000089558)[-1],
                                               Mean_sample = data_sample_MSV000089558 %>% column_to_rownames("SampleID") %>% colMeans(), 
                                               SD_sample =  data_sample_MSV000089558 %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_sample = SD_sample/Mean_sample) %>% left_join(annotations_MSV000089558, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature_MSV000089558) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                                         Precursor_MZ, Mean_sample, SD_sample, CV_sample) %>% 
  dplyr::filter(Mean_sample > 0) %>% arrange(desc(Mean_sample))


# MSV000092652
data_blank_MSV000092652 <- data_MSV000092652_feces %>% dplyr::filter(str_detect(pattern = "Blank", SampleID))
data_sixmix_MSV000092652 <- data_MSV000092652_feces %>% dplyr::filter(str_detect(pattern = "Mix", SampleID)) %>%
  dplyr::filter(!(str_detect(pattern = "_00", SampleID)))
data_sample_MSV000092652 <- data_MSV000092652_feces %>% dplyr::filter(!(str_detect(pattern = "Blank|Mix", SampleID)))

# Blank
blanks_feature_info_MSV000092652 <- data.frame(Feature = colnames(data_blank_MSV000092652)[-1],
                                               Mean_blank = data_blank_MSV000092652 %>% column_to_rownames("SampleID") %>% colMeans(), 
                                               SD_blank =  data_blank_MSV000092652 %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_blank = SD_blank/Mean_blank) %>% left_join(annotations_MSV000092652, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature_MSV000092652) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                                         Precursor_MZ, Mean_blank, SD_blank, CV_blank) %>% 
  dplyr::filter(Mean_blank > 0) %>% arrange(desc(Mean_blank))

# QCmix
sixmix_feature_info_MSV000092652 <- data.frame(Feature = colnames(data_sixmix_MSV000092652)[-1],
                                               Mean_sixmix = data_sixmix_MSV000092652 %>% column_to_rownames("SampleID") %>% colMeans(), 
                                               SD_sixmix = data_sixmix_MSV000092652 %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_sixmix = SD_sixmix/Mean_sixmix) %>% left_join(annotations_MSV000092652, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature_MSV000092652) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                                         Precursor_MZ, Mean_sixmix, SD_sixmix, CV_sixmix) %>% 
  dplyr::filter(Mean_sixmix > 0) %>% arrange(desc(Mean_sixmix))

# Sample
sample_feature_info_MSV000092652 <- data.frame(Feature = colnames(data_sample_MSV000092652)[-1],
                                               Mean_sample = data_sample_MSV000092652 %>% column_to_rownames("SampleID") %>% colMeans(), 
                                               SD_sample =  data_sample_MSV000092652 %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_sample = SD_sample/Mean_sample) %>% left_join(annotations_MSV000092652, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature_MSV000092652) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                                         Precursor_MZ, Mean_sample, SD_sample, CV_sample) %>% 
  dplyr::filter(Mean_sample > 0) %>% arrange(desc(Mean_sample))


#################
# Blank removal #
#################

# MSV000089558
# Features to be removed Sample/Blank < 5
feature_to_remove_MSV000089558 <- blanks_feature_info_MSV000089558 %>% 
  left_join(sample_feature_info_MSV000089558) %>%
  dplyr::filter(Mean_blank > 0) %>% 
  dplyr::mutate(Sample_Blank = Mean_sample/Mean_blank) %>% 
  dplyr::filter(Sample_Blank < 5 | is.na(Sample_Blank))

# Data with blank removal
data_clean_MSV000089558 <- data_MSV000089558_feces %>% 
  dplyr::select(-c(feature_to_remove_MSV000089558$Feature)) %>%
  column_to_rownames("SampleID") %>%
  select_if(~sum(.) != 0) %>% rownames_to_column("SampleID")

# Features to be removed Sample/QCmix < 5
feature_to_remove_mix_MSV000089558 <- sixmix_feature_info_MSV000089558 %>% 
  left_join(sample_feature_info_MSV000089558) %>% 
  dplyr::filter(Mean_sixmix > 0) %>% 
  dplyr::mutate(Sample_Mix = Mean_sample/Mean_sixmix) %>% 
  dplyr::filter(Sample_Mix < 5 | is.na(Sample_Mix)) %>% 
  dplyr::filter(!(Feature %in% feature_to_remove_MSV000089558$Feature))

# Data with qc blank removal
data_clean2_MSV000089558 <- data_clean_MSV000089558 %>% 
  dplyr::select(-c(feature_to_remove_mix_MSV000089558$Feature))

# Remove feature before 0.2 minutes and after 8 minutes
feature_to_rt_MSV000089558 <- info_feature_MSV000089558_complete %>% 
  dplyr::filter(RT < 0.2 | RT > 8) %>%
  dplyr::filter(!(Feature %in% feature_to_remove_MSV000089558$Feature)) %>%
  dplyr::filter(!(Feature %in% feature_to_remove_mix_MSV000089558$Feature)) %>%
  dplyr::filter(Feature %in% colnames(data_clean2_MSV000089558))

# Final cleaned table
data_clean3_MSV000089558 <- data_clean2_MSV000089558 %>% 
  dplyr::select(-c(feature_to_rt_MSV000089558$Feature))


# MSV000092652
# Features to be removed Sample/Blank < 5
feature_to_remove_MSV000092652 <- blanks_feature_info_MSV000092652 %>% 
  left_join(sample_feature_info_MSV000092652) %>%
  dplyr::filter(Mean_blank > 0) %>% 
  dplyr::mutate(Sample_Blank = Mean_sample/Mean_blank) %>% 
  dplyr::filter(Sample_Blank < 5 | is.na(Sample_Blank))

# Data with blank removal
data_clean_MSV000092652 <- data_MSV000092652_feces %>% 
  dplyr::select(-c(feature_to_remove_MSV000092652$Feature)) %>%
  column_to_rownames("SampleID") %>%
  select_if(~sum(.) != 0) %>% rownames_to_column("SampleID")

# Features to be removed Sample/QCmix < 5
feature_to_remove_mix_MSV000092652 <- sixmix_feature_info_MSV000092652 %>% 
  left_join(sample_feature_info_MSV000092652) %>% 
  dplyr::filter(Mean_sixmix > 0) %>% 
  dplyr::mutate(Sample_Mix = Mean_sample/Mean_sixmix) %>% 
  dplyr::filter(Sample_Mix < 5 | is.na(Sample_Mix)) %>% 
  dplyr::filter(!(Feature %in% feature_to_remove_MSV000092652$Feature))

# Data with qc blank removal
data_clean2_MSV000092652 <- data_clean_MSV000092652 %>% 
  dplyr::select(-c(feature_to_remove_mix_MSV000092652$Feature))

# Remove feature before 0.2 minutes and after 8 minutes
feature_to_rt_MSV000092652 <- info_feature_MSV000092652_complete %>% 
  dplyr::filter(RT < 0.2 | RT > 8) %>%
  dplyr::filter(!(Feature %in% feature_to_remove_MSV000092652$Feature)) %>%
  dplyr::filter(!(Feature %in% feature_to_remove_mix_MSV000092652$Feature)) %>%
  dplyr::filter(Feature %in% colnames(data_clean2_MSV000092652))

# Final cleaned table
data_clean3_MSV000092652 <- data_clean2_MSV000092652 %>% 
  dplyr::select(-c(feature_to_rt_MSV000092652$Feature))


#################
# PCA - Cleaned #
#################

# MSV000089558
PCA_raw <- mixOmics::pca(data_clean3_MSV000089558 %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
                           column_to_rownames("SampleID"), ncomp = 2, center = TRUE, scale = TRUE)
PCA_raw_scores <- data.frame(PCA_raw$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_filter, by = c("SampleID" = "Sample"))

i <- "AMP"

PCA_raw_plot <- PCA_raw_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - ", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_raw$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_raw$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_raw_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# MSV000092652
PCA_raw <- mixOmics::pca(data_clean3_MSV000092652 %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
                           column_to_rownames("SampleID"), ncomp = 2, center = TRUE, scale = TRUE)
PCA_raw_scores <- data.frame(PCA_raw$variates$X) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_filter, by = c("SampleID" = "Sample"))

i <- "AMP"

PCA_raw_plot <- PCA_raw_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - ", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_raw$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_raw$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_raw_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

#####################
# Keep only samples #
#####################

# MSV000089558
data_sample_MSV000089558 <- data_clean3_MSV000089558 %>% 
  dplyr::filter(!(str_detect(pattern = "Mix|Blank", SampleID))) %>%
  dplyr::filter(SampleID != "ME_C10_B1") # removed because clustered separately from other animals receiving PBS

# RCLR transformation
data_sample_clr_MSV000089558 <- decostand(data_sample_MSV000089558 %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sample_clr_MSV000089558 %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores_MSV000089558 <- data.frame(PCA_whole$variates$X) %>% rownames_to_column("Sample") %>% 
  left_join(metadata_filter) %>% dplyr::mutate_at("Days_Post_Birth", as.factor)

PCA_MSV000089558_plots <- list()

for (i in c("Days_Post_Birth", "AMP", "Cage")) {
  
  PCA_plot <- PCA_whole_scores_MSV000089558 %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores_MSV000089558 %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_MSV000089558_plots[[i]] <- PCA_plot
  
}

PCA_MSV000089558_plots_final <- wrap_plots(PCA_MSV000089558_plots, nrow = 1)

i <- "Days_Post_Birth"

PCA_plot_hypo1_time <- PCA_whole_scores_MSV000089558 %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, title = "PCA - Prepartum", legend = "none",
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_scores_MSV000089558 %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6)) + scale_color_viridis_d()

#ggsave(filename = "PCA_hypo1_time.svg", plot = PCA_plot_hypo1_time, device = "svg", dpi = "retina", width = 2.5, height = 2.5)

# PERMANOVA
dist_metabolites_MSV000089558 <- vegdist(data_sample_clr_MSV000089558, method = "euclidean")
disper_amp <- betadisper(dist_metabolites_MSV000089558, PCA_whole_scores_MSV000089558$AMP)
anova(disper_amp)
permanova_MSV000089558 <- adonis2(dist_metabolites_MSV000089558 ~ AMP + Days_Post_Birth + Cage, 
                                  PCA_whole_scores_MSV000089558, na.action = na.omit, by = "terms")

# No AMP effect at Day -3 --> animal did not receive AMP
# AMP effect can be noticed between days -2 and 2
# Little effect seems to be present from day 7 onward
# PERMANOVA confirms AMP and time effect


# MSV000092652
data_sample_MSV000092652 <- data_clean3_MSV000092652 %>% 
  dplyr::filter(!(str_detect(pattern = "Mix|Blank", SampleID)))

# RCLR transformation
data_sample_clr_MSV000092652 <- decostand(data_sample_MSV000092652 %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sample_clr_MSV000092652 %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores_MSV000092652 <- data.frame(PCA_whole$variates$X) %>% rownames_to_column("Sample") %>% 
  left_join(metadata_filter) %>% dplyr::mutate_at("Days_Post_Birth", as.factor)

PCA_MSV000092652_plots <- list()

for (i in c("Days_Post_Birth", "AMP", "Cage")) {
  
  PCA_plot <- PCA_whole_scores_MSV000092652 %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores_MSV000092652 %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_MSV000092652_plots[[i]] <- PCA_plot
  
}

PCA_MSV000092652_plots_final <- wrap_plots(PCA_MSV000092652_plots, nrow = 1)

i <- "Days_Post_Birth"

PCA_plot_hypo2_time <- PCA_whole_scores_MSV000092652 %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, title = "PCA - Postpartum", legend = "none",
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_scores_MSV000092652 %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6)) + scale_color_viridis_d()

#ggsave(filename = "PCA_hypo2_time.svg", plot = PCA_plot_hypo2_time, device = "svg", dpi = "retina", width = 2.5, height = 2.5)

# PERMANOVA
dist_metabolites_MSV000092652 <- vegdist(data_sample_clr_MSV000092652, method = "euclidean")
disper_amp <- betadisper(dist_metabolites_MSV000092652, PCA_whole_scores_MSV000092652$AMP)
anova(disper_amp)
permanova_MSV000092652 <- adonis2(dist_metabolites_MSV000092652 ~ AMP + Days_Post_Birth + Cage, 
                                  PCA_whole_scores_MSV000092652, na.action = na.omit, by = "terms")

# No AMP effect at Day 1 --> animals did not receive AMP
# AMP effect observable between days 2 and 3
# Little effect present from day 7 onward
# Animal separation based on treatment is less apparent compared to antepartum
# PERMANOVA confirms AMP and time effect


##########
# WEIGHT #
##########

# Check weight at weaning. Hypo3 is not taken into consideration in this case.
# Hypo1 was re-run. Keep only second batch since the first batch of animals had a
# very low weight --> usual weight of animals at PND21 is ~12g

weight_plot <- data_weight %>% dplyr::filter(Hypothesis != 3) %>% 
  dplyr::filter(!(Hypothesis == 1 & Batch == 1)) %>% 
  ggboxplot(x = "AMP", y = "Weight", facet.by = "Hypothesis", add = "jitter", 
            add.params = list(color = "AMP", alpha = 0.8, size = 1),
            palette = c("#287DAB", "#E5BF86"), title = "Pups weight at weaning",
            xlab = FALSE, ylab = "Weight (g)", legend = "none") + ylim(0, 25) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

#ggsave(plot = weight_plot, filename = "weight_pups.svg", device = "svg", dpi = "retina", width = 3, height = 2)

model1 <- data_weight %>% dplyr::filter(Hypothesis != 3) %>%
  dplyr::filter(!(Hypothesis == 1 & Batch == 1)) %>%
  dplyr::filter(Hypothesis == 1) %>%
  lm(formula = Weight ~ AMP + Sex + Cage)

summary(model1)

model2 <- data_weight %>% dplyr::filter(Hypothesis != 3) %>%
  dplyr::filter(!(Hypothesis == 1 & Batch == 1)) %>%
  dplyr::filter(Hypothesis == 2) %>%
  lm(formula = Weight ~ AMP + Sex + Cage)

summary(model2)


##################
# AMP DETECTTION #
##################

amp_MSV000089558 <- data_sample_MSV000089558 %>% 
  dplyr::filter(SampleID %in% (metadata_filter %>% dplyr::filter(AMP == "TRUE"))$Sample) %>%
  dplyr::mutate(Row_Sum = rowSums(select(., -SampleID), na.rm = TRUE)) %>%
  dplyr::select(SampleID, "13912", "10689", "17242", Row_Sum) %>%
  left_join(metadata_filter, by = c("SampleID" = "Sample")) %>%
  dplyr::mutate(LogAmp = log(`13912` +1)) %>%
  dplyr::mutate(Amp_RA = `13912`/Row_Sum)

amp_MSV000089558_plot <- amp_MSV000089558 %>% 
  ggboxplot(x = "Days_Post_Birth", y = "Amp_RA", add = "jitter", 
            legend = "none", ylab = "Relative Abundance", title = "Fecal AMP - Prepartum",
            xlab = "Day", add.params = list(color = "AMP", alpha = 0.6, size = 1)) +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

amp_MSV000092652 <- data_sample_MSV000092652 %>% 
  dplyr::filter(SampleID %in% (metadata_filter %>% dplyr::filter(AMP == "TRUE"))$Sample) %>%
  dplyr::mutate(Row_Sum = rowSums(select(., -SampleID), na.rm = TRUE)) %>%
  dplyr::select(SampleID, "13409", "9515", "17000", Row_Sum) %>%
  left_join(metadata_filter, by = c("SampleID" = "Sample")) %>%
  dplyr::mutate(LogAmp = log(`13409` + 1)) %>%
  dplyr::mutate(Amp_RA = `13409`/Row_Sum)

amp_MSV000092652_plot <- amp_MSV000092652 %>% 
  ggboxplot(x = "Days_Post_Birth", y = "Amp_RA", add = "jitter", 
            legend = "none", ylab = "Relative Abundance", title = "Fecal AMP - Postpartum",
            xlab = "Day", add.params = list(color = "AMP", alpha = 0.6, size = 1)) +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))
# no AMP detected in one sample at day 3 (ME_C6_B3)

# Combine AMP plots 
combined_amp_plot <- wrap_plots(amp_MSV000089558_plot, amp_MSV000092652_plot)

#ggsave(filename = "AMP_detection.svg", plot = combined_amp_plot, device = "svg", dpi = "retina", height = 2, width = 4)

# As AMP was administered on PND-2 and PND-1 (GD17 and GD18) in the antepartum cohort and on 
# PND2 and PND3 in the postpartum cohort, extract features affected on those 
# specific timepoints as it was done for the microbiome analysis

prepartum_sample <- metadata %>% dplyr::filter(Hypothesis == 1 & Batch == 1 & Days_Post_Birth %in% c(-2, -1))
postpartum_sample <- metadata %>% dplyr::filter(Hypothesis == 2 & Batch == 1 & Days_Post_Birth %in% c(2, 3))

# Filter data
data_sample_MSV000089558_amp <- data_sample_MSV000089558 %>% dplyr::filter(SampleID %in% prepartum_sample$Sample)
data_sample_MSV000092652_amp <- data_sample_MSV000092652 %>% dplyr::filter(SampleID %in% postpartum_sample$Sample)


#############################
# AMP EFFECT - MSV000089558 #
#############################

# RCLR transformation
data_sample_clr_MSV000089558_amp <- decostand(data_sample_MSV000089558_amp %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sample_clr_MSV000089558_amp %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores_MSV000089558_amp <- data.frame(PCA_whole$variates$X) %>% rownames_to_column("Sample") %>% 
  left_join(metadata) %>% dplyr::mutate_at("Days_Post_Birth", as.factor)

PCA_MSV000089558_amp_plots <- list()

for (i in c("Days_Post_Birth", "AMP", "Cage")) {
  
  PCA_plot <- PCA_whole_scores_MSV000089558_amp %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, 
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores_MSV000089558_amp %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_MSV000089558_plots[[i]] <- PCA_plot
  
}

PCA_MSV000089558_plots_final <- wrap_plots(PCA_MSV000089558_plots, nrow = 1)

i <- "AMP"

PCA_plot_MSV000089558_amp <- PCA_whole_scores_MSV000089558_amp %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, shape = "Days_Post_Birth",
            title = "PCA - Prepartum Exposure", palette = c("#287DAB", "#E5BF86"), legend ="none",
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_scores_MSV000089558_amp %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

#ggsave(filename = "PCA_Hypo1_amp.svg", plot = PCA_plot_MSV000089558_amp, device = "svg", dpi = "retina", width = 2, height = 2)

# PERMANOVA
dist_metabolites_MSV000089558_amp <- vegdist(data_sample_clr_MSV000089558_amp, method = "euclidean")
disper_amp_ealry <- betadisper(dist_metabolites_MSV000089558_amp, PCA_whole_scores_MSV000089558_amp$AMP)
anova(disper_amp_ealry)
permanova_MSV000089558_amp <- adonis2(dist_metabolites_MSV000089558_amp ~ AMP + Days_Post_Birth + Cage, 
                                      PCA_whole_scores_MSV000089558_amp, na.action = na.omit, by = "terms")

# PLS-DA
PLSDA_MSV000089558_amp <- mixOmics::plsda(data_sample_clr_MSV000089558_amp %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                          PCA_whole_scores_MSV000089558_amp$AMP, ncomp = 2, scale = TRUE)
PLSDA_scores_MSV000089558_amp <- data.frame(PLSDA_MSV000089558_amp$variates$X) %>% 
  rownames_to_column("Sample") %>% left_join(metadata_filter)

PLSDA_plot_MSV000089558_amp <- PLSDA_scores_MSV000089558_amp %>%
  ggscatter(x = "comp1", y = "comp2", color = "AMP", alpha = 0.6, title = "PLSDA - Prepartum Exposure",
            xlab = paste("Component 1 (", round(PLSDA_MSV000089558_amp$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_MSV000089558_amp$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_scores_MSV000089558_amp %>% group_by(AMP) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = AMP), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_MSV000089558_amp <- mixOmics::plotLoadings(PLSDA_MSV000089558_amp, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_MSV000089558_amp <- mixOmics::perf(PLSDA_MSV000089558_amp, validation = "loo") 
#plot(perf_plsda_MSV000089558_amp, legend = FALSE)

VIPs_MSV000089558_amp <- as.data.frame(mixOmics::vip(PLSDA_MSV000089558_amp))
VIPs_MSV000089558_amp_filter <- dplyr::filter(VIPs_MSV000089558_amp, VIPs_MSV000089558_amp$comp1 > 1)
VIPs_MSV000089558_amp_filter$ID <- rownames(VIPs_MSV000089558_amp_filter)
VIPs_MSV000089558_amp_select <- VIPs_MSV000089558_amp_filter %>% dplyr::select(ID, comp1)
VIPs_MSV000089558_amp_Load <- VIPs_MSV000089558_amp_select %>% 
  left_join(Loadings_MSV000089558_amp, by = c("ID" = "rowname")) %>% 
  left_join(info_feature_MSV000089558_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1)) %>%
  dplyr::mutate(dataset = paste("MSV000089558", ID, sep = "_"))

#write_csv(x = VIPs_MSV000089558_amp_Load, file = "Table_S5.csv")


#############################
# AMP EFFECT - MSV000092652 #
#############################

# RCLR transformation
data_sample_clr_MSV000092652_amp <- decostand(data_sample_MSV000092652_amp %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sample_clr_MSV000092652_amp %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores_MSV000092652_amp <- data.frame(PCA_whole$variates$X) %>% rownames_to_column("Sample") %>% 
  left_join(metadata)  %>% dplyr::mutate_at("Days_Post_Birth", as.factor)

PCA_MSV000092652_amp_plots <- list()

for (i in c("Days_Post_Birth", "AMP", "Cage")) {
  
  PCA_plot <- PCA_whole_scores_MSV000092652_amp %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores_MSV000092652_amp %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_MSV000092652_plots[[i]] <- PCA_plot
  
}

PCA_MSV000092652_plots_final <- wrap_plots(PCA_MSV000092652_plots, nrow = 1)

i <- "AMP"

PCA_plot_MSV000092652_amp <- PCA_whole_scores_MSV000092652_amp %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, shape = "Days_Post_Birth",
            title = "PCA - Prepartum Exposure", palette = c("#287DAB", "#E5BF86"), legend ="none",
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_scores_MSV000092652_amp %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

#ggsave(filename = "PCA_Hypo2_amp.svg", plot = PCA_plot_MSV000092652_amp, device = "svg", dpi = "retina", width = 2, height = 2)

# PERMANOVA
dist_metabolites_MSV000092652_amp <- vegdist(data_sample_clr_MSV000092652_amp, method = "euclidean")
disper_amp_ealry <- betadisper(dist_metabolites_MSV000092652_amp, PCA_whole_scores_MSV000092652_amp$AMP)
anova(disper_amp_ealry)
permanova_MSV000092652_amp <- adonis2(dist_metabolites_MSV000092652_amp ~ AMP + Days_Post_Birth + Cage, 
                                      PCA_whole_scores_MSV000092652_amp, na.action = na.omit, by = "terms")

# For the PLS-DA model exclude the three samples clustering separately to 
# maximize feature detection. Samples will be included in the final ratio calculation

data_sample_clr_MSV000092652_amp_filter <- data_sample_clr_MSV000092652_amp %>% 
  rownames_to_column("SampleID") %>%
  dplyr::filter(!(SampleID %in% c("ME_C6_B3", "ME_C1_B3", "ME_C5_B2"))) %>%
  column_to_rownames("SampleID")

PCA_whole_scores_MSV000092652_amp_filter <- PCA_whole_scores_MSV000092652_amp %>%
  dplyr::filter(!(Sample %in% c("ME_C6_B3", "ME_C1_B3", "ME_C5_B2")))

# PLS-DA
PLSDA_MSV000092652_amp <- mixOmics::plsda(data_sample_clr_MSV000092652_amp_filter %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                          PCA_whole_scores_MSV000092652_amp_filter$AMP, ncomp = 2, scale = TRUE)
PLSDA_scores_MSV000092652_amp <- data.frame(PLSDA_MSV000092652_amp$variates$X) %>% 
  rownames_to_column("Sample") %>% left_join(metadata_filter)

PLSDA_plot_MSV000092652_amp <- PLSDA_scores_MSV000092652_amp %>%
  ggscatter(x = "comp1", y = "comp2", color = "AMP", alpha = 0.6, title = "PLSDA - Hypo 1 - AMP",
            xlab = paste("Component 1 (", round(PLSDA_MSV000092652_amp$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_MSV000092652_amp$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic()) +
  geom_point(data = PLSDA_scores_MSV000092652_amp %>% group_by(AMP) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = AMP), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_MSV000092652_amp <- mixOmics::plotLoadings(PLSDA_MSV000092652_amp, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib)

#perf_plsda_MSV000092652_amp <- mixOmics::perf(PLSDA_MSV000092652_amp, validation = "loo") 
#plot(perf_plsda_MSV000092652_amp, legend = FALSE)

VIPs_MSV000092652_amp <- as.data.frame(mixOmics::vip(PLSDA_MSV000092652_amp))
VIPs_MSV000092652_amp_filter <- dplyr::filter(VIPs_MSV000092652_amp, VIPs_MSV000092652_amp$comp1 > 1)
VIPs_MSV000092652_amp_filter$ID <- rownames(VIPs_MSV000092652_amp_filter)
VIPs_MSV000092652_amp_select <- VIPs_MSV000092652_amp_filter %>% dplyr::select(ID, comp1)
VIPs_MSV000092652_amp_Load <- VIPs_MSV000092652_amp_select %>% 
  left_join(Loadings_MSV000092652_amp, by = c("ID" = "rowname")) %>% 
  left_join(info_feature_MSV000092652_complete, by = c("ID" = "Feature")) %>% arrange(desc(comp1)) %>%
  dplyr::mutate(dataset = paste("MSV000092652", ID, sep = "_"))

#write_csv(x = VIPs_MSV000092652_amp_Load, file = "Table_S6.csv")


############################
# Check ratios across time #
############################

MSV000089558_up_features <- VIPs_MSV000089558_amp_Load %>% dplyr::filter(GroupContrib == TRUE) 
MSV000089558_down_features <- VIPs_MSV000089558_amp_Load %>% dplyr::filter(GroupContrib == FALSE) 

MSV000092652_up_features <- VIPs_MSV000092652_amp_Load %>% dplyr::filter(GroupContrib == TRUE) 
MSV000092652_down_features <- VIPs_MSV000092652_amp_Load %>% dplyr::filter(GroupContrib == FALSE) 

MSV000089558_ratio <- data_sample_MSV000089558 %>% 
  dplyr::select(SampleID, MSV000089558_up_features$ID, MSV000089558_down_features$ID) %>%
  dplyr::mutate(AMP_Yes = rowSums(select(., MSV000089558_up_features$ID))) %>%
  dplyr::mutate(AMP_No = rowSums(select(., MSV000089558_down_features$ID))) %>%
  dplyr::mutate(Ratio = log(AMP_Yes/AMP_No)) %>% 
  dplyr::select(SampleID, Ratio) %>%
  left_join(metadata_filter, by = c("SampleID" = "Sample"))

MSV000092652_ratio <- data_sample_MSV000092652 %>% 
  dplyr::select(SampleID, MSV000092652_up_features$ID, MSV000092652_down_features$ID) %>%
  dplyr::mutate(AMP_Yes = rowSums(select(., MSV000092652_up_features$ID))) %>%
  dplyr::mutate(AMP_No = rowSums(select(., MSV000092652_down_features$ID))) %>%
  dplyr::mutate(Ratio = log(AMP_Yes/AMP_No)) %>%
  dplyr::select(SampleID, Ratio) %>%
  left_join(metadata_filter, by = c("SampleID" = "Sample"))

plot_MSV000089558_ratio <- MSV000089558_ratio %>% 
  mutate(Day = case_when(Days_Post_Birth == "-3" ~ 1, Days_Post_Birth == "-2" ~ 2, Days_Post_Birth == "-1" ~ 3,
                         Days_Post_Birth == "1" ~ 4, Days_Post_Birth == "2" ~ 5, Days_Post_Birth == "7" ~ 10,
                         Days_Post_Birth == "14" ~ 17, Days_Post_Birth == "21" ~ 24)) %>% 
  mutate_at("Day", as.numeric) %>%
  ggscatter(x = "Day", y = "Ratio", color = "AMP", add = "loess", ylab = "Ln(AMP/PBS)",
            legend = "none", palette = c("#287DAB", "#E5BF86"), alpha = 0.5,
            title = "Differental Ratio - Prepartum") + 
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 10, 17, 24)) + 
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

plot_MSV000092652_ratio <- MSV000092652_ratio %>% 
  dplyr::filter(!(SampleID %in% c("ME_C6_B3", "ME_C5_B2"))) %>%
  mutate(Day = case_when(Days_Post_Birth == "1" ~ 4, Days_Post_Birth == "2" ~ 5, Days_Post_Birth == "3" ~ 6, 
                         Days_Post_Birth == "7" ~ 10, Days_Post_Birth == "14" ~ 17, Days_Post_Birth == "21" ~ 24)) %>% 
  mutate_at("Day", as.numeric) %>%
  ggscatter(x = "Day", y = "Ratio", color = "AMP", add = "loess", ylab = "Ln(AMP/PBS)",
            legend = "none", palette = c("#287DAB", "#E5BF86"), alpha = 0.5, add.params = list(span = 1.25),
            title = "Differental Ratio - Postpartum") +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 10, 17, 24), limits = c(1, 24)) + 
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

combined_ratio_plots <- wrap_plots(plot_MSV000089558_ratio, plot_MSV000092652_ratio, nrow = 2)

#ggsave(plot = combined_ratio_plots, filename = "Ratio_metabolomics_time.svg", device = "svg", dpi = "retina", width = 3, height = 3)


# Test each timepoint
p_value_MSV000089558 <- MSV000089558_ratio %>%
  dplyr::select(Ratio, Days_Post_Birth, AMP) %>%
  group_by(Days_Post_Birth) %>%
  nest() %>%
  dplyr::mutate(t_test = map(data, ~t.test(Ratio ~ AMP, data = .x)),
                p_value = map_dbl(t_test, ~.x$p.value)) %>%
  ungroup() %>%
  dplyr::mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
  dplyr::select(Days_Post_Birth, p_value, p_adjust) %>% arrange(Days_Post_Birth)

p_value_MSV000092652 <- MSV000092652_ratio %>%
  dplyr::filter(!(SampleID %in% c("ME_C6_B3", "ME_C5_B2"))) %>%
  dplyr::select(Ratio, Days_Post_Birth, AMP) %>%
  group_by(Days_Post_Birth) %>%
  nest() %>%
  dplyr::mutate(t_test = map(data, ~t.test(Ratio ~ AMP, data = .x)),
                p_value = map_dbl(t_test, ~.x$p.value)) %>%
  ungroup() %>%
  dplyr::mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
  dplyr::select(Days_Post_Birth, p_value, p_adjust) %>% arrange(Days_Post_Birth)


###################################
# Trace molecules across datasets #
###################################

clustering_mn <- read.delim("data_metabolomics/clusterinfo.tsv") %>% 
  dplyr::filter(str_detect(pattern = "MSV000089558|MSV000092652", X.Filename)) %>%
  dplyr::select(-8)
colnames(clustering_mn) <- c("clusterIDX", "filename", "specIDX", "Scan", "mz", "charge", "rt")

clustering_mn_match <- clustering_mn %>% 
  dplyr::mutate(dataset = case_when(str_detect(pattern = "MSV000089558", filename) ~ paste("MSV000089558", Scan, sep = "_"),
                                    str_detect(pattern = "MSV000092652", filename) ~ paste("MSV000092652", Scan, sep = "_")))

network_pairs <- read.delim("data_metabolomics/merged_pairs.tsv")
network_components_info <- read.delim("data_metabolomics/clustersummary_with_network.tsv") %>% 
  dplyr::select(1,3,7,8)

network_pairs_info <- network_pairs %>% 
  dplyr::filter(abs(DeltaMZ) < 0.02) %>% # cause run with orbitrap
  left_join(clustering_mn_match %>% dplyr::select(clusterIDX, dataset, mz, rt), by = c("CLUSTERID1" = "clusterIDX")) %>% 
  left_join(clustering_mn_match %>% dplyr::select(clusterIDX, dataset, mz, rt), by = c("CLUSTERID2" = "clusterIDX")) %>%
  dplyr::mutate(DeltaRT = abs(rt.x - rt.y)) %>% 
  dplyr::filter(DeltaRT < 0.3) %>% # same method was used so RT should be similar
  dplyr::filter(Cosine > 0.7) %>% 
  dplyr::mutate(SameDataset = (gsub("_.*$", "", dataset.x) == gsub("_.*$", "", dataset.y))) %>%
  dplyr::filter(SameDataset == FALSE) %>% # I want to remove matches within the same dataset
  dplyr::filter(rt.x > 0.2 & rt.y > 0.2) %>%
  dplyr::filter(rt.x < 8 & rt.y < 8) %>% 
  dplyr::filter(mz.x < 1500 & mz.y < 1500) %>%
  dplyr::select(dataset.x, dataset.y, mz.x, mz.y, DeltaMZ, rt.x, rt.y, DeltaRT, Cosine, CLUSTERID1, CLUSTERID2)

network_pairs_info %>% ggdensity("DeltaMZ")
network_pairs_info %>% ggdensity("DeltaRT")
network_pairs_info %>% ggdensity(x = "Cosine")

network_pairs_info_interest <- network_pairs_info %>% 
  dplyr::filter(dataset.x %in% VIPs_MSV000089558_amp_Load$dataset) %>%
  dplyr::filter(dataset.y %in% VIPs_MSV000092652_amp_Load$dataset) %>%
  distinct(dataset.x, .keep_all = TRUE) %>% 
  distinct(dataset.y, .keep_all = TRUE) %>%
  dplyr::mutate(index_new = seq_len(n()))

# Extract overlap
VIPs_MSV000089558_amp_overlap <- network_pairs_info_interest %>% dplyr::select(dataset.x, index_new) %>%
  left_join(VIPs_MSV000089558_amp_Load, by = c("dataset.x" = "dataset"))
VIPs_MSV000092652_amp_overlap <- network_pairs_info_interest %>% dplyr::select(dataset.y, index_new) %>%
  left_join(VIPs_MSV000092652_amp_Load, by = c("dataset.y" = "dataset"))

combined_overlap_table <- VIPs_MSV000089558_amp_overlap %>% 
  dplyr::select(dataset.x, index_new, comp1, GroupContrib) %>%
  inner_join(VIPs_MSV000092652_amp_overlap %>% 
               dplyr::select(dataset.y, index_new, comp1, GroupContrib, mz, RT, 
                             Compound_Name, rev_cosine, ba_lib, syn_lib), by = "index_new") %>%
  dplyr::select(2,8,9,1,5,3,6,4,7,10,11,12,13) %>%
  dplyr::mutate(Concordance = case_when(GroupContrib.x == GroupContrib.y ~ "Yes",
                                        TRUE ~ "No")) %>%
  dplyr::mutate(mean_vip = rowMeans(across(c(comp1.x, comp1.y)), na.rm = TRUE)) %>%
  arrange(desc(mean_vip))

# Generate Upset plot
list_overlap <- list(
  `Prepartum PBS` = (VIPs_MSV000089558_amp_overlap %>% dplyr::filter(GroupContrib == "FALSE"))$index_new,
  `Prepartum AMP` =  (VIPs_MSV000089558_amp_overlap %>% dplyr::filter(GroupContrib == "TRUE"))$index_new,
  `Postpartum PBS` = (VIPs_MSV000092652_amp_overlap %>% dplyr::filter(GroupContrib == "FALSE"))$index_new,
  `Postpartum AMP` = (VIPs_MSV000092652_amp_overlap %>% dplyr::filter(GroupContrib == "TRUE"))$index_new)

overlap_upset <- UpSetR::upset(fromList(list_overlap), nsets = 8, nintersects = NA,
                               point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE,  
                               sets = c("Prepartum PBS", "Prepartum AMP", "Postpartum PBS", "Postpartum AMP"),
                               queries = list(list(query = intersects, params = list("Prepartum PBS", 
                                                                                     "Postpartum PBS"), color = "#287DAB", active = T), 
                                              list(query = intersects, params = list("Prepartum AMP", 
                                                                                     "Postpartum AMP"),color = "#E5BF86", active = T)))


######################
# CONCORDAT FEATURES #
######################

concordant <- combined_overlap_table %>%
  dplyr::filter(Concordance == "Yes")

# Generate an mgf of concordant spectra to be classified via CANOPUS
#dda <- Spectra("data_metabolomics/gnps_MSV000089558.mgf", source = MsBackendMgf())
#dda_ids <- data.frame(ID = dda@backend@spectraData@listData$FEATURE_ID) %>%
#  dplyr::mutate(Interest = ID %in% concordant$Feature.x)
#dda_filtered <- dda[dda_ids$Interest]
#export(dda_filtered, MsBackendMgf(), file = "data_metabolomics/concordant_mprint.mgf", exportTitle = FALSE)

# Add CANOPUS output
canopus_concordace <- read_tsv("data_metabolomics/canopus_concordance_mprint.tsv") %>%
  dplyr::mutate(Canopus = case_when(`NPC#class Probability` > 0.6 ~ `NPC#class`,
                                    `NPC#class Probability` < 0.6 & `NPC#superclass Probability` > 0.6 ~ `NPC#superclass`,
                                    TRUE ~ `NPC#pathway`)) %>%
  dplyr::select(mappingFeatureId, Canopus) %>%
  dplyr::mutate(mappingFeatureId = as.character(mappingFeatureId))

concordant_info <- concordant %>% 
  dplyr::mutate(mappingFeatureId = gsub("MSV000089558_", "", dataset.x)) %>%
  left_join(canopus_concordace)

#write_csv(x = concordant_info, file = "Table_S7.csv")

# Pie chart for NPC Pathway
concordant_amp_yes <- concordant_info %>% dplyr::filter(`GroupContrib.x` == "TRUE") %>%
  group_by(Canopus) %>% summarise(count_yes = n()) %>% arrange(desc(count_yes)) %>%
  dplyr::mutate(Canopus = case_when(count_yes < 10 ~ "zOther", TRUE ~ as.character(Canopus))) %>%
  group_by(Canopus) %>% summarise(count_yes = sum(count_yes)) %>%
  dplyr::mutate(ratio = count_yes/sum(count_yes)) %>%
  arrange(desc(ratio))

concordant_amp_no <- concordant_info %>% dplyr::filter(`GroupContrib.x` == "FALSE") %>% 
  group_by(Canopus) %>% summarise(count_no = n()) %>% arrange(desc(count_no)) %>%
  dplyr::mutate(Canopus = case_when(count_no < 10 ~ "zOther", TRUE ~ as.character(Canopus))) %>%
  group_by(Canopus) %>% summarise(count_no = sum(count_no)) %>%
  dplyr::mutate(ratio = count_no/sum(count_no)) %>%
  arrange(desc(ratio))

concordant_comb <- concordant_amp_yes %>% full_join(concordant_amp_no, by = "Canopus") %>%
  replace_na(list(count_yes = 0, count_no = 0))

pie_amp_yes <- concordant_comb %>% ggpie(x = "count_yes", fill = "Canopus", legend = "right") + scale_fill_viridis_d()
pie_amp_no <- concordant_comb %>% ggpie(x = "count_no", fill = "Canopus", legend = "right") + scale_fill_viridis_d()

pie_combined <- ggarrange(pie_amp_yes, pie_amp_no, common.legend = TRUE, legend = "right")

#ggsave(plot = pie_combined, filename = "Canopus_pies.svg", device = "svg")


##############
# CARNITINES #
##############

# Extract concordant carnitines
carn_conc <- concordant_info %>% 
  dplyr::filter(str_detect(pattern = "Carnitine", syn_lib) | str_detect(pattern = "carnitines", Canopus)) %>%
  dplyr::mutate(id = gsub(pattern = "MSV000089558_", "", dataset.x)) %>% dplyr::select(-index_new)

##############
# BILE ACIDS #
##############

# Use massql filter
no_ba <- ba_massql_pre %>% 
  dplyr::filter(query_validation == "Did not pass any selected query")

concordant_ba <- concordant_info %>% 
  dplyr::filter(str_detect(pattern = "pentanoic|cholic|bile", Compound_Name) | str_detect(pattern = "Cholane", Canopus)) %>%
  dplyr::filter(!(mappingFeatureId %in% no_ba$`#Scan#`)) %>%
  dplyr::mutate(mappingFeatureId = as.numeric(mappingFeatureId)) %>%
  left_join(ba_massql_pre %>% dplyr::select(1:2), by = c("mappingFeatureId" = "#Scan#")) %>%
  dplyr::filter(!is.na(query_validation))

# Generate ratio overtime
ba_up_MSV000089558 <- concordant_ba %>% dplyr::filter(GroupContrib.x == "TRUE") %>%
  dplyr::mutate(dataset.x = gsub("MSV000089558_", "", dataset.x))
ba_down_MSV000089558 <- concordant_ba %>% dplyr::filter(GroupContrib.x == "FALSE") %>%
  dplyr::mutate(dataset.x = gsub("MSV000089558_", "", dataset.x))

ratio_ba_MSV000089558 <- data_sample_MSV000089558 %>% 
  dplyr::select(SampleID, (concordant_ba %>% dplyr::mutate(dataset.x = gsub("MSV000089558_", "", dataset.x)))$`dataset.x`) %>%
  dplyr::mutate(AMP_Yes = rowSums(select(., ba_up_MSV000089558$dataset.x)) + 1100.115) %>%
  dplyr::mutate(AMP_No = rowSums(select(., ba_down_MSV000089558$dataset.x))) %>%
  dplyr::mutate(Ratio = log(AMP_No/AMP_Yes)) %>%
  left_join(metadata, by = c("SampleID" = "Sample")) %>% 
  mutate(Day = case_when(Days_Post_Birth == "-3" ~ 1, Days_Post_Birth == "-2" ~ 2, Days_Post_Birth == "-1" ~ 3,
                         Days_Post_Birth == "1" ~ 4, Days_Post_Birth == "2" ~ 5, Days_Post_Birth == "7" ~ 10,
                         Days_Post_Birth == "14" ~ 17, Days_Post_Birth == "21" ~ 24)) %>% 
  mutate_at("Day", as.numeric)

plot_ba_ratio_MSV000089558 <- ratio_ba_MSV000089558 %>%
  ggscatter(x = "Day", y = "Ratio", color = "AMP", add = "loess", ylab = "Ln(TopPBS/TopAMP)",
            legend = "none", palette = c("#287DAB", "#E5BF86"), alpha = 0.5,
            title = "Bile Acids Ratio - Prepartum") + 
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 10, 17, 24)) + 
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

ba_up_MSV000092652 <- concordant_ba %>% dplyr::filter(GroupContrib.y == "TRUE") %>%
  dplyr::mutate(dataset.y = gsub("MSV000092652_", "", dataset.y))
ba_down_MSV000092652 <- concordant_ba %>% dplyr::filter(GroupContrib.y == "FALSE") %>%
  dplyr::mutate(dataset.y = gsub("MSV000092652_", "", dataset.y))

ratio_ba_MSV000092652 <- data_sample_MSV000092652 %>% 
  dplyr::select(SampleID, (concordant_ba %>% dplyr::mutate(dataset.y = gsub("MSV000092652_", "", dataset.y)))$`dataset.y`) %>%
  dplyr::mutate(AMP_Yes = rowSums(select(., ba_up_MSV000092652$dataset.y)) + 872.1348) %>%
  dplyr::mutate(AMP_No = rowSums(select(., ba_down_MSV000092652$dataset.y))) %>%
  dplyr::mutate(Ratio = log(AMP_No/AMP_Yes)) %>%
  left_join(metadata, by = c("SampleID" = "Sample")) %>%
  mutate(Day = case_when(Days_Post_Birth == "1" ~ 4, Days_Post_Birth == "2" ~ 5, Days_Post_Birth == "3" ~ 6, 
                         Days_Post_Birth == "7" ~ 10, Days_Post_Birth == "14" ~ 17, Days_Post_Birth == "21" ~ 24)) %>% 
  mutate_at("Day", as.numeric)

plot_ba_ratio_MSV000092652 <- ratio_ba_MSV000092652 %>%
  ggscatter(x = "Day", y = "Ratio", color = "AMP", add = "loess", ylab = "Ln(TopPBS/TopAMP)",
            legend = "none", palette = c("#287DAB", "#E5BF86"), alpha = 0.5,
            title = "Bile Acids Ratio - Postpartum") + 
  scale_x_continuous(breaks = c(1, 4, 5, 6, 10, 17, 24), limits = c(1, 24)) + 
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

plot_ratio_ba_combined <- wrap_plots(plot_ba_ratio_MSV000089558, plot_ba_ratio_MSV000092652, nrow = 2)

#ggsave(plot = plot_ratio_ba_combined, filename = "Ratio_ba_time.svg", device = "svg", dpi = "retina", width = 3, height = 3)


# Test each timepoint 
p_value_MSV000089558_ba <- ratio_ba_MSV000089558 %>%
  dplyr::select(Ratio, Days_Post_Birth, AMP) %>%
  group_by(Days_Post_Birth) %>%
  nest() %>%
  dplyr::mutate(t_test = map(data, ~t.test(Ratio ~ AMP, data = .x)),
                p_value = map_dbl(t_test, ~.x$p.value)) %>%
  ungroup() %>%
  dplyr::mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
  dplyr::select(Days_Post_Birth, p_value, p_adjust) %>% arrange(Days_Post_Birth)

p_value_MSV000092652 <- ratio_ba_MSV000092652 %>%
  dplyr::filter(!(SampleID %in% c("ME_C6_B3", "ME_C5_B2"))) %>%
  dplyr::select(Ratio, Days_Post_Birth, AMP) %>%
  group_by(Days_Post_Birth) %>%
  nest() %>%
  dplyr::mutate(t_test = map(data, ~t.test(Ratio ~ AMP, data = .x, var.equal = FALSE)),
                p_value = map_dbl(t_test, ~.x$p.value)) %>%
  ungroup() %>%
  dplyr::mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
  dplyr::select(Days_Post_Birth, p_value, p_adjust) %>% arrange(Days_Post_Birth)


################
# Plot Adducts #
################
adducts <- read_csv("data_metabolomics/data_adducts.csv") %>%
  dplyr::filter(`RT (min)` < 4.4 & `RT (min)` > 4)

# Reshape to long format
adducts_long <- adducts %>%
  pivot_longer(cols = starts_with("Intensity_mz_"),
               names_to = "mz",
               values_to = "Intensity") %>%
  dplyr::mutate(mz = str_replace(mz, "Intensity_mz_", ""),
                mz = factor(mz, levels = c("678.3508", "695.3790", "700.3328")))

# Define offset step sizes
vertical_step <- 1e6        # vertical spacing (Intensity)
horizontal_step <- 0.02     # horizontal spacing (RT)

# Apply offset for 3D effect
adducts_long_offset <- adducts_long %>%
  mutate(mz_order = factor(mz, levels = c("678.3508", "700.3328", "695.3790")),
         mz_group = as.numeric(mz_order),
         Intensity_offset = Intensity + vertical_step * (mz_group - 1),
         RT_offset = `RT (min)` + horizontal_step * (mz_group - 1))

# Plot with both vertical and horizontal offset
plot_adduuct <- ggplot(adducts_long_offset, aes(x = RT_offset, y = Intensity_offset, color = mz)) +
  geom_line(size = 0.5) +
  labs(x = "Retention Time", y = "Intensity",
       title = "Adducts") +
  scale_y_continuous(labels = scales::scientific) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 7),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#ggsave(plot = plot_adduuct, filename = "plot_adducts.svg", device = "svg", dpi = "retina", width = 5, height = 3)

neg_add <- read_csv("data_metabolomics/neg_adduct.csv") %>%
  dplyr::filter(`Retention time` < 4.4 & `Retention time` > 4)

plot_neg <- ggplot(neg_add, aes(x = `Retention time`, y = `Base peak intensity`)) +
  geom_line(color = "black", size = 0.5) + scale_fill_viridis_d() +
  labs(x = "Retention Time", y = "Intensity", title = "neg") +
  theme_minimal() + scale_y_continuous(labels = scales::scientific) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0.5),  
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")  +
  theme(plot.title = element_text(size = 8), axis.title = element_text(size = 7),
        axis.text = element_text(size = 5))

#ggsave(plot = plot_neg, filename = "plot_neg_adduct.svg", device = "svg", dpi = "retina", width = 5, height = 3)
