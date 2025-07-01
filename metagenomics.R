setwd("~/Desktop/Manuscript_AMP_Perinatal")

library(tidyverse)
library(phyloseq)
library(vegan)
library(caret)
library(ggpubr)
library(ALDEx2)
library(patchwork)
library(UpSetR)
library(KEGGREST)

# Function to fetch pathway information for a single KO
fetch_pathway_info <- function(ko_id) {
  pathway_info <- keggGet(ko_id)
  if (length(pathway_info) == 0 || is.null(pathway_info[[1]]$PATHWAY)) {
    return(data.frame(ENTRY = ko_id, PATHWAY_MAP = NA))
  }
  pathways <- pathway_info[[1]]$PATHWAY
  df <- data.frame(
    ENTRY = ko_id,
    PATHWAY_MAP = names(pathways),
    Description = unname(pathways),
    stringsAsFactors = FALSE
  )
  return(df)
}


# Read data
ogu_table <- read_tsv("data_metagenomics/ogu.tsv")
path_table <- read_tsv("data_metagenomics/pathway.tsv")
metadata <- read_csv("data_metagenomics/metadata_fixed.csv") %>%
  dplyr::mutate(Day = gsub("B15|B16", "B14", Day)) %>%
  dplyr::mutate(Day = gsub("B22|B23", "B21", Day))
lineage_table <- read.delim("data_metagenomics/lineages.txt", header = FALSE)

# I will work only with samples collected from the dams and collected 
# in one single batch as I have previously explored this data and found, also in the
# metabolomics data, a huge variation between batches difficult to correct. 
# I will also work only on Hypo1 and Hypo2, since Hypo3 will not be considered.

metadata_filter <- metadata %>% 
  dplyr::filter(Mother_Infant != "Infant") %>%
  dplyr::filter(Hypothesis != 3) %>%
  dplyr::filter(Batch != 2)

# Fix column names
colnames(ogu_table) <- gsub("15345.", "", colnames(ogu_table))
colnames(ogu_table) <- gsub("\\.", "_", colnames(ogu_table))
colnames(path_table) <- gsub("15345.", "", colnames(path_table))
colnames(path_table) <- gsub("\\.", "_", colnames(path_table))

# Keep samples of interest
ogu_table_filter <- ogu_table %>% column_to_rownames("OTUID") %>%
  dplyr::select(metadata_filter$SampleID)
path_table_filter <- path_table %>% column_to_rownames("OTUID") %>%
  dplyr::select(metadata_filter$SampleID)

# Fix taxonomy
lineage_table_split <- lineage_table %>% 
  separate(V2, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "; ")

lineage_table_split <- lineage_table_split %>%
  mutate(Phylum = ifelse(Phylum == "p__", paste(Kingdom, Phylum, sep = "_"), Phylum)) %>%
  mutate(Class = ifelse(Class == "c__", paste(Phylum, Class, sep = "_"), Class)) %>%
  mutate(Order = ifelse(Order == "o__", paste(Class, Order, sep = "_"), Order)) %>%
  mutate(Family = ifelse(Family == "f__", paste(Order, Family, sep = "_"), Family)) %>%
  mutate(Genus = ifelse(Genus == "g__", paste(Family, Genus, sep = "_"), Genus)) %>%
  mutate(Species = ifelse(Species == "s__", paste(Genus, Species, sep = "_"), Species)) %>%
  dplyr::filter(V1 %in% rownames(ogu_table_filter)) %>% arrange(V1) %>%
  column_to_rownames("V1") %>% as.matrix()

# Generate PS object
ps <- phyloseq(otu_table(ogu_table_filter, taxa_are_rows = TRUE), tax_table(lineage_table_split))


# Sequencing depth
sequencing_depth <- sample_sums(ps) %>% as.data.frame() %>%
  rownames_to_column("SampleID") %>% arrange(SampleID) %>%
  dplyr::mutate(Group = case_when(str_detect(pattern = "MO", SampleID) ~ "Sample", TRUE ~ "Control"))
colnames(sequencing_depth)[2] <- "Seq_depth"

seq_depth_sample <- sequencing_depth %>% dplyr::filter(Group == "Sample") %>% arrange(desc(Seq_depth))
seq_depth_sample %>% summarise(median_depth = median(Seq_depth)) # ~2,000,000

# Remove samples with less than 500,000 reads
sample_to_remove <- sequencing_depth %>% 
  dplyr::filter(Group == "Sample" & Seq_depth < 500000) %>% 
  left_join(metadata_filter)

table_species <- as.data.frame(ps@otu_table) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("SampleID") %>% 
  dplyr::filter(!(SampleID %in% sample_to_remove$SampleID))

path_table_final <- path_table_filter %>% 
  dplyr::select(-sample_to_remove$SampleID)


# Extract controls
table_ctl <- table_species %>% dplyr::filter(!(str_detect(SampleID, "MO"))) %>% 
  arrange(SampleID) %>% select_if(~ any(. != 0)) %>% column_to_rownames("SampleID")
table_ctl_sum <- table_ctl %>% rowSums() %>% as.data.frame()

table_zymo <- table_ctl %>% rownames_to_column("SampleID") %>% 
  dplyr::filter(str_detect(SampleID, "ZymoMock")) %>% 
  select_if(~ any(. != 0)) %>% column_to_rownames("SampleID")
table_zymo_sum <- table_zymo %>% colSums() %>% as.data.frame() %>% 
  rownames_to_column("OGU") %>% left_join((rownames_to_column(as.data.frame(lineage_table_split), "OGU"))) 
colnames(table_zymo_sum)[2] <- "count"

top_zymo <- table_zymo_sum %>% arrange(desc(count)) %>%
  head(500) %>% group_by(Species) %>% 
  summarise(total = sum(count)) %>% arrange(desc(total))
# top detected are correct and total ~97% accuracy in classification

blank_ctl <- table_ctl %>% rownames_to_column("SampleID") %>% 
  dplyr::filter(str_detect(SampleID, "Blank")) %>% 
  column_to_rownames("SampleID") %>% 
  dplyr::select(where(~ sum(.) > 100)) %>% 
  summarise(across(everything(), sum)) %>%
  t() %>% as.data.frame() %>% arrange(desc(V1)) %>%
  rownames_to_column("OGU") %>% 
  left_join(lineage_table_split %>% as.data.frame() %>% rownames_to_column("OGU"))
# traces of salmonella enterica are found in blanks. plate contamination
# due to the zymo_mock? Remove associated OGUs (G004127505, G000195995, G002032965)


# Check samples
species_sample <- table_species %>% dplyr::filter(str_detect(SampleID, "MO")) %>% 
  select_if(~ any(. != 0)) %>% column_to_rownames("SampleID") %>% 
  dplyr::select(-c("G004127505", "G000195995", "G002032965"))

# RCLR Transformation
species_clr <- species_sample %>% decostand(method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(species_clr %>% select_at(vars(-one_of(caret::nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_filter)

PCA_whole_scores$Day <- factor(PCA_whole_scores$Day, 
                               levels = c("E16", "E17", "E18", "B1", "B2", "B3", 
                                          "B7", "B10", "B14", "B15", "B16", "B21", 
                                          "B22", "B23"))

PCA_whoe_plots <- list()

for (i in c("Box", "Cage", "Day", "Amp", "Hypothesis")) {
  
  PCA_plot <- PCA_whole_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Microbiome"),
              xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_whoe_plots[[i]] <- PCA_plot
  
}

PCA_whoe_plots_final <- wrap_plots(PCA_whoe_plots, nrow = 3)

# There is a strong separation between the two cohorts (observed also in metabolomics)
# so work on them separately. Ampicillin effect seems evident. Two samples from
# Hypo2 are clustering with Hypo1 (MO_C1_B14_Hypo2_1, MO_C7_2_B21_Hypo2_1) but
# they appeared to be from later timepoints and not AMP treated.

# PERMANOVA
dist_metabolites <- vegdist(species_clr, method = "euclidean")
disper <- betadisper(dist_metabolites, PCA_whole_scores$Amp)
anova(disper)
permanova <- adonis2(dist_metabolites ~ Hypothesis + Amp + Day + Box, PCA_whole_scores, 
                     na.action = na.omit, by = "terms")


# Consider each hypothesis/cohort on its own given the strong separation
hypo1 <- metadata_filter %>% dplyr::filter(Hypothesis == 1) %>%
  dplyr::filter(SampleID %in% rownames(species_sample))
hypo2 <- metadata_filter %>% dplyr::filter(Hypothesis == 2) %>% 
  dplyr::filter(SampleID %in% rownames(species_sample))


##########
# HYPO 1 #
##########
hypo1_samples <- species_sample %>% rownames_to_column("SampleID") %>%
  dplyr::filter(SampleID %in% hypo1$SampleID) %>% 
  dplyr::filter(!(SampleID %in% c("MO_C1_B21_Hypo1_1", "MO_C10_B21_Hypo1_1"))) %>% #outliers
  column_to_rownames("SampleID")

hypo1_path <- path_table_final %>% t() %>% as.data.frame() %>% 
  rownames_to_column("SampleID") %>% dplyr::filter(SampleID %in% hypo1$SampleID) %>% 
  dplyr::filter(!(SampleID %in% c("MO_C1_B21_Hypo1_1", "MO_C10_B21_Hypo1_1"))) %>%
  column_to_rownames("SampleID")

# Remove bacterial species present in less than 10% of samples
prop_zeros <- colMeans(hypo1_samples == 0) %>% as.data.frame()
colnames(prop_zeros)[1] <- "Proportion"
cols_to_remove <- prop_zeros %>% dplyr::filter(Proportion > 0.9) %>% rownames_to_column("ID")
hypo1_samples_filter <- hypo1_samples %>% dplyr::select(-cols_to_remove$ID)

prop_zeros_path <- colMeans(hypo1_path == 0) %>% as.data.frame()
colnames(prop_zeros_path)[1] <- "Proportion"
cols_to_remove_path <- prop_zeros_path %>% dplyr::filter(Proportion > 0.9) %>% rownames_to_column("ID")
hypo1_path_filter <- hypo1_path %>% dplyr::select(-cols_to_remove_path$ID)

# Remove OGUs accounting for less that 0.0001% in any given sample
hypo1_ogu_ra <- hypo1_samples_filter %>% 
  dplyr::mutate(TotalRead = rowSums(.)) %>%
  dplyr::mutate(across(-TotalRead, ~ . / TotalRead)) %>%
  dplyr::select(where(~ any(. > 0.0001))) %>%
  dplyr::select(-TotalRead)

hypo1_samples_final <- hypo1_samples_filter %>% dplyr::select(colnames(hypo1_ogu_ra))

# RCLR transformation
hypo1_clr <- hypo1_samples_final %>% decostand(method = "rclr")

# PCA
PCA_hypo1 <- mixOmics::pca(hypo1_clr %>% select_at(vars(-one_of(caret::nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_hypo1_scores <- data.frame(PCA_hypo1$variates$X) %>% rownames_to_column("SampleID") %>% left_join(hypo1)
PCA_hypo1_scores$Day <- factor(PCA_hypo1_scores$Day, levels = c("E16", "E17", "E18", "B1", 
                                                                "B2", "B7", "B14", "B15", 
                                                                "B16", "B21", "B22", "B23"))

PCA_hypo1_plots <- list()

for (i in c("Box", "Cage", "Day", "Amp")) {
  
  PCA_hypo1_plot <- PCA_hypo1_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Hypo 1 _ All Samples"),
              xlab = paste("PC1 (", round(PCA_hypo1$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_hypo1$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_hypo1_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  
  PCA_hypo1_plots[[i]] <- PCA_hypo1_plot
  
}

PCA_hypo1_plots_final <- wrap_plots(PCA_hypo1_plots, nrow = 2)

i <- "Day"

PCA_plot_hypo1_time <- PCA_hypo1_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, title = "PCA - Antepartum", legend = "none",
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_hypo1_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6)) + scale_color_viridis_d()

#ggsave(filename = "FigS1A_ante.svg", plot = PCA_plot_hypo1_time, device = "svg", dpi = "retina", width = 2.5, height = 2.5)

# PERMANOVA
dist_metabolites <- vegdist(hypo1_clr, method = "euclidean")
disper <- betadisper(dist_metabolites, PCA_hypo1_scores$Amp)
anova(disper)
permanova <- adonis2(dist_metabolites ~ Amp + Day + Cage, PCA_hypo1_scores, na.action = na.omit, by = "terms")


# Look at AMP administration period only --> combine multiple timepoints given low number of samples
hypo1_amp <- hypo1_samples_final %>% rownames_to_column("SampleID") %>% 
  dplyr::filter(SampleID %in% (hypo1 %>% dplyr::filter(Day %in% c("E17", "B1", "B2")))$SampleID) %>%
  column_to_rownames("SampleID")

hypo1_path <- hypo1_path_filter %>% rownames_to_column("SampleID") %>% 
  dplyr::filter(SampleID %in% (hypo1 %>% dplyr::filter(Day %in% c("E17", "B1", "B2")))$SampleID) %>%
  column_to_rownames("SampleID")

# Remove species present in less than 10% of samples
prop_zeros <- colMeans(hypo1_amp == 0)
cols_to_remove <- names(prop_zeros[prop_zeros > 0.9]) %>% as.data.frame()
colnames(cols_to_remove) <- "OGU"
hypo1_amp_filter <- hypo1_amp %>% dplyr::select(-cols_to_remove$OGU)

# RCLR transformation
hypo1_amp_clr <- hypo1_amp_filter %>% decostand(method = "rclr")

# PCA
PCA_hypo1_amp <- mixOmics::pca(hypo1_amp_clr %>% select_at(vars(-one_of(caret::nearZeroVar(., names = TRUE)))),
                               ncomp = 2, center = TRUE, scale = TRUE)
PCA_hypo1_amp_scores <- data.frame(PCA_hypo1_amp$variates$X) %>% rownames_to_column("SampleID") %>% left_join(hypo1)

i <- "Amp"

PCA_hypo1_plot <- PCA_hypo1_amp_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, legend = "none",
            title = paste("PCA - Antepartum"), palette = c("#287DAB", "#E5BF86"), shape = "Day",
            xlab = paste("PC1 (", round(PCA_hypo1_amp$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_hypo1_amp$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_hypo1_amp_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

#ggsave(filename = "Fig2A_ante.svg", plot = PCA_hypo1_plot, device = "svg", dpi = "retina", width = 2, height = 2)

# PERMANOVA
dist_metabolites <- vegdist(hypo1_amp_clr, method = "euclidean")
disper <- betadisper(dist_metabolites, PCA_hypo1_amp_scores$Amp)
anova(disper)
permanova <- adonis2(dist_metabolites ~ Amp + Day + Cage, PCA_hypo1_amp_scores, na.action = na.omit, by = "terms")


# Alpha diversity
ps_hypo1_dam_amp <- phyloseq(otu_table(as.matrix(hypo1_samples), taxa_are_rows = FALSE), 
                             tax_table(lineage_table_split))

ps_hypo1_dam_amp_rar <- rarefy_even_depth(ps_hypo1_dam_amp, 
                                          sample.size = min(sample_sums(ps_hypo1_dam_amp)), 
                                          rngseed = 1234,
                                          replace = FALSE, 
                                          trimOTUs = TRUE)

hypo1_dam_pre_alpha <- ps_hypo1_dam_amp_rar %>%
  estimate_richness(split = TRUE, measures = c("Observed", "Shannon", "Simpson", "Fisher", "Chao1")) %>%
  mutate(SampleID = sample_names(ps_hypo1_dam_amp_rar))

hypo1_dam_pre_alpha_info <- hypo1_dam_pre_alpha %>% left_join(hypo1)

hypo1_dam_pre_shannon_plot <- hypo1_dam_pre_alpha_info %>% 
  dplyr::filter(Day %in% c("E17", "B1", "B2")) %>%
  ggboxplot(x = "Amp", y = "Shannon", add = "jitter", title = "Antepartum AMP",
            ylab = "Shannon Diversity Index", xlab = FALSE, legend = "none",
            add.params = list(color = "Amp", alpha = 0.6), 
            palette = c("#287DAB", "#E5BF86")) + ylim(0,4.5) + stat_compare_means() +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

#ggsave(filename = "Fig2B_ante.svg", plot = hypo1_dam_pre_shannon_plot, device = "svg", dpi = "retina", width = 1.5, height = 2)

# Longitudinal Shannon
hypo1_shannon_long <- hypo1_dam_pre_alpha_info %>% 
  mutate(Day = case_when(Day == "E16" ~ 1, Day == "E17" ~ 2, Day == "E18" ~ 3,
                         Day == "B1" ~ 4, Day == "B2" ~ 5, Day == "B7" ~ 10,
                         Day == "B14" ~ 17, Day == "B21" ~ 24)) %>% 
  mutate_at("Day", as.numeric) %>%
  ggscatter(x = "Day", y = "Shannon", color = "Amp", add = "loess", ylab = "Shannon Diversity Index",
            legend = "none", palette = c("#287DAB", "#E5BF86"), alpha = 0.5,
            title = "Alpha Diversity - Prepartum Exposure") + ylim(0,4) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 10, 17, 24)) + 
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Test timepoints
p_value_hypo1_alpha <- hypo1_dam_pre_alpha_info %>%
  dplyr::select(Shannon, Day, Amp) %>%
  dplyr::mutate(Day = case_when(Day == "E18" ~ "E17",
                                Day == "B1" ~ "B2",
                                TRUE ~ Day)) %>%
  group_by(Day) %>% nest() %>%
  dplyr::mutate(test = purrr::map(data, function(df) wilcox.test(Shannon ~ Amp, data = df)),
                p_value = purrr::map_dbl(test, ~.x$p.value)) %>%
  ungroup() %>% dplyr::mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
  dplyr::select(Day, p_value, p_adjust) %>% arrange(Day)


# ALDEx2 - Differential abundance 
set.seed(123)
cond_amp_hypo1 <- as.vector(PCA_hypo1_amp_scores$Amp)
ALDE_amp_hypo1 <- aldex.clr(t(hypo1_amp_filter), cond_amp_hypo1, denom = "all") 
ALDE_amp_hypo1_test <- aldex.ttest(ALDE_amp_hypo1) 
ALDE_amp_hypo1_eff <- aldex.effect(ALDE_amp_hypo1, useMC = TRUE)

# Store results
res_whole_amp_hypo1 <- data.frame(rownames(ALDE_amp_hypo1_eff), ALDE_amp_hypo1_eff, ALDE_amp_hypo1_test)
res_whole_amp_hypo1_filter <- res_whole_amp_hypo1 %>% dplyr::filter(wi.eBH < 0.05)
colnames(res_whole_amp_hypo1_filter)[1] <- "OGU"

res_whole_amp_hypo1_filter_species <- res_whole_amp_hypo1_filter %>% 
  left_join(lineage_table_split %>% as.data.frame() %>% rownames_to_column("OGU")) %>%
  arrange(effect)

res_whole_amp_hypo1_filter_species %>% dplyr::filter(effect > 0) %>%
  group_by(Genus) %>% summarise(count = n()) %>% arrange(desc(count))

#write_csv(x = res_whole_amp_hypo1_filter_species, file = "Table_S1.csv")

# Make summary figure
genera_amp_hypo1 <- res_whole_amp_hypo1_filter_species %>% 
  dplyr::mutate(Presence = case_when(effect >= 0 ~ "Amp_Yes", TRUE ~ "Amp_No")) %>%
  group_by(Genus) %>% count(Presence) %>% arrange(desc(n)) %>% 
  dplyr::filter(n > 0) %>% dplyr::filter(Genus != "g__")

hypo1_pbs_top <- res_whole_amp_hypo1_filter_species %>% 
  arrange(effect) %>% dplyr::filter(effect < 0) %>% head(n = 15) %>%
  dplyr::select(effect, Species) %>% dplyr::mutate(Species = gsub("s__", "", Species))
hypo1_amp_top <- res_whole_amp_hypo1_filter_species %>% 
  arrange(desc(effect)) %>% dplyr::filter(effect > 0) %>% head(n = 10) %>%
  dplyr::select(effect, Species) %>% dplyr::mutate(Species = gsub("s__", "", Species))

hypo1_top <- rbind(hypo1_pbs_top, hypo1_amp_top) %>% arrange(effect)
hypo1_top$Species <- factor(hypo1_top$Species, levels = unique(hypo1_top$Species[order(hypo1_top$effect)]))

# Create the lollipop plot
hypo1_diff <- ggplot(hypo1_top, aes(x = Species, y = effect)) +
  geom_segment(aes(xend = Species, y = 0, yend = effect), color = "gray") +
  geom_point(aes(color = effect), size = 2) +
  scale_color_gradient2(low = "#287DAB", mid = "white", high = "#E5BF86", midpoint = 0) +
  coord_flip() + labs(title = "Antepartum", x = "Species", y = "Effect") +
  theme_minimal() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6), legend.position = "none")

#ggsave(plot = hypo1_diff, filename = "Fig2C_ante.svg", device = "svg", height = 2.8, width = 2.2)

# Plot ratio
ogu_amp_hypo1_yes_hypo1 <- res_whole_amp_hypo1_filter %>% 
  dplyr::filter(effect > 0) %>% arrange(desc(effect))
ogu_amp_hypo1_no_hypo1 <- res_whole_amp_hypo1_filter %>% 
  dplyr::filter(effect < 0) %>% arrange(effect)

ratio_hypo1 <- hypo1_samples_filter %>% dplyr::select(res_whole_amp_hypo1_filter$OGU) %>%
  dplyr::mutate(AMP_Yes = rowSums(dplyr::select(., ogu_amp_hypo1_yes_hypo1$OGU))) %>%
  dplyr::mutate(AMP_No = rowSums(dplyr::select(., ogu_amp_hypo1_no_hypo1$OGU))) %>%
  dplyr::mutate(Ratio = log(AMP_Yes/AMP_No)) %>% dplyr::filter(!(is.infinite(Ratio))) %>%
  rownames_to_column("SampleID") %>% left_join(hypo1)

plot_hypo1_ratio <- ratio_hypo1 %>% 
  mutate(Day = case_when(Day == "E16" ~ 1, Day == "E17" ~ 2, Day == "E18" ~ 3,
                         Day == "B1" ~ 4, Day == "B2" ~ 5, Day == "B7" ~ 10,
                         Day == "B14" ~ 17, Day == "B21" ~ 24)) %>% 
  mutate_at("Day", as.numeric) %>%
  ggscatter(x = "Day", y = "Ratio", color = "Amp", add = "loess", ylab = "Ln(AMP/PBS)",
            legend = "none", palette = c("#287DAB", "#E5BF86"), alpha = 0.5,
            title = "Prepartum Exposure") + 
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 10, 17, 24)) + 
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Test each timepoint 
p_value_hypo1 <- ratio_hypo1 %>%
  dplyr::select(Ratio, Day, Amp) %>%
  dplyr::mutate(Day = case_when(Day == "E18" ~ "E17",
                                Day == "B1" ~ "B2",
                                TRUE ~ Day)) %>%
  group_by(Day) %>% nest() %>%
  dplyr::mutate(test = purrr::map(data, function(df) wilcox.test(Ratio ~ Amp, data = df)),
                p_value = purrr::map_dbl(test, ~.x$p.value)) %>%
  ungroup() %>% dplyr::mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
  dplyr::select(Day, p_value, p_adjust) %>% arrange(Day)

# ALDEx2 - Pathway
all(PCA_hypo1_amp_scores$SampleID == rownames(hypo1_path))
cond_amp_hypo1 <- as.vector(PCA_hypo1_amp_scores$Amp)
ALDE_amp_hypo1_path <- aldex.clr(t(hypo1_path), cond_amp_hypo1, denom = "all") 
ALDE_amp_hypo1_test_path <- aldex.ttest(ALDE_amp_hypo1_path) 
ALDE_amp_hypo1_eff_path <- aldex.effect(ALDE_amp_hypo1_path, useMC = TRUE)

# Store results
res_whole_amp_hypo1_path <- data.frame(rownames(ALDE_amp_hypo1_eff_path), ALDE_amp_hypo1_eff_path, ALDE_amp_hypo1_test_path)
res_whole_amp_hypo1_path_filter <- res_whole_amp_hypo1_path %>% dplyr::filter(wi.eBH < 0.05) %>%
  dplyr::arrange(desc(effect))
colnames(res_whole_amp_hypo1_path_filter)[1] <- "Pathway"

pathway_hypo1_list <- lapply(res_whole_amp_hypo1_path_filter$Pathway, fetch_pathway_info)
pathway_hypo1_combined <- do.call(rbind, pathway_hypo1_list)
pathway_hypo1_info <- pathway_hypo1_combined[, c("ENTRY", "Description")]

res_whole_amp_hypo1_path_final <- res_whole_amp_hypo1_path_filter %>%
  left_join(pathway_hypo1_info, by = c("Pathway" = "ENTRY"))

#write_csv(x = res_whole_amp_hypo1_path_final, file = "hypo1_differential_pathways.csv")


##########
# HYPO 2 #
##########
hypo2_samples <- species_sample %>% rownames_to_column("SampleID") %>%
  dplyr::filter(SampleID %in% hypo2$SampleID) %>% 
  column_to_rownames("SampleID")

hypo2_path <- path_table_final %>% t() %>% as.data.frame() %>% 
  rownames_to_column("SampleID") %>% dplyr::filter(SampleID %in% hypo2$SampleID) %>%
  column_to_rownames("SampleID")

# Remove species present in less than 10% of samples
prop_zeros <- colMeans(hypo2_samples == 0) %>% as.data.frame()
colnames(prop_zeros)[1] <- "Proportion"
cols_to_remove <- prop_zeros %>% dplyr::filter(Proportion > 0.9) %>% rownames_to_column("ID")
hypo2_samples_filter <- hypo2_samples %>% dplyr::select(-cols_to_remove$ID)

prop_zeros_path <- colMeans(hypo2_path == 0) %>% as.data.frame()
colnames(prop_zeros_path)[1] <- "Proportion"
cols_to_remove_path <- prop_zeros_path %>% dplyr::filter(Proportion > 0.9) %>% rownames_to_column("ID")
hypo2_path_filter <- hypo2_path %>% dplyr::select(-cols_to_remove_path$ID)

# Remove OGUs accounting for less that 0.0001% in any given sample
hypo2_ogu_ra <- hypo2_samples_filter %>% 
  dplyr::mutate(TotalRead = rowSums(.)) %>%
  dplyr::mutate(across(-TotalRead, ~ . / TotalRead)) %>%
  dplyr::select(where(~ any(. > 0.0001))) %>%
  dplyr::select(-TotalRead)

hypo2_samples_final <- hypo2_samples_filter %>% dplyr::select(colnames(hypo2_ogu_ra))

# RCLR transformation
hypo2_clr <- hypo2_samples_final %>% decostand(method = "rclr")

# PCA
PCA_hypo2 <- mixOmics::pca(hypo2_clr %>% select_at(vars(-one_of(caret::nearZeroVar(., names = TRUE)))),
                           ncomp = 2, center = TRUE, scale = TRUE)
PCA_hypo2_scores <- data.frame(PCA_hypo2$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(hypo2)

PCA_hypo2_scores$Day <- factor(PCA_hypo2_scores$Day, levels = c("B1", "B2", "B3", "B7", "B14", "B21"))

PCA_hypo2_plots <- list()

for (i in c("Cage", "Day", "Amp")) {
  
  PCA_hypo2_plot <- PCA_hypo2_scores %>%
    ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
              title = paste("PCA - Hypo 2 _ All Samples"),
              xlab = paste("PC1 (", round(PCA_hypo2$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
              ylab = paste("PC2 (", round(PCA_hypo2$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
              ggtheme = theme_classic()) +
    geom_point(data = PCA_hypo2_scores %>% group_by((!!sym(i))) %>% 
                 summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
    theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
          axis.text = element_text(size = 4)) + coord_fixed()
  
  PCA_hypo2_plots[[i]] <- PCA_hypo2_plot
  
}

PCA_hypo2_plots_final <- wrap_plots(PCA_hypo2_plots, nrow = 2)
# Antibiotic effect is evident

i <- "Day"

PCA_plot_hypo2_time <- PCA_hypo2_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, title = "PCA - Postpartum", legend = "none",
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_hypo2_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6)) + scale_color_viridis_d()

#ggsave(filename = "FigS1A_post.svg", plot = PCA_plot_hypo2_time, device = "svg", dpi = "retina", width = 2.5, height = 2.5)

# PERMANOVA
dist_metabolites <- vegdist(hypo2_clr, method = "euclidean")
disper <- betadisper(dist_metabolites, PCA_hypo2_scores$Amp)
anova(disper)
permanova <- adonis2(dist_metabolites ~ Amp + Day + Cage, PCA_hypo2_scores, na.action = na.omit, by = "terms")


# Look at administration period only for Amp (combined multiple timepoints given low n)
hypo2_amp <- hypo2_samples_final %>% rownames_to_column("SampleID") %>% 
  dplyr::filter(SampleID %in% (hypo2 %>% dplyr::filter(Day %in% c("B2", "B3")))$SampleID) %>%
  column_to_rownames("SampleID")

hypo2_path <- hypo2_path_filter %>% rownames_to_column("SampleID") %>% 
  dplyr::filter(SampleID %in% (hypo2 %>% dplyr::filter(Day %in% c("B2", "B3")))$SampleID) %>%
  column_to_rownames("SampleID")

# Remove species present in less than 5% of samples
prop_zeros <- colMeans(hypo2_amp == 0)
cols_to_remove <- names(prop_zeros[prop_zeros > 0.9]) %>% as.data.frame()
colnames(cols_to_remove) <- "OGU"
hypo2_amp_filter <- hypo2_amp %>% dplyr::select(-cols_to_remove$OGU)

# RCLR transformation
hypo2_amp_clr <- hypo2_amp_filter %>% decostand(method = "rclr")

# PCA
PCA_hypo2_amp <- mixOmics::pca(hypo2_amp_clr, ncomp = 2, center = TRUE, scale = TRUE)
PCA_hypo2_amp_scores <- data.frame(PCA_hypo2_amp$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(hypo2)

i <- "Amp"

PCA_hypo2_plot <- PCA_hypo2_amp_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, legend = "none",
            title = paste("PCA - Postpartum"), palette = c("#287DAB", "#E5BF86"), shape = "Day",
            xlab = paste("PC1 (", round(PCA_hypo2_amp$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_hypo2_amp$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_hypo2_amp_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Sample MO_9_2_B3_Hypo2_1 must be FALSE. error in the metadata. Fix

PCA_hypo2_amp_scores <- PCA_hypo2_amp_scores %>% 
  dplyr::mutate(Amp = case_when(SampleID == "MO_9_2_B3_Hypo2_1" ~ "FALSE",
                                TRUE ~ Amp))

PCA_hypo2_plot <- PCA_hypo2_amp_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6, legend = "none",
            title = paste("PCA - Postpartum"), palette = c("#287DAB", "#E5BF86"), shape = "Day",
            xlab = paste("PC1 (", round(PCA_hypo2_amp$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_hypo2_amp$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_hypo2_amp_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Fix also hypo2 and metadata_filter
hypo2 <- hypo2 %>% dplyr::mutate(Amp = case_when(SampleID == "MO_9_2_B3_Hypo2_1" ~ "FALSE",
                                                 TRUE ~ Amp))
metadata_filter <- metadata_filter %>%  dplyr::mutate(Amp = case_when(SampleID == "MO_9_2_B3_Hypo2_1" ~ "FALSE",
                                                                      TRUE ~ Amp))

#ggsave(filename = "Fig2A_post.svg", plot = PCA_hypo2_plot, device = "svg", dpi = "retina", width = 2, height = 2)

# PERMANOVA
dist_metabolites <- vegdist(hypo2_amp_clr, method = "euclidean")
disper <- betadisper(dist_metabolites, PCA_hypo2_amp_scores$Amp)
anova(disper)
permanova <- adonis2(dist_metabolites ~ Amp + Day + Cage, PCA_hypo2_amp_scores, na.action = na.omit, by = "terms")


# Alpha diversity
ps_hypo2_dam_amp <- phyloseq(otu_table(as.matrix(hypo2_samples), taxa_are_rows = FALSE), 
                             tax_table(lineage_table_split))

ps_hypo2_dam_amp_rar <- rarefy_even_depth(ps_hypo2_dam_amp, 
                                          sample.size = min(sample_sums(ps_hypo2_dam_amp)), 
                                          rngseed = 1234,
                                          replace = FALSE, 
                                          trimOTUs = TRUE)

hypo2_dam_pre_alpha <- ps_hypo2_dam_amp_rar %>%
  estimate_richness(split = TRUE, measures = c("Observed", "Shannon", "Simpson", "Fisher", "Chao1")) %>%
  mutate(SampleID = sample_names(ps_hypo2_dam_amp_rar))

hypo2_dam_pre_alpha_info <- hypo2_dam_pre_alpha %>% left_join(hypo2) %>%
  dplyr::mutate(Amp = factor(Amp, levels = c("FALSE", "TRUE")))

hypo2_dam_pre_shannon_plot <- hypo2_dam_pre_alpha_info %>% 
  dplyr::filter(Day %in% c("B2", "B3")) %>%
  ggboxplot(x = "Amp", y = "Shannon", add = "jitter", title = "Alpha Diversity",
            ylab = "Shannon Diversity Index", xlab = FALSE, legend = "none",
            add.params = list(color = "Amp", alpha = 0.6), 
            palette = c("#287DAB", "#E5BF86")) + ylim(0,4.5) +
  stat_compare_means() +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

#ggsave(filename = "Fig2B_post.svg", plot = hypo2_dam_pre_shannon_plot, device = "svg", dpi = "retina", width = 1.5, height = 2)

# Longitudinal Shannon
hypo2_shannon_long <- hypo2_dam_pre_alpha_info %>%
  mutate(Day = case_when(Day == "B1" ~ 4, Day == "B2" ~ 5, Day == "B3" ~ 6, Day == "B7" ~ 10,
                         Day == "B14" ~ 17, Day == "B21" ~ 24)) %>% 
  mutate_at("Day", as.numeric) %>%
  ggscatter(x = "Day", y = "Shannon", color = "Amp", add = "loess", ylab = "Shannon Diversity Index",
            legend = "none", palette = c("#287DAB", "#E5BF86"), alpha = 0.5,
            title = "Alpha - Postpartum Exposure") + ylim(0,4) +
  scale_x_continuous(breaks = c(1, 4, 5, 6, 10, 17, 24), limits = c(1, 24)) + 
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

# Test each timepoint
p_value_hypo2_alpha <- hypo2_dam_pre_alpha_info %>%
  dplyr::select(Shannon, Day, Amp) %>%
  group_by(Day) %>% nest() %>%
  dplyr::mutate(t_test = purrr::map(data, function(df) t.test(Shannon ~ Amp, data = df)),
                p_value = purrr::map_dbl(t_test, ~.x$p.value)) %>%
  ungroup() %>% dplyr::mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
  dplyr::select(Day, p_value, p_adjust) %>%
  arrange(Day)

# Combine alpha longitudinal plots
alpha_long <- wrap_plots(hypo1_shannon_long, hypo2_shannon_long, nrow = 2)

#ggsave(plot = alpha_long, filename = "FigS1B.svg", device = "svg", dpi = "retina", width = 3, height = 4)


# ALDEx2 - Differential abundance
set.seed(123)
cond_amp_hypo2 <- as.vector(PCA_hypo2_amp_scores$Amp)
ALDE_amp_hypo2 <- aldex.clr(t(hypo2_amp_filter), cond_amp_hypo2, denom = "all") 
ALDE_amp_hypo2_test <- aldex.ttest(ALDE_amp_hypo2) 
ALDE_amp_hypo2_eff <- aldex.effect(ALDE_amp_hypo2, useMC = TRUE)

# Store results
res_whole_amp_hypo2 <- data.frame(rownames(ALDE_amp_hypo2_eff), ALDE_amp_hypo2_eff, ALDE_amp_hypo2_test)
res_whole_amp_hypo2_filter <- res_whole_amp_hypo2 %>% dplyr::filter(wi.eBH < 0.05)
colnames(res_whole_amp_hypo2_filter)[1] <- "OGU"

res_whole_amp_hypo2_filter_species <- res_whole_amp_hypo2_filter %>% 
  left_join(lineage_table_split %>% as.data.frame() %>% rownames_to_column("OGU")) %>%
  arrange(effect)

res_whole_amp_hypo2_filter_species %>% dplyr::filter(effect > 0) %>%
  group_by(Genus) %>% summarise(count = n()) %>% arrange(desc(count))

#write_csv(x = res_whole_amp_hypo2_filter_species, file = "TableS2.csv")

# Make a summary figure
genera_amp_hypo2 <- res_whole_amp_hypo2_filter_species %>% 
  dplyr::mutate(Presence = case_when(effect >= 0 ~ "Amp_Yes", TRUE ~ "Amp_No")) %>%
  group_by(Genus) %>% count(Presence) %>% arrange(desc(n)) %>% 
  dplyr::filter(n > 0) %>% dplyr::filter(Genus != "g__")

hypo2_pbs_top <- head(res_whole_amp_hypo2_filter_species %>% arrange(effect), n = 15) %>%
  dplyr::select(effect, Species) %>% dplyr::mutate(Species = gsub("s__", "", Species))
hypo2_amp_top <- head(res_whole_amp_hypo2_filter_species %>% arrange(desc(effect)), n = 10) %>%
  dplyr::select(effect, Species) %>% dplyr::mutate(Species = gsub("s__", "", Species))

hypo2_top <- rbind(hypo2_pbs_top, hypo2_amp_top) %>% arrange(effect)
hypo2_top$Species <- factor(hypo2_top$Species, levels = unique(hypo2_top$Species[order(hypo2_top$effect)]))

# Create the lollipop plot
hypo2_diff <- ggplot(hypo2_top, aes(x = Species, y = effect)) +
  geom_segment(aes(xend = Species, y = 0, yend = effect), color = "gray") +
  geom_point(aes(color = effect), size = 2) +
  scale_color_gradient2(low = "#287DAB", mid = "white", high = "#E5BF86", midpoint = 0) +
  coord_flip() + labs(title = "Postpartum", x = "Species", y = "Effect") +
  theme_minimal() +
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6), legend.position = "none")

#ggsave(plot = hypo2_diff, filename = "Fig2C_post.svg", device = "svg", height = 2.8, width = 2.2)

# Plot ratio
ogu_amp_hypo2_yes_hypo2 <- res_whole_amp_hypo2_filter %>% dplyr::filter(effect > 0)
ogu_amp_hypo2_no_hypo2 <- res_whole_amp_hypo2_filter %>% dplyr::filter(effect < 0)

ratio_hypo2 <- hypo2_samples_filter %>% dplyr::select(res_whole_amp_hypo2_filter$OGU) %>%
  dplyr::mutate(AMP_Yes = rowSums(dplyr::select(., ogu_amp_hypo2_yes_hypo2$OGU))) %>%
  dplyr::mutate(AMP_No = rowSums(dplyr::select(., ogu_amp_hypo2_no_hypo2$OGU))) %>%
  dplyr::mutate(Ratio = log(AMP_Yes/AMP_No)) %>%
  rownames_to_column("SampleID") %>% left_join(hypo2)

plot_hypo2_ratio <- ratio_hypo2  %>%
  mutate(Day = case_when(Day == "B1" ~ 4, Day == "B2" ~ 5, Day == "B3" ~ 6, Day == "B7" ~ 10,
                         Day == "B14" ~ 17, Day == "B21" ~ 24)) %>% 
  mutate_at("Day", as.numeric) %>%
  ggscatter(x = "Day", y = "Ratio", color = "Amp", add = "loess", ylab = "Ln(AMP/PBS)",
            legend = "none", palette = c("#287DAB", "#E5BF86"), alpha = 0.5,
            title = "Postpartum Exposure") +
  scale_x_continuous(breaks = c(1, 4, 5, 6, 10, 17, 24), limits = c(1, 24)) + 
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

p_value_hypo2 <- ratio_hypo2 %>%
  dplyr::select(Ratio, Day, Amp) %>%
  group_by(Day) %>% nest() %>%
  dplyr::mutate(test = purrr::map(data, function(df) t.test(Ratio ~ Amp, data = df)),
                p_value = purrr::map_dbl(test, ~.x$p.value)) %>%
  ungroup() %>% dplyr::mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
  dplyr::select(Day, p_value, p_adjust) %>% arrange(Day)

# Combined ratio plots
combined_ratio_time <- wrap_plots(plot_hypo1_ratio, plot_hypo2_ratio, nrow = 2)

#ggsave(plot = combined_ratio_time, filename = "Fig2E.svg", device = "svg", dpi = "retina", width = 3, height = 4)

# ALDEx2 - Pathway
all(PCA_hypo2_amp_scores$SampleID == rownames(hypo2_path))
cond_amp_hypo2 <- as.vector(PCA_hypo2_amp_scores$Amp)
ALDE_amp_hypo2_path <- aldex.clr(t(hypo2_path), cond_amp_hypo2, denom = "all") 
ALDE_amp_hypo2_test_path <- aldex.ttest(ALDE_amp_hypo2_path) 
ALDE_amp_hypo2_eff_path <- aldex.effect(ALDE_amp_hypo2_path, useMC = TRUE)

# Store results
res_whole_amp_hypo2_path <- data.frame(rownames(ALDE_amp_hypo2_eff_path), ALDE_amp_hypo2_eff_path, ALDE_amp_hypo2_test_path)
res_whole_amp_hypo2_path_filter <- res_whole_amp_hypo2_path %>% dplyr::filter(wi.eBH < 0.05) %>%
  dplyr::arrange(desc(effect))
colnames(res_whole_amp_hypo2_path_filter)[1] <- "Pathway"

pathway_hypo2_list <- lapply(res_whole_amp_hypo2_path_filter$Pathway, fetch_pathway_info)
pathway_hypo2_combined <- do.call(rbind, pathway_hypo2_list)
pathway_hypo2_info <- pathway_hypo2_combined[, c("ENTRY", "Description")]

res_whole_amp_hypo2_path_final <- res_whole_amp_hypo2_path_filter %>%
  left_join(pathway_hypo2_info, by = c("Pathway" = "ENTRY"))

#write_csv(x = res_whole_amp_hypo2_path_final, file = "hypo2_differential_pathways.csv")


##########################
# Check overlapping OGUs #
##########################

importnat_ogus <- res_whole_amp_hypo1_filter_species %>% 
  full_join(res_whole_amp_hypo2_filter_species, by = "OGU")

# Extract overlap
hypo1_ogu_overlap <- res_whole_amp_hypo1_filter_species %>% 
  dplyr::select(OGU, effect) %>% dplyr::filter(OGU %in% importnat_ogus$OGU)
hypo2_ogu_overlap <- res_whole_amp_hypo2_filter_species %>% 
  dplyr::select(OGU, effect) %>% dplyr::filter(OGU %in% importnat_ogus$OGU)

# Generate Upset plot
list_overlap <- list(
  `Prepartum PBS` = (res_whole_amp_hypo1_filter_species %>% dplyr::filter(effect < 0))$OGU,
  `Prepartum AMP` =  (res_whole_amp_hypo1_filter_species %>% dplyr::filter(effect > 0))$OGU,
  `Postpartum PBS` = (res_whole_amp_hypo2_filter_species %>% dplyr::filter(effect < 0))$OGU,
  `Postpartum AMP` = (res_whole_amp_hypo2_filter_species %>% dplyr::filter(effect > 0))$OGU)

importnat_upset <- UpSetR::upset(fromList(list_overlap), nsets = 4, nintersects = 8,
                                 point.size = 1.5, line.size = 1, text.scale = 1, keep.order = TRUE,  
                                 sets = c("Prepartum PBS", "Prepartum AMP", "Postpartum PBS", "Postpartum AMP"),
                                 queries = list(list(query = intersects, params = list("Prepartum PBS", 
                                                                                       "Postpartum PBS"), color = "#287DAB", active = T), 
                                                list(query = intersects, params = list("Prepartum AMP", 
                                                                                       "Postpartum AMP"),color = "#E5BF86", active = T)))


# Check overlapping OGUs and directionality
overlap_ogus <- res_whole_amp_hypo1_filter_species %>% 
  inner_join(res_whole_amp_hypo2_filter_species, by = "OGU") %>%
  dplyr::filter((effect.x > 0 & effect.y > 0) | (effect.x < 0 & effect.y < 0))  %>%
  dplyr::select(OGU, effect.x, effect.y, Species.x) %>%
  dplyr::mutate(effect_mean = rowMeans(cbind(effect.x, effect.y))) %>%
  arrange(effect_mean)

#write_csv(x = overlap_ogus, file = "Table_S3.csv")

# Plot ratio over time of overlapping features
overlap_amp_yes <- overlap_ogus %>% dplyr::filter(effect_mean > 0)
overlap_amp_no <- overlap_ogus %>% dplyr::filter(effect_mean < 0)

hypo1_overlap <- hypo1_samples_filter %>% dplyr::select(overlap_ogus$OGU)
hypo2_overlap <- hypo2_samples_filter %>% dplyr::select(overlap_ogus$OGU)

hypo_combined <- rbind(hypo1_overlap, hypo2_overlap) %>%
  dplyr::mutate(AMP_Yes = rowSums(dplyr::select(., overlap_amp_yes$OGU))) %>%
  dplyr::mutate(AMP_No = rowSums(dplyr::select(., overlap_amp_no$OGU))) %>%
  dplyr::mutate(Ratio = log(AMP_Yes/AMP_No)) %>%
  rownames_to_column("SampleID") %>% left_join(metadata_filter) %>%
  dplyr::select(SampleID, Ratio, Cage, Day, Amp, Hypothesis) %>%
  mutate(Day_num = case_when(Day == "E16" ~ 1, Day == "E17" ~ 2, Day == "E18" ~ 3,
                             Day == "B1" ~ 4, Day == "B2" ~ 5, Day == "B3" ~ 6, 
                             Day == "B7" ~ 10, Day == "B14" ~ 17, Day == "B21" ~ 24)) %>% 
  mutate_at("Day_num", as.numeric) %>%
  dplyr::mutate(Condition = paste(Amp, Hypothesis, sep = "_"))

mean_diff <- mean((hypo_combined %>% dplyr::filter(Hypothesis == 1 & Day == "E16"))$Ratio) / mean((hypo_combined %>% dplyr::filter(Hypothesis == 2 & Day == "B1"))$Ratio)

hypo_combined_plot <- hypo_combined %>%
  dplyr::mutate(Rato_norm = case_when(Hypothesis == 2 ~ Ratio * mean_diff, TRUE ~ Ratio)) %>%
  ggscatter(x = "Day_num", y = "Rato_norm", color = "Condition", add = "loess",
            title = "Ratio of Shared Driving Species", alpha = 0.3, legend = "none",
            xlab = "Day", ylab = "Ratio(TopAMP/TopPBS)", 
            palette = c("#009dc2","#22678d", "#e6ce8c","#ba992a")) + 
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 10, 17, 24), limits = c(1, 25)) + 
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

#ggsave(plot = hypo_combined_plot, filename = "Fig2F.svg", device = "svg", dpi = "retina", width = 2.5, height = 2)

# Test each timepoint
p_value_combined <- hypo_combined %>%
  dplyr::filter(!is.infinite(Ratio)) %>%
  dplyr::mutate(Day = case_when(Hypothesis == "1" & Day  == "E18" ~ "E17",
                                Hypothesis == "1" & Day  == "B1" ~ "B2",
                                TRUE ~ Day)) %>%
  group_by(Hypothesis, Day) %>% nest() %>%
  dplyr::mutate(t_test = purrr::map(data, function(df) t.test(Ratio ~ Amp, data = df)),
                p_value = purrr::map_dbl(t_test, ~.x$p.value)) %>%
  ungroup() %>%
  dplyr::mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
  dplyr::select(Hypothesis, Day, p_value, p_adjust) %>%
  arrange(Day)

# Check overlapping pathways
overlap_path <- res_whole_amp_hypo1_path_final %>% 
  inner_join(res_whole_amp_hypo2_path_final, by = "Pathway") %>%
  dplyr::mutate(effect_mean = rowMeans(cbind(effect.x, effect.y))) %>%
  arrange(desc(effect_mean)) %>%
  dplyr::select(Pathway, effect.x, effect.y, Description.x, effect_mean)

#write_csv(x = overlap_path, file = "TableS4.csv")
