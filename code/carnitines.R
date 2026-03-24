setwd("~/Desktop/Manuscript_AMP_Perinatal")

library(tidyverse)
library(readxl)
library(ggpubr)
library(Spectra)
library(MsBackendMgf)

# Read carnitine table
carn <- read_excel("data_metabolomics/supplemenatry_tables.xlsx", sheet = 8)

# Homologous series
carn_0 <- carn %>% 
  dplyr::filter(str_detect(pattern = ":0", Chain)) %>%
  dplyr::filter(!(str_detect(pattern = "OH|DC", Chain)))

carn_0 %>% ggscatter(x = "RT", y = "mz", add = "reg.line", label = "Name") + stat_cor()

carn_1 <- carn %>% 
  dplyr::filter(str_detect(pattern = ":1", Chain)) %>%
  dplyr::filter(!(str_detect(pattern = "OH|DC", Chain)))

carn_1 %>% ggscatter(x = "RT", y = "mz", add = "reg.line", label = "Name") + stat_cor()

carn_2 <- carn %>% 
  dplyr::filter(str_detect(pattern = ":2", Chain)) %>%
  dplyr::filter(!(str_detect(pattern = "OH|DC", Chain)))

carn_2 %>% ggscatter(x = "RT", y = "mz", add = "reg.line", label = "Name") + stat_cor()

carn_3 <- carn %>% 
  dplyr::filter(str_detect(pattern = ":3", Chain)) %>%
  dplyr::filter(!(str_detect(pattern = "OH|DC", Chain)))

carn_3 %>% ggscatter(x = "RT", y = "mz", add = "reg.line", label = "Name") + stat_cor()

carn_4 <- carn %>% 
  dplyr::filter(str_detect(pattern = ":4", Chain)) %>%
  dplyr::filter(!(str_detect(pattern = "OH|DC", Chain)))

carn_4 %>% ggscatter(x = "RT", y = "mz", add = "reg.line", label = "Name") + stat_cor()

carn_5 <- carn %>% 
  dplyr::filter(str_detect(pattern = ":5", Chain)) %>%
  dplyr::filter(!(str_detect(pattern = "OH|DC", Chain)))

carn_5 %>% ggscatter(x = "RT", y = "mz", add = "reg.line", label = "Name") + stat_cor()


# Combined plot
carn_combined <- carn %>%
  dplyr::filter(str_detect(Chain, ":[0-5]")) %>%
  dplyr::filter(!str_detect(Chain, "OH|DC")) %>%
  dplyr::mutate(Series = str_extract(Chain, ":[0-5]"))

carn_combined_plot <- carn_combined %>% 
  ggscatter(x = "RT", y = "mz", color = "Series", palette = "jco",           
            add = "reg.line", conf.int = FALSE, label = "Name", 
            repel = TRUE, legend = "right", font.label = list(size = 6)) +            
  labs(title = "Carnitine Homologous Series: m/z vs RT",
       x = "Retention Time (RT)",
       y = "m/z")+ 
  theme(plot.title = element_text(size = 8),axis.title = element_text(size = 7),
        axis.text = element_text(size = 6))

#ggsave(plot = carn_combined_plot, filename = "carnitne_series.svg", device = "svg", dpi = "retina", width = 5, height = 5)


# Generate an mgf of concordant spectra to be classified via CANOPUS
dda <- Spectra("data_metabolomics/gnps_MSV000089558.mgf", source = MsBackendMgf())
dda_ids <- data.frame(ID = dda@backend@spectraData@listData$FEATURE_ID) %>%
  dplyr::mutate(Interest = ID %in% carn$Index)
dda_filtered <- dda[dda_ids$Interest]
export(dda_filtered, MsBackendMgf(), file = "data_metabolomics/carnitine.mgf", exportTitle = FALSE)
