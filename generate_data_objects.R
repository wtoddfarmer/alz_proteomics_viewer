library(readxl)
library(dplyr)
library(tidyr)


## read in the spreadsheet containing the data form the paper and fix the names.  
johnson_xlsx <- readxl::read_xlsx("Johnson2022_41593_2021_999_MOESM3_ESM.xlsx", sheet = 3, skip = 4, trim_ws = TRUE) %>% 
  setNames(gsub("\n", "_", names(.))) %>%
  setNames(gsub(" ", "", names(.))) %>%
  setNames(gsub("\r", "", names(.))) %>%
  setNames(gsub("-", "_", names(.))) 


## read in the metadata spreadsheet containing the data form the paper and save it to a n R object.  
johnson_sample_meta <-  readxl::read_xlsx("Johnson2022_41593_2021_999_MOESM3_ESM.xlsx", sheet = 3, trim_ws = TRUE, range ="S3:SL5"  )
colnames(johnson_sample_meta) <- johnson_sample_meta[2,]
johnson_sample_meta <- johnson_sample_meta[-c(2),]
johnson_sample_meta <- johnson_sample_meta %>%
 pivot_longer(cols = everything(), names_to = "case", values_to = "diagnosis") %>% 
  mutate(diagnosis = factor(diagnosis, levels = c("Control", "AsymAD", "AD")))

saveRDS(johnson_sample_meta, file = "alz_view/sample_meta.RDS")


## Extract the meaningful summary columns from the published spreadsheet and save as a R object.  

johnson_summary <- johnson_xlsx %>%
  dplyr::select(1:6, 507:516) %>%
  select(Symbol,
         UniprotAccession,
         Description,
         ANOVA_PValue,
         HolmAdj._PValue_ADvs.Control,
         HolmAdj._PValue_AsymADvs.Control,
         HolmAdj._PValue_ADvs.AsymAD) %>% 
  rename(padj_ADvs.Control=HolmAdj._PValue_ADvs.Control,
         padj_AsymADvs.Control=HolmAdj._PValue_AsymADvs.Control,
         padj_ADvs.AsymAD=HolmAdj._PValue_ADvs.AsymAD)

saveRDS(johnson_summary, file = "alz_view/summary.RDS")


## Extract the raw data columns from the published spreadsheet and save as an R object.  
expression_data <- johnson_xlsx %>%
  dplyr::select(3:4, 19:506)
saveRDS(expression_data, file = "alz_view/expression.RDS")

