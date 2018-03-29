# Read using read_cell_seg_data
source("Functions.R")

# Only needs to be ran once 
# Preprocess
## Cell_seg_data
folderoffolders <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports")
targetdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Cell_seg_only")
origindir <- c("/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Extra/Nahla_Analysis")

options(stringsAsFactors = F)
files <- list.files(targetdir, pattern = ".txt$")
setwd(targetdir)
lists <- lapply(files, read_cell_seg_data)
combined_df_ <- rbindlist(lists, fill = T)

df.combined <- combined_df_ %>% 
  dplyr:: select(-contains("TMA")) %>% 
  dplyr::select(-matches("Lab.ID")) %>% 
  dplyr::select(-matches("inForm.2.2.6004.14667")) %>%
  dplyr:: select(-matches("Path")) %>%
  dplyr::select(-matches("Distance.from.Process.Region.Edge..pixels.")) %>%
  dplyr::select(-matches("Process.Region.ID"))%>%
  dplyr::select(-matches("Category.Region.ID")) %>%
  dplyr::select(-matches("Total.Cells")) %>% 
  dplyr:: select(-contains("Autofluorescence")) %>% 
  dplyr::select(-contains("Axis")) %>% 
  dplyr::select(-contains("..percent")) %>% 
  dplyr:: select(-contains("Compactness")) %>% 
  dplyr::select(-contains("Nuclei"))

df1 <- df.combined 

fix_tissue_cat1 <- function(df){
  df$`Tissue Category` %<>%
    gsub("TUMOUR", "Tumour", .) %>%
    gsub("tumour", "Tumour", .) %>%
    gsub("tumor", "Tumour", .)   %>% 
    gsub("STROMA", "Stroma", .) %>%
    gsub("stroma", "Stroma", .) %>%
    gsub("stromal", "Stroma", .) %>%
    gsub("Stromal", "Stroma", .) %>%
    gsub("background", "Background", .) %>%
    gsub("BACKGROUND", "Background", .)  %>%
    gsub("Background1", "Background", .)
  return(df$`Tissue Category`)
}


df1$`Tissue Category` <- fix_tissue_cat1(df1)
df1$Phenotype <- as.character(df1$Phenotype)
df1$Phenotype[df1$Phenotype == ""] <- "DAPI"
df1$Phenotype <- fix_pheno(df1)
df2 <- filter(df1, `Tissue Category` != "Background")

df2$`Sample Name` <- as.factor(df2$`Sample Name`)

# Fix the raw files (for Radius)
txt_list <- list()
c <- 1
for(i in levels(df2$`Sample Name`)){
  name <- basename(i)
  cat('Processing', i, '\n')
  txt_list[[i]] <- droplevels(subset(df2, `Sample Name` == i))
  c <- c + 1}

for(i in names(txt_list)){
  setwd("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Export_Fixed/Cell_seg_only")
  print(i)
  write.table(txt_list[[i]], paste0(i,".txt"), quote = F, sep = "\t", row.names = F)
}
