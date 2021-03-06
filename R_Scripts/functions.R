# FUNCTIONS
## Packages
library(devtools)
library(tidyverse)
# devtools::install_github("PerkinElmer/phenoptr", build_vignettes = T)
required <- c("tidyverse",
              # "fields",
              "reshape2",
              "ggpubr",
              "phenoptr",
              "ggbiplot",
              "MASS",
              "corrplot",
              "pscl",
              "factoextra",
              "randomForest",
              "devtools",
              "data.table",
              "magrittr")

for (lib in required)
{
  if (!require(lib, character.only = T))
  {
    install.packages(lib)
    suppressMessages(library(lib, character.only = T, quietly = T))
  }
}

# Comparisons and Colours
"%!in%" <- function(x,y)!('%in%'(x,y))


cbcols <- c("Type1" = "#999999",
            "Type2" = "#56B4E9",
            "Type3" = "#E69F00",
            "Type4" = "#009E73")

# Preprocess
## Should set directories first
set_directories <- function(folderoffolders, targetdir, targetdir1, origindir){
  print("Setting up Directories")
  setwd(folderoffolders)
  thousand.folders <- list.dirs(full.names = T)
  print("Processing Cell seg data")
  filelist1 <- sapply(thousand.folders[-1], function(x)
  {list.files(x, pattern = "]_cell_seg_data.txt$", full.names = T)
  })
  filelist = unlist(filelist1)
  # Create a collated cell_seg_only folder
  subDir <- "Cell_seg_only"
  if (!file.exists(subDir)){
    dir.create(file.path(folderoffolders, subDir))
  }
  # Copy from the origin, to the target
  print("Copying Files to New Directory")
  file.copy(from = filelist, to = targetdir, recursive = F, 
            copy.mode = T)
  # New directory for CSVs
  setwd(targetdir)
  subDir1 <- "CSV"
  if (!file.exists(subDir1)){
    dir.create(file.path(targetdir, subDir1))
  }
  
  print("Processing Cell seg summary data")
  filelist2 <- sapply(thousand.folders[-1], function(x)
  {list.files(x, pattern = "]_cell_seg_summary_data.txt$", full.names = T)
  })
  filelist2a = unlist(filelist2)
  # Create a collated cell_seg_only folder
  subDir1a <- "Cell_seg_summary_only"
  if (!file.exists(subDir1)){
    dir.create(file.path(folderoffolders, subDir1a))
  }
  
  # Copy from the origin, to the target
  print("Copying Files to New Directory")
  file.copy(from = filelist2a, to = targetdir1, recursive = F, 
            copy.mode = T)
  # New directory for CSVs
  setwd(targetdir1)
  subDir1 <- "CSV"
  if (!file.exists(subDir1)){
    dir.create(file.path(targetdir, subDir1))
  }
}

## Cell_seg_data
preprocess_cell_seg <- function(folderoffolders, targetdir, subTdir, origindir){
  print("Setting up Directories")
  setwd(folderoffolders)
  thousand.folders <- list.dirs(full.names = T)
  filelist1 <- sapply(thousand.folders[-1], function(x)
  {list.files(x, pattern = "]_cell_seg_data.txt$", full.names = T)
  })
  filelist = unlist(filelist1)
  # Create a collated cell_seg_only folder
  subDir <- "Cell_seg_only"
  if (!file.exists(subDir)){
    dir.create(file.path(folderoffolders, subDir))
  }
  # Copy from the origin, to the target
  print("Copying Files to New Directory")
  file.copy(from = filelist, to = targetdir, recursive = F, 
            copy.mode = T)
  # New directory for CSVs
  setwd(targetdir)
  subDir1 <- "CSV"
  
  if (!file.exists(subDir1)){
    dir.create(file.path(targetdir, subDir1))
  }
  filelist <- list.files(targetdir, pattern = ".txt$")
  print("Start conversion to CSV")
  ## Cycle through converting to CSV
  i <- 1
  
  for (i in 1:length(filelist)) {
    cur.input.file <- filelist[i]
    print(cur.input.file)
    cur.output.file <- paste0("CSV/",cur.input.file, ".csv") 
    print(cur.output.file)
    print(paste("Working on file:", cur.input.file))
    data <- read.delim(cur.input.file, header = T)
    write.table(data, file = cur.output.file, sep = ",", col.names = T, row.names = F)}
  # Bind all together
  options(stringsAsFactors = F)
  files <- list.files(subTdir, pattern = ".csv$")
  setwd(subTdir)
  lists <- lapply(files, read.csv, header = T)
  print("Binding together")
  combined.df <- rbindlist(lists, fill = T)
  return(combined.df)
  setwd(origindir)
}

## cell_seg_summary
preprocess_cell_seg_summary <- function(folderoffolders, targetdir, subTdir, origindir){
  print("Setting up Directories")
  setwd(folderoffolders)
  thousand.folders <- list.dirs(full.names = T)
  filelist1 <- sapply(thousand.folders[-1], function(x)
  {list.files(x, pattern = "]_cell_seg_data_summary.txt$", full.names = T)
  })
  filelist = unlist(filelist1)
  # Create a collated cell_seg_only folder
  subDir <- "Cell_seg_summary_only"
  if (!file.exists(subDir)){
    dir.create(file.path(folderoffolders, subDir))
  }
  # Copy from the origin, to the target
  print("Copying to new directory")
  file.copy(from = filelist, to = targetdir, recursive = F, 
            copy.mode = T)
  # New directory for CSVs
  setwd(targetdir)
  subDir1 <- "CSV"
  
  if (!file.exists(subDir1)){
    dir.create(file.path(targetdir, subDir1))
  }
  filelist <- list.files(targetdir, pattern = ".txt$")
  print("Start CSV conversion")
  ## Cycle through converting to CSV
  i <- 1
  
  for (i in 1:length(filelist)) {
    cur.input.file <- filelist[i]
    print(cur.input.file)
    cur.output.file <- paste0("CSV/",cur.input.file, ".csv") 
    print(cur.output.file)
    print(paste("Working on file:", cur.input.file))
    data <- read.delim(cur.input.file, header = T)
    write.table(data, file = cur.output.file, sep = ",", col.names = T, row.names = F)}
  # Bind all together
  options(stringsAsFactors = F)
  files <- list.files(subTdir, pattern = ".csv$")
  setwd(subTdir)
  lists <- lapply(files, read.csv, header = T)
  print("Binding together")
  combined_summary_df <- rbindlist(lists, fill = T)
  writeCsvD(combined_summary_df)
  setwd(origindir)
}

# Cleancellseg
clean_merged_cell_seg <- function(df){
  print("Remove Columns")
  df.combined <- df %>% 
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
  print("Fix Columns")
  df.combined$Confidence <- as.numeric(sub("%", "", df.combined$Confidence, fixed = T))
  df1 <- df.combined 
  df1$Tissue.Category <- fix_tissue_cat(df1)
  df1$Phenotype <- as.character(df1$Phenotype)
  df1$Phenotype[df1$Phenotype == ""] <- "DAPI"
  df1$Phenotype <- fix_pheno(df1)
  print("Remove Background")
  df2 <- filter(df1, Tissue.Category != "Background")
  df2a <- df2
  print("Return cleaned data")
  writeCsvO(df2a)
}

# cleancellseg_batch2
clean_merged_cell_seg1 <- function(df){
  print("Remove Columns")
  df.combined <- df %>% 
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
  print("Fix Columns")
  df.combined$Confidence <- as.numeric(sub("%", "", df.combined$Confidence, fixed = T))
  df1 <- df.combined 
  df1$Tissue.Category <- fix_tissue_cat(df1)
  df1$Phenotype <- as.character(df1$Phenotype)
  df1$Phenotype[df1$Phenotype == ""] <- "DAPI"
  df1$Phenotype <- fix_pheno(df1)
  print("Remove Background")
  df2 <- filter(df1, Tissue.Category != "Background")
  df2b <- df2
  print("Return cleaned data")
  writeCsvO(df2b)
}

# Read and bind
Bind_them <- function(folder){
  print("Listing files")
  options(stringsAsFactors = F)
  files <- list.files(folder, pattern = ".csv$")
  setwd(folder)
  print("Reading files")
  lists <- lapply(files, read.csv, header = T)
  print("Binding files")
  combined.df <- rbindlist(lists)
  setwd(c("/Users/jlm650/OneDrive/UoB/PhD/1st_Year/Projects/5_Extra/Nahla_Analysis"))
  print("Writing bound file")
  writeCsvD(combined)
}

# Fixing discrepancies
## Fix tissue.Category column
fix_tissue_cat <- function(df){
  df$Tissue.Category %<>%
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
  return(df$Tissue.Category)
}

## Fix Phenotype Spellings
fix_pheno <- function(df){
  df$Phenotype %<>%
    gsub("cd4", "CD4", .) %>%
    gsub("cd68", "CD68", .) %>%
    gsub("dapi", "DAPI", .)   %>% 
    gsub("nuclei", "DAPI", .) %>%
    gsub("Nuclei", "DAPI", .) %>%
    gsub("nucleui", "DAPI", .) %>%
    gsub("foxp3", "FOXP3", .) %>%
    gsub("Foxp3", "FOXP3", .) %>%
    gsub("FOXp3", "FOXP3", .) %>%
    gsub("CD63", "CD68", .) %>%
    gsub("^$", "DAPI", .)
  return(df$Phenotype)
}

# Find diameter
find_diameter <- function(df){
  radiussq <- df$Entire.Cell.Area..pixels./pi
  radius <- sqrt(radiussq)
  diameter <- radius*2
  return(diameter)
}

# Find combinations of phenotypes
find_pheno_comb <- function(Phenotypes){
  Phenotype1 <- Phenotypes
  Phenotype2 <- Phenotypes
  combin <- as.data.frame(cbind(Phenotype1, Phenotype2))
  combin$Phenotype1 <- as.factor(Phenotype1)
  combin$Phenotype2 <- as.factor(Phenotype2)
  combinations <- expand.grid(levels(combin$Phenotype1), levels(combin$Phenotype2))
  uniqcomb <- combinations[!combinations$Var1==combinations$Var2,]
  pairs = list()
  for(i in seq_len(nrow(uniqcomb))){
    here <- as.vector(unlist(uniqcomb[i,]))
    pairs[[length(pairs)+1]] <- here
  }
  return(pairs)
}

# Factor a list of given columns
factorthese <- function(df, somecolumns){
  Fctr <- names(df) %in% somecolumns
  df[,Fctr] <- lapply(df[,Fctr], function(column) as.factor(as.character(column)))
  return(df)
}

# Trim whitespace
trim <- function (x){
  gsub("^\\s+|\\s+$", "", x)}

# returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)

# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)


# Convert CSV to TSV
csvtotsv <- function(bp, op){
  print("Reading")
  df1 <- read.csv(bp)
  names(df1) <- gsub("\\.", " ", names(df1))
  print("Writing")
  write.table(df1, file = op, row.names = F, sep="\t")
}

# Write a csv for raw data
writeCsvD <- function(df){
  fn <- deparse(substitute(df))
  write.csv(df, file = paste0("/Users/jlm650/OneDrive/UoB/PhD/1st_Year/Projects/5_Extra/Nahla_Analysis/Data/", fn, ".csv"), row.names = F)}

# Write a csv for non raw data
writeCsvO <- function(df){
  fn <- deparse(substitute(df))
  write.csv(df, file = paste0("/Users/jlm650/OneDrive/UoB/PhD/1st_Year/Projects/5_Extra/Nahla_Analysis//Output/",fn,".csv"), row.names = F)}

# Function to remove the addition of X. column at the start of a dataframe
Correct_Colnames <- function(df) {
  delete.columns <- grep("(^X$)|(^X\\.)(\\d+)($)", colnames(df), perl=T)
  if (length(delete.columns) > 0) {
    row.names(df) <- as.character(df[, grep("^X$", colnames(df))])
    df <- df[,-delete.columns]
    colnames(df) <- gsub("^X", "",  colnames(df))
  }
  return(df)
}

# Make a colour transparent (use in a vector of colours, do not make new object)
makeTransparent<-function(someColor, alpha=25){
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

# Add regression formula
lm_eqn <- function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}

# gg plot for comparing Ranks across subtype
ggSubtype <- function(df, YTitle, MainTitle){
  p <- ggplot(df,aes(x = Subtype, y = Rank)) +
    geom_boxplot(
      alpha = 0.5,
      width = 0.2)+ 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width",
                alpha = 0.8) +
    scale_fill_manual(values = cbcols) +
    labs(x = "Subtype",
         y = YTitle,
         title = MainTitle
    )+ 
    geom_dotplot(binaxis='y',
                 method="histodot",
                 stackdir='center',
                 binwidth = 20,
                 position=position_jitter(0.1),
                 alpha=0,
                 dotsize=0.4)+
    theme_bw()+
    theme(
      axis.text = element_text(size=16))
  sig <- p+stat_compare_means(comparisons = my_comparisons,
                              label = "p.signif")
  return(sig)}

# Find the correlation between variable and PC
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}

# Calculate the contribution of variable to the PC
contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}

# GGdist
ggDist <- function(df, YTitle, MainTitle){
  p <- ggplot(df,aes(x = Subtype, y = Max)) +
    geom_boxplot(
      alpha = 0.5,
      width = 0.2)+ 
    geom_violin(aes(Subtype, fill = Subtype),
                scale = "width",
                alpha = 0.8) +
    scale_color_manual(values = cbcols) +
    labs(x = "Subtype",
         y = , YTitle, 
         main = MainTitle
    )+ 
    geom_dotplot(binaxis='y',
                 method="histodot",
                 stackdir='center',
                 binwidth = 20,
                 position=position_jitter(0.1),
                 alpha=0,
                 dotsize=0.4)+
    theme_bw()+
    theme(
      axis.text = element_text(size=16)) +
    theme(legend.direction = 'horizontal', 
          legend.position = 'top')
  sig <- p+stat_compare_means(comparisons = my_comparisons,
                              label = "p.signif")
  return(sig)}





