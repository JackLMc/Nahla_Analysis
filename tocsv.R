# Source Functions + Packages
source("Functions.R")

# First batch
# Cell_seg
## Set the directories
folderoffolders <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports")
targetdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Cell_seg_only")
subTdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Cell_seg_only/CSV")
origindir <- c("/Users/jlm650/OneDrive/UoB/PhD/1st_Year/Projects/5_Extra/Nahla_Analysis/")

# Preprocess the folders, convert .txt to .csv and bind all
combined.df <- preprocess_cell_seg(folderoffolders, targetdir, subTdir, origindir)
writeCsvD(combined.df)

# Write to allow reading
setwd(origindir)
combined.df <- read.csv("Data/combined.df.csv")

# # Number of images per slide (cell seg)
# CD1a <- factorthese(combined.df, c("Sample.Name", "Tissue.Category", "Slide.ID", "Phenotype"))
# 
# Number_of_Images_CS <- data.frame(Slide = character(),
#                      Number_of_Fields_CS = double(),
#                      stringsAsFactors = F)
# c <- 1
# for(i in levels(CD1a$Slide.ID)){
#   print(i)
#   work <- droplevels(subset(CD1a, Slide.ID == i))
#   nom <- nlevels(work$Sample.Name)
#   Number_of_Images_CS[c, "Slide"] <- i
#   Number_of_Images_CS[c, "Number_of_Fields_CS"] <- nom
#   c <- c + 1
# }
# writeCsvO(Number_of_Images_CS)

clean_merged_cell_seg(combined.df)
df2a <- read.csv("Output/df2a.csv")
df3 <- df2a[, c("Slide.ID", "Sample.Name", "Tissue.Category", "Phenotype", "Cell.ID", "Cell.X.Position", "Cell.Y.Position")]
writeCsvO(df3)

# Batch 2
folderoffolders <- c("/Volumes/ResearchData/Vectra/Vectra3/chris folder FEDOR/Nahla new results/")
targetdir <- c("/Volumes/ResearchData/Vectra/Vectra3/chris folder FEDOR/Nahla new results/Cell_seg_only")
subTdir <- c("/Volumes/ResearchData/Vectra/Vectra3/chris folder FEDOR/Nahla new results/Cell_seg_only/CSV")
origindir <- c("/Users/jlm650/OneDrive/UoB/PhD/1st_Year/Projects/5_Extra/Nahla_Analysis/")

# Preprocess the folders, convert .txt to .csv and bind all
combined.df1 <- preprocess_cell_seg(folderoffolders, targetdir, subTdir, origindir)
writeCsvD(combined.df1)

# Write to allow reading
setwd(origindir)
combined.df1 <- read.csv("Data/combined.df1.csv")

clean_merged_cell_seg1(combined.df1)
df2b <- read.csv("Output/df2b.csv")
df4 <- df2b[, c("Slide.ID", "Sample.Name", "Tissue.Category", "Phenotype", "Cell.ID", "Cell.X.Position", "Cell.Y.Position")]
writeCsvO(df4)











# Fix the raw files (for Radius)
txt_list <- list()
c <- 1
for(i in levels(df2a$Sample.Name)){
  name <- basename(i)
  cat('Processing', i, '\n')
  txt_list[[i]] <- droplevels(subset(df2a, Sample.Name == i))
  c <- c + 1}

for(i in names(txt_list)){
  setwd("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Export_Fixed/Cell_seg_only")
  print(i)
  write.table(txt_list[[i]], paste0(i,".txt"), quote = F, sep = "\t", row.names = F)
}


# Cell_seg_summary
## Scan through folders and separate cell_seg_data_summary.txt files
## Set new directories
targetdir1 <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Cell_seg_summary_only")
subTdir1 <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Cell_seg_summary_only/CSV")

preprocess_cell_seg_summary(folderoffolders, targetdir1, subTdir1, origindir)
combined_summary_df <- read.csv("Data/combined_summary_df.csv")

# # Number of images per slide
# CD1 <- factorthese(combined_summary_df, c("Sample.Name", "Tissue.Category", "Slide.ID", "Phenotype"))
# 
# Number_of_Images <- data.frame(Slide = character(),
#                      Number_of_Fields_CSS = double(),
#                      stringsAsFactors = F)
# c <- 1
# for(i in levels(CD1$Slide.ID)){
#   print(i)
#   work <- droplevels(subset(CD1, Slide.ID == i))
#   nom <- nlevels(work$Sample.Name)
#   Number_of_Images[c, "Slide"] <- i
#   Number_of_Images[c, "Number_of_Fields_CSS"] <- nom
#   c <- c + 1
# }
# writeCsvO(Number_of_Images)

# Fix discrepancies in spelling
combined_summary_df$Tissue.Category <- fix_tissue_cat(combined_summary_df)
combined_summary_df$Phenotype <- fix_pheno(combined_summary_df)

writeCsvD(combined_summary_df)
