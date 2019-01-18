# Source Functions + Packages
library(UsefulFunctions)
library(tidyverse)
library(data.table)

# First batch
# Cell_seg
## Set the directories
folderoffolders <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Batch_1")
targetdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Batch_1/Cell_seg_only")
subTdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Batch_1/Cell_seg_only/CSV")
origindir <- c("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/")

# Preprocess the folders, convert .txt to .csv and bind all
combined.df <- preprocess_cell_seg(folderoffolders, targetdir, subTdir, origindir)
write.csv(file = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Data/Batch1_single.csv", x = combined.df, row.names = F)

# Write to allow reading
setwd(origindir)
Batch1_single <- read.csv("Data/Batch1_single.csv")

df1 <- clean_merged_cell_seg(Batch1_single)
Batch1_single <- df1[, c("Slide.ID", "Sample.Name", "Tissue.Category", "Phenotype", "Cell.ID", "Cell.X.Position", "Cell.Y.Position")]
write.csv(file = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Output/Batch1_single.csv", x = Batch1_single, row.names = F)

# Batch 2
folderoffolders <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Batch_2")
targetdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Batch_2/Cell_seg_only")
subTdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Batch_2/Cell_seg_only/CSV")
origindir <- c("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/")

# Preprocess the folders, convert .txt to .csv and bind all
combined.df <- preprocess_cell_seg(folderoffolders, targetdir, subTdir, origindir)

# Write to allow reading
setwd(origindir)
Batch2_single <- read.csv("Data/Batch2_single.csv")

df1 <- clean_merged_cell_seg(Batch2_single)
Batch2_single <- df1[, c("Slide.ID", "Sample.Name", "Tissue.Category", "Phenotype", "Cell.ID", "Cell.X.Position", "Cell.Y.Position")]
write.csv(file = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Output/Batch2_single.csv", x = Batch2_single, row.names = F)

# Number of images per slide (cell seg)
# CD1a <- factorthese(Batch1_single, c("Sample.Name", "Tissue.Category", "Slide.ID", "Phenotype"))
# 
# Number_of_Images <- data.frame(Slide = character(),
#                      Number_of_Fields_sin = double(),
#                      stringsAsFactors = F)
# c <- 1
# for(i in levels(CD1a$Slide.ID)){
#   print(i)
#   work <- droplevels(subset(CD1a, Slide.ID == i))
#   nom <- nlevels(work$Sample.Name)
#   Number_of_Images_CS[c, "Slide"] <- i
#   Number_of_Images_CS[c, "Number_of_Fields_sin"] <- nom
#   c <- c + 1
# }
# 
# CD1a <- factorthese(Batch2_single, c("Sample.Name", "Tissue.Category", "Slide.ID", "Phenotype"))
# 
# Number_of_Images_CS <- data.frame(Slide = character(),
#                                   Number_of_Fields_sin = double(),
#                                   stringsAsFactors = F)
# c <- 1
# for(i in levels(CD1a$Slide.ID)){
#   print(i)
#   work <- droplevels(subset(CD1a, Slide.ID == i))
#   nom <- nlevels(work$Sample.Name)
#   Number_of_Images[c, "Slide"] <- i
#   Number_of_Images[c, "Number_of_Fields_sin"] <- nom
#   c <- c + 1
# }
# 
# single_pics <- rbind(Number_of_Images_CS, Number_of_Images)
# write.csv(file = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Output/Single_pic_check.csv", x = single_pics, row.names = F)


# Source Functions + Packages
library(UsefulFunctions)
library(tidyverse)
library(data.table)

# First batch
# Cell_seg_summary
## Set the directories
folderoffolders <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Batch_1")
targetdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Batch_1/Cell_seg_summary_only")
subTdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Batch_1/Cell_seg_summary_only/CSV")
origindir <- c("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/")

# Preprocess the folders, convert .txt to .csv and bind all
combined.df <- preprocess_cell_seg_summary(folderoffolders, targetdir, subTdir, origindir)
write.csv(file = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Data/Batch1_summary.csv", x = combined.df, row.names = F)

# Write to allow reading
Batch1_summary <- read.csv("Data/Batch1_summary.csv")
Batch1_summary_clean <- Batch1_summary[, c("Sample.Name", "Tissue.Category", "Phenotype",
                                           "Total.Cells", "Tissue.Category.Area..pixels.", "Cell.Density..per.megapixel.",
                                           "Slide.ID", "Confidence")]
Batch1_summary_clean$Confidence <- as.numeric(sub("%", "", Batch1_summary_clean$Confidence, fixed = T))
df1 <- Batch1_summary_clean
df1$Tissue.Category <- fix_tissue_cat(df1)
df1$Phenotype <- as.character(df1$Phenotype)
df1$Phenotype[df1$Phenotype == ""] <- "DAPI"
df1$Phenotype <- fix_pheno(df1)
df2 <- filter(df1, Tissue.Category != "Background")

Batch1_summary <- df2
write.csv(file = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Output/Batch1_summary.csv", x = Batch1_summary, row.names = F)

# Batch 2
folderoffolders <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Batch_2")
targetdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Batch_2/Cell_seg_summary_only")
subTdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Batch_2/Cell_seg_summary_only/CSV")
origindir <- c("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/")

# Preprocess the folders, convert .txt to .csv and bind all
combined.df <- preprocess_cell_seg_summary(folderoffolders, targetdir, subTdir, origindir)
write.csv(file = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Data/Batch2_summary.csv", x = combined.df, row.names = F)

# Write to allow reading
Batch2_summary <- read.csv("Data/Batch2_summary.csv")
Batch2_summary_clean <- Batch2_summary[, c("Sample.Name", "Tissue.Category", "Phenotype",
                                           "Total.Cells", "Tissue.Category.Area..pixels.", "Cell.Density..per.megapixel.",
                                           "Slide.ID", "Confidence")]
Batch2_summary_clean$Confidence <- as.numeric(sub("%", "", Batch2_summary_clean$Confidence, fixed = T))
df1 <- Batch2_summary_clean
df1$Tissue.Category <- fix_tissue_cat(df1)
df1$Phenotype <- as.character(df1$Phenotype)
df1$Phenotype[df1$Phenotype == ""] <- "DAPI"
df1$Phenotype <- fix_pheno(df1)
df2 <- filter(df1, Tissue.Category != "Background")

Batch2_summary <- df2
write.csv(file = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Output/Batch2_summary.csv", x = Batch2_summary, row.names = F)



# Bind into big summary file
summary <- rbind(Batch1_summary, Batch2_summary)
write.csv("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Output/Summary_data.csv", x = summary)


# Number of images per slide (summary)
# CD1 <- factorthese(summary, c("Sample.Name", "Tissue.Category", "Slide.ID", "Phenotype"))
# 
# Number_of_Images <- data.frame(Slide = character(),
#                                Number_of_Fields_CSS = double(),
#                                stringsAsFactors = F)
# c <- 1
# for(i in levels(CD1$Slide.ID)){
#   print(i)
#   work <- droplevels(subset(CD1, Slide.ID == i))
#   nom <- nlevels(work$Sample.Name)
#   Number_of_Images[c, "Slide"] <- i
#   Number_of_Images[c, "Number_of_Fields_CSS"] <- nom
#   c <- c + 1
# }
# 
# write.csv(file = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Output/Summary_pic_check.csv", x = Number_of_Images, row.names = F)

