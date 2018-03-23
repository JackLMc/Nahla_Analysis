# Source Functions + Packages
source("Functions.R")

# Only needs to be ran once 
# Cell_seg
## Set the directories
folderoffolders <- c("/Volumes/ResearchData/Vectra/Vectra3/Jack/CRC_Project/MSS_hiCIRC")
targetdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Jack/CRC_Project/MSS_hiCIRC/Cell_seg_only")
subTdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Jack/CRC_Project/MSS_hiCIRC/Cell_seg_only/CSV")
origindir <- c("/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC")

# Preprocess the folders, convert .txt to .csv and bind all
preprocess_cell_seg(folderoffolders, targetdir, subTdir, origindir)

# Write to allow reading
combined.df <- read.csv("Data/combined.df.csv")
clean_merged_cell_seg(combined.df)


# Cell_seg_summary
## Scan through folders and separate cell_seg_data_summary.txt files
## Set new directories
targetdir1 <- c("/Volumes/ResearchData/Vectra/Vectra3/Jack/CRC_Project/MSS_hiCIRC/Cell_seg_summary_only")
subTdir1 <- c("/Volumes/ResearchData/Vectra/Vectra3/Jack/CRC_Project/MSS_hiCIRC/Cell_seg_summary_only/CSV")

preprocess_cell_seg_summary(folderoffolders, targetdir1, subTdir1, origindir)
