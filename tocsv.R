# Source Functions + Packages
source("Functions.R")

# Only needs to be ran once 
# Cell_seg
## Set the directories
folderoffolders <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports")
targetdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Cell_seg_only")
subTdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Cell_seg_only/CSV")
origindir <- c("/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Extra/Nahla_Analysis")

# Preprocess the folders, convert .txt to .csv and bind all
preprocess_cell_seg(folderoffolders, targetdir, subTdir, origindir)

# Write to allow reading
combined.df <- read.csv("Data/combined.df.csv")
clean_merged_cell_seg(combined.df)


# Cell_seg_summary
## Scan through folders and separate cell_seg_data_summary.txt files
## Set new directories
targetdir1 <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Cell_seg_summary_only")
subTdir1 <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Cell_seg_summary_only/CSV")

preprocess_cell_seg_summary(folderoffolders, targetdir1, subTdir1, origindir)

CD1 <- factorthese

Number_of_Images <- data.frame(Slide = character(),
                     Number_of_Fields = double(),
                     stringsAsFactors = F)
c <- 1
for(i in levels(CD1$Slide.ID)){
  print(i)
  work <- droplevels(subset(CD1, Slide.ID == i))
  nom <- nlevels(work$Sample.Name)
  Number_of_Images[c, "Slide"] <- i
  Number_of_Images[c, "Number_of_Fields"] <- nom
  c <- c + 1
}
writeCsvO(Number_of_Images)
