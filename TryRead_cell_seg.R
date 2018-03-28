# Read using read_cell_seg_data
source("Functions.R")

# Only needs to be ran once 

# Preprocess
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
  file.copy(from = filelist, to = targetdir, recursive = FALSE, 
            copy.mode = TRUE)

  filelist <- list.files(targetdir, pattern = ".txt$")
  # Bind all together
  options(stringsAsFactors = F)
  files <- list.files(targetdir, pattern = ".csv$")
  setwd(subTdir)
  lists <- lapply(files, read_cell_seg_data, header = TRUE)
  print("Binding together")
  combined.df <- rbindlist(lists)
  return(combined.df)
  setwd(origindir)
}

options(stringsAsFactors = F)
files <- list.files(targetdir, pattern = ".txt$")
setwd(targetdir)
lists <- lapply(files, read_cell_seg_data)
print("Binding together")
combined.df <- rbindlist(lists)




# Cell_seg
## Set the directories
folderoffolders <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports")
targetdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Cell_seg_only")
subTdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Exports/Cell_seg_only/CSV")
origindir <- c("/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Extra/Nahla_Analysis")

# Preprocess the folders, convert .txt to .csv and bind all
preprocess_cell_seg(folderoffolders, targetdir, subTdir, origindir)
setwd(origindir)

# Write to allow reading
combined.df <- read.csv("Data/combined.df.csv")
