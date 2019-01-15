
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
    write.table(data, file = cur.output.file, sep = ",", col.names = TRUE, row.names = FALSE)}
  # Bind all together
  options(stringsAsFactors = F)
  files <- list.files(subTdir, pattern = ".csv$")
  setwd(subTdir)
  lists <- lapply(files, read.csv, header = TRUE)
  print("Binding together")
  combined.df <- rbindlist(lists)
  #writeCsvD(combined.df)
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
  file.copy(from = filelist, to = targetdir, recursive = FALSE, 
            copy.mode = TRUE)
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
    write.table(data, file = cur.output.file, sep = ",", col.names = TRUE, row.names = FALSE)}
  # Bind all together
  options(stringsAsFactors = F)
  files <- list.files(subTdir, pattern = ".csv$")
  setwd(subTdir)
  lists <- lapply(files, read.csv, header = TRUE)
  print("Binding together")
  combined_summary_df <- rbindlist(lists)
  #writeCsvD(combined_summary_df)
  return(combined_summary_df)
  setwd(origindir)
}





folderoffolders <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/new_export_directory")
targetdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/new_export_directory/Cell_seg_only")
subTdir <- c("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/new_export_directory/Cell_seg_only/CSV")
origindir <- c("/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Extra/Nahla_Analysis")

# Preprocess the folders, convert .txt to .csv and bind all
check <- preprocess_cell_seg(folderoffolders, targetdir, subTdir, origindir)
check <- preprocess_cell_seg_summary(folderoffolders, targetdir, subTdir, origindir)

ncol(check)
setwd(origindir)
