# Source Functions + Packages
library(UsefulFunctions)
library(tidyverse)
library(data.table)

# Read data / source("tocsv.R)
# Reading in data and cleaning up ----
cell_seg_summary <- read.csv("Output/Summary_data.csv")
cell_seg_summary$Tissue.Category <- as.factor(cell_seg_summary$Tissue.Category)

# Remove unneeded stuff
CS1 <- droplevels(subset(cell_seg_summary, Tissue.Category != "Background"))
CS1a <- droplevels(subset(CS1, Phenotype != "All"))
keep <- c("Sample.Name",
          "Tissue.Category",
          "Tissue.Category.Area..pixels.",
          "Slide.ID",
          "Phenotype",
          "Total.Cells",
          "Cell.Density..per.megapixel.")
CS1b <- CS1a[,(names(CS1a) %in% keep)]
CS1b$Sample.Name <- as.factor(CS1b$Sample.Name)
CS1b$Tissue.Category <- as.factor(CS1b$Tissue.Category)
CS1b$Phenotype <- as.factor(CS1b$Phenotype)

# Create a unique row
CS1b$uniq <- as.factor(paste(CS1b$Slide.ID, CS1b$Tissue.Category, CS1b$Phenotype, sep = "/"))
CS1c <- CS1b

CS1c[is.na(CS1c$Total.Cells)] <- 0
CS1c$Total.Cells <- as.numeric(CS1c$Total.Cells)


output <- data.frame(Combine = character(),
                     Total.Cells = double(),
                     stringsAsFactors = FALSE)

c <- 1
for(i in levels(CS1c$uniq)){
  print(paste("Working on:", i))
  Working <- droplevels(subset(CS1c, uniq == i))
  TotalC <- sum(Working$Total.Cells)
  output[c, "Combine"] <- i
  output[c, "Total.Cells"] <- TotalC
  c <- c + 1}


head(output)



# Check for missing phenotypes from pictures and create dataframe ----
# List missing phenotypes from pictures
missing_pheno <- list()

for(i in levels(CS1b$Sample.Name)){
  print(i)
  work <- droplevels(subset(CS1b, Sample.Name == i))
  pat_pheno <- as.character(levels(work$Phenotype))
  all_pheno <- as.character(levels(CS1b$Phenotype))
  see <- all_pheno[!(all_pheno %in% pat_pheno)]
  missing <- ifelse((see != 0), c(see), "None")
  missing_pheno[[i]] <- missing
}

# Create dataframe and change column names
missing_df <- plyr::ldply(missing_pheno, data.frame)
colnames(missing_df) <- c("Sample.Name", "Phenotype")
head(missing_df)
# Replicate each element three times (each tissue.category - All, Epi, Stroma)
df.expanded <- missing_df[rep(seq_len(nrow(missing_df)), each = 3),]
Tissue.Category <- c("All", "Epithelium", "Stroma")
Cell.Density..per.megapixel. <- 0 
Total.Cells <- 0
missing_data <- as.data.frame(cbind(df.expanded, Tissue.Category, Total.Cells, Cell.Density..per.megapixel.))
missing_data$Slide.ID <- word(missing_data$Sample.Name, 1, sep = "_")


head(missing_data)





# Create a new dataframe for average density
CS1c[is.na(CS1c$Cell.Density..per.megapixel.)] <- 0
CS1c$Cell.Density..per.megapixel. <- as.numeric(CS1c$Cell.Density..per.megapixel.)
output <- data.frame(Combine = character(),
                     Average.cell.Density = double(),
                     stringsAsFactors = F)


# check <- data.frame(Combine = character(),
#                    num_image = double(),
#                    stringsAsFactors = F)
# 
# c <- 1
# for(i in levels(CS1c$uniq)){
#   work <- droplevels(subset(CS1c, uniq == i))
#   num_samp <- nlevels(work$Sample.Name)
#   check[c, "Combine"] <- i
#   check[c, "num_image"] <- num_samp
#   c <- c + 1
# }


c <- 1
for(i in levels(CS1c$uniq)){
  print(paste("Working on:", i))
  Working <- droplevels(subset(CS1c, uniq == i))
  TotalC <- mean(Working$Cell.Density..per.megapixel.)
  output[c, "Combine"] <- i
  output[c, "Average.cell.Density"] <- TotalC
  c <- c + 1}

# Separate and factor
CS1d <- output %>%
  separate(Combine, c("Slide.ID", "Tissue.Category", "Phenotype"), "/")
CS1d$Slide.ID <- as.factor(CS1d$Slide.ID)
CS1d$Tissue.Category <- as.factor(CS1d$Tissue.Category)
CS1d$Phenotype <- as.factor(CS1d$Phenotype)

CS1d$Rank <- rank(CS1d$Average.cell.Density)

# Paste together Tissue.Category and Phenotype
CS1d$Parameter <- as.factor(paste("Average Density of ", CS1d$Phenotype, " Positive Cells in the ", CS1d$Tissue.Category, sep = ""))
levels(CS1d$Parameter)

CS1d$Slide.ID <- gsub("  5PLEX|5plex Fedor| 5PLEX| breast tumour CD4 CD8 CD20 CD68 FOXP3", "", CS1d$Slide.ID) 
CS1d$Slide.ID <- trim.trailing(CS1d$Slide.ID) %>% as.factor()


clinical <- read.csv("./Data/Nahla_clinical.csv")
colnames(clinical)[colnames(clinical) == "Lablels.as.in.inform.cases"] <- "Slide.ID"
clinical <- droplevels(subset(clinical, Slide.ID != ""))
clinical$ER[clinical$ER == "Negative"] <- "negative"
clinical$HER2[clinical$HER2 == "Positive"] <- "positive"
clinical$Path.Response[clinical$Path.Response == "pCR"] <- "PCR"
clinical <- droplevels(clinical)

# 
data_clin <- merge(CS1d, clinical, by = "Slide.ID") %>% droplevels()
data_clin$Slide.ID <- as.factor(data_clin$Slide.ID)

# levels(clinical$Slide.ID)[!('%in%'(levels(clinical$Slide.ID), levels(CS1d$Slide.ID)))]
# levels(CS1d$Slide.ID)[grep("S073561", levels(CS1d$Slide.ID))]
head(data_clin)

data_clin1 <- data_clin[, c("Slide.ID", "Average.cell.Density", "Parameter", "Grade", "HER2", "ER", "Path.Response")]
data_clin2 <- spread(data_clin1, key = "Parameter", value = "Average.cell.Density") %>% column_to_rownames(., var = "Slide.ID")

head(data_clin2)

pca1 <- data_clin2 %>% dplyr:: select(-contains("All")) %>% dplyr:: select(-contains("DAPI"))
pca1[is.na(pca1)] <- 0


prin_comp <- prcomp(pca1[, names(pca1) != "Grade" & names(pca1) != "HER2" & names(pca1) != "ER" & names(pca1) != "Path.Response"], scale. = T)
Subtype <- pca1[, "Grade"] %>% as.factor()


library(ggbiplot)
g <- ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
         groups = Subtype, ellipse = T,
         circle = T,
         var.axes = F,
         choices = c(1, 2)) +
  #scale_color_manual(values = cbcols) +
  theme_bw() + 
  theme(legend.direction = 'horizontal', 
               legend.position = 'top')
ggsave("PCA of Immune Cell Abundance (Density per Mpix).pdf" , plot = g, device = "pdf",
       path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Cell_Density/",
       height = 6, width = 6, units = 'in', dpi = 600)





# Plot
cbcols <- c("#999999",
            "#56B4E9",
            "#E69F00",
            "#009E73")

data_clin3 <- data_clin2 %>% gather(contains("Positive Cells in the"), key = "Parameter", value = "Average.cell.Density")
data_clin4 <- factorthese(data_clin3, c("Grade", "HER2", "ER", "Path.Response", "Parameter"))



for(i in levels(data_clin4$Parameter)){
  print(i)
  df <- droplevels(subset(data_clin4, Parameter == i))
  p <- ggplot(df, aes(x = ER, y = Average.cell.Density)) +
    geom_boxplot(alpha = 0.5, width = 0.2) +
    geom_violin(aes(ER, fill = ER), scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) + labs(x = "ER", y = i) +
    geom_dotplot(binaxis = "y",
                 method = "histodot",
                 stackdir = "center",
                 binwidth = 20,
                 position = position_jitter(0.1),
                 alpha = 0,
                 dotsize = 0.4)+
    theme_bw()+
    theme(axis.text = element_text(size = 16)) +
    stat_compare_means(comparisons = list(c("negative", "positive")), label = "p.signif")
  filen <- paste0(i,".pdf")
  ggsave(filen, plot = p, device = "pdf",
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Cell_Density/ER",
         height = 5, width = 5, units = "in")
}


for(i in levels(data_clin4$Parameter)){
  print(i)
  df <- droplevels(subset(data_clin4, Parameter == i))
  p <- ggplot(df, aes(x = HER2, y = Average.cell.Density)) +
    geom_boxplot(alpha = 0.5, width = 0.2) +
    geom_violin(aes(HER2, fill = HER2), scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) + labs(x = "HER2", y = i) +
    geom_dotplot(binaxis = "y",
                 method = "histodot",
                 stackdir = "center",
                 binwidth = 20,
                 position = position_jitter(0.1),
                 alpha = 0,
                 dotsize = 0.4)+
    theme_bw()+
    theme(axis.text = element_text(size = 16)) +
    stat_compare_means(comparisons = list(c("negative", "positive")), label = "p.signif")
  filen <- paste0(i,".pdf")
  ggsave(filen, plot = p, device = "pdf",
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Cell_Density/HER2",
         height = 5, width = 5, units = "in")
}

for(i in levels(data_clin4$Parameter)){
  print(i)
  df <- droplevels(subset(data_clin4, Parameter == i))
  p <- ggplot(df, aes(x = Grade, y = Average.cell.Density)) +
    geom_boxplot(alpha = 0.5, width = 0.2) +
    geom_violin(aes(Grade, fill = Grade), scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) + labs(x = "Grade", y = i) +
    geom_dotplot(binaxis = "y",
                 method = "histodot",
                 stackdir = "center",
                 binwidth = 20,
                 position = position_jitter(0.1),
                 alpha = 0,
                 dotsize = 0.4)+
    theme_bw()+
    theme(axis.text = element_text(size = 16)) +
    stat_compare_means(comparisons = list(c("1", "2"),
                                          c("2", "3"),
                                          c("1", "3")), label = "p.signif")
  filen <- paste0(i,".pdf")
  ggsave(filen, plot = p, device = "pdf",
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Cell_Density/Grade",
         height = 5, width = 5, units = "in")
}


