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

# Replicate each element three times (each tissue.category - All, Epi, Stroma)
df.expanded <- missing_df[rep(seq_len(nrow(missing_df)), each = 3),]
Tissue.Category <- c("All", "Epithelium", "Stroma")
Cell.Density..per.megapixel. <- 0 
Total.Cells <- 0
missing_data <- as.data.frame(cbind(df.expanded, Tissue.Category, Total.Cells, Cell.Density..per.megapixel.))
missing_data$Slide.ID <- word(missing_data$Sample.Name, 1, sep = "_")


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

data_clin1 <- data_clin[, c("Slide.ID", "Average.cell.Density", "Parameter", "Grade", "HER2", "ER", "Path.Response", "survival.in.months", "Dead.vs.alive")]
data_clin2 <- spread(data_clin1, key = "Parameter", value = "Average.cell.Density") %>% column_to_rownames(., var = "Slide.ID")

head(data_clin2)

# pca1 <- data_clin2 %>% dplyr:: select(-contains("All")) %>% dplyr:: select(-contains("DAPI"))
# pca1[is.na(pca1)] <- 0
# 
# 
# prin_comp <- prcomp(pca1[, names(pca1) != "Grade" & names(pca1) != "HER2" & names(pca1) != "ER" & names(pca1) != "Path.Response"], scale. = T)
# Subtype <- pca1[, "Grade"] %>% as.factor()
# 
# 
# library(ggbiplot)
# g <- ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
#          groups = Subtype, ellipse = T,
#          circle = T,
#          var.axes = F,
#          choices = c(1, 2)) +
#   #scale_color_manual(values = cbcols) +
#   theme_bw() + 
#   theme(legend.direction = 'horizontal', 
#                legend.position = 'top')
# ggsave("PCA of Immune Cell Abundance (Density per Mpix).pdf" , plot = g, device = "pdf",
#        path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Cell_Density/",
#        height = 6, width = 6, units = 'in', dpi = 600)





# Plot
cbcols <- c("#999999",
            "#56B4E9",
            "#E69F00",
            "#009E73")

data_clin3 <- data_clin2 %>% gather(contains("Positive Cells in the"), key = "Parameter", value = "Average.cell.Density")
data_clin4 <- factorthese(data_clin3, c("Grade", "HER2", "ER", "Path.Response", "Parameter"))


## ER
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
ER <- list()
for(i in levels(data_clin4$Parameter)){
  print(i)
  work <- droplevels(subset(data_clin4, Parameter == i))
  x <- compare_means(Average.cell.Density ~ ER, data = work, method = "wilcox.test")
  ER[[i]] <- x
}
ER1 = do.call(rbind, ER)
ER1$id <- rep(names(ER), sapply(ER, nrow))

write.csv(ER1, file = "./Figures/Cell_Density/Stats/ER_stats.csv", row.names = F)

## HER2
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
HER2 <- list()
for(i in levels(data_clin4$Parameter)){
  print(i)
  work <- droplevels(subset(data_clin4, Parameter == i))
  x <- compare_means(Average.cell.Density ~ HER2, data = work, method = "wilcox.test")
  HER2[[i]] <- x
}
HER2a = do.call(rbind, HER2)

HER2a$id <- rep(names(HER2), sapply(HER2, nrow))
write.csv(HER2a, file = "./Figures/Cell_Density/Stats/HER2_stats.csv", row.names = F)


## Grade
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

grade <- list()
for(i in levels(data_clin4$Parameter)){
  print(i)
  work <- droplevels(subset(data_clin4, Parameter == i))
  x <- compare_means(Average.cell.Density ~ Grade, data = work, method = "wilcox.test")
  grade[[i]] <- x
}
grade1 = do.call(rbind, grade)
grade1$id <- rep(names(grade), sapply(grade, nrow))

write.csv(grade1, file = "./Figures/Cell_Density/Stats/grade_stats.csv", row.names = F)


## Path Response
levels(data_clin4$Path.Response)

for(i in levels(data_clin4$Parameter)){
  print(i)
  df <- droplevels(subset(data_clin4, Parameter == i))
  p <- ggplot(df, aes(x = Path.Response, y = Average.cell.Density)) +
    geom_boxplot(alpha = 0.5, width = 0.2) +
    geom_violin(aes(Path.Response, fill = Path.Response), scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) + labs(x = "Path.Response", y = i) +
    geom_dotplot(binaxis = "y",
                 method = "histodot",
                 stackdir = "center",
                 binwidth = 20,
                 position = position_jitter(0.1),
                 alpha = 0,
                 dotsize = 0.4)+
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    stat_compare_means(comparisons = list(c("no", "PCR"),
                                          c("PCR", "PPR"),
                                          c("no", "PPR")), label = "p.signif")
  filen <- paste0(i,".pdf")
  ggsave(filen, plot = p, device = "pdf",
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Cell_Density/Path.Response",
         height = 5, width = 5, units = "in")
}

PResp <- list()
for(i in levels(data_clin4$Parameter)){
  print(i)
  work <- droplevels(subset(data_clin4, Parameter == i))
  x <- compare_means(Average.cell.Density ~ Path.Response, data = work, method = "wilcox.test")
  PResp[[i]] <- x
}
PResp1 = do.call(rbind, PResp)
PResp1$id <- rep(names(PResp), sapply(PResp, nrow))

write.csv(PResp1, file = "./Figures/Cell_Density/Stats/PResp_stats.csv", row.names = F)

# ------------------------------------------------------------------------------------------

# Survival analysis
library(survival)
library(survminer)

data_clin4$OS_STATUS <- ifelse((data_clin4$Dead.vs.alive == "alive"), T, F)

for(i in levels(data_clin4$Parameter)){
  print(i)
  setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Cell_Density/Survival/")
  working <- droplevels(subset(data_clin4, Parameter == i))
  med <- median(working$Average.cell.Density)
  working$Category <- ifelse((working$Average.cell.Density >= med), "High", "Low")
  os.surv.fit <- survfit(Surv(survival.in.months, OS_STATUS) ~ Category,
                         data = working)
  filen <- paste0(i, ".pdf")
  pdf(filen, height = 6, width = 6)
  print(ggsurvplot(os.surv.fit, data = working, risk.table = F, pval = T))
  dev.off()
}

# ------------------------------------------------------------------------------------------





# TSPAN
TSPAN <- read.csv("./Data/TSPAN6_July.csv")
colnames(TSPAN)[colnames(TSPAN) == "Number"] <- "Slide.ID"
TSPAN$Slide.ID <- trim.trailing(TSPAN$Slide.ID)

TSPAN <- TSPAN[, c("Slide.ID", "Tetraspanin.6..membrane.",
                   "Tetraspanin.6..cytoplasm.",
                   "TSpan6.Pos.vs.Neg..memb.", "TSpan6.POS.vs.Neg..Cyt.")]

data_clin <- merge(CS1d, TSPAN, by = "Slide.ID") %>% droplevels()
data_clin$Slide.ID <- as.factor(data_clin$Slide.ID)

# levels(clinical$Slide.ID)[!('%in%'(levels(clinical$Slide.ID), levels(CS1d$Slide.ID)))]
# levels(CS1d$Slide.ID)[grep("S073561", levels(CS1d$Slide.ID))]
head(data_clin)
colnames(data_clin)

data_clin1 <- data_clin[, c("Slide.ID", "Average.cell.Density", "Parameter",
                            "Tetraspanin.6..membrane.",  "Tetraspanin.6..cytoplasm.", "TSpan6.Pos.vs.Neg..memb.",
                            "TSpan6.POS.vs.Neg..Cyt.")]
# Plot
cbcols <- c("#999999",
            "#56B4E9",
            "#E69F00",
            "#009E73")

data_clin4 <- factorthese(data_clin1, c("Parameter",
                                        "Tetraspanin.6..membrane.",  "Tetraspanin.6..cytoplasm.", "TSpan6.Pos.vs.Neg..memb.",
                                        "TSpan6.POS.vs.Neg..Cyt."))

data_TSPAN <- data_clin4
data_TSPAN$TSpan6.Pos.vs.Neg..memb. <- tolower(data_TSPAN$TSpan6.Pos.vs.Neg..memb.)
data_TSPAN$TSpan6.POS.vs.Neg..Cyt. <- tolower(data_TSPAN$TSpan6.POS.vs.Neg..Cyt.)

### Membrane
for(i in levels(data_TSPAN$Parameter)){
  print(i)
  df <- droplevels(subset(data_TSPAN, Parameter == i))
  p <- ggplot(df, aes(x = TSpan6.Pos.vs.Neg..memb., y = Average.cell.Density)) +
    geom_boxplot(alpha = 0.5, width = 0.2) +
    geom_violin(aes(TSpan6.Pos.vs.Neg..memb., fill = TSpan6.Pos.vs.Neg..memb.), scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) + labs(x = "TSpan6.Pos.vs.Neg..memb.", y = i) +
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
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Cell_Density/TSPAN_mem/pos_neg/",
         height = 5, width = 5, units = "in")
}

TSPAN <- list()
for(i in levels(data_TSPAN$Parameter)){
  print(i)
  work <- droplevels(subset(data_TSPAN, Parameter == i))
  x <- compare_means(Average.cell.Density ~ TSpan6.Pos.vs.Neg..memb., data = work, method = "wilcox.test")
  TSPAN[[i]] <- x
}
TSPAN1 = do.call(rbind, TSPAN)
TSPAN1$id <- rep(names(TSPAN), sapply(TSPAN, nrow))

write.csv(TSPAN1, file = "./Figures/Cell_Density/Stats/TSPAN_membrane_PN_stats.csv", row.names = F)

#### Score
data_TSPAN$Tetraspanin.6..membrane. <- as.factor(data_TSPAN$Tetraspanin.6..membrane.)

for(i in levels(data_TSPAN$Parameter)){
  print(i)
  df <- droplevels(subset(data_TSPAN, Parameter == i))
  p <- ggplot(df, aes(x = Tetraspanin.6..membrane., y = Average.cell.Density)) +
    geom_boxplot(alpha = 0.5, width = 0.2) +
    geom_violin(aes(Tetraspanin.6..membrane., fill = Tetraspanin.6..membrane.), scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) + labs(x = "Tetraspanin.6..membrane.", y = i) +
    geom_dotplot(binaxis = "y",
                 method = "histodot",
                 stackdir = "center",
                 binwidth = 20,
                 position = position_jitter(0.1),
                 alpha = 0,
                 dotsize = 0.4)+
    theme_bw()+
    theme(axis.text = element_text(size = 16)) +
    stat_compare_means(comparisons = list(c("0", "1"), c("0", "2"), c("0", "3"),
                                          c("1", "2"), c("1", "3"), c("2", "3")), label = "p.signif")
  filen <- paste0(i,".pdf")
  ggsave(filen, plot = p, device = "pdf",
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Cell_Density/TSPAN_mem/score/",
         height = 5, width = 5, units = "in")
}

TSPAN <- list()
for(i in levels(data_TSPAN$Parameter)){
  print(i)
  work <- droplevels(subset(data_TSPAN, Parameter == i))
  x <- compare_means(Average.cell.Density ~ Tetraspanin.6..membrane., data = work, method = "wilcox.test")
  TSPAN[[i]] <- x
}
TSPAN1 = do.call(rbind, TSPAN)
TSPAN1$id <- rep(names(TSPAN), sapply(TSPAN, nrow))

write.csv(TSPAN1, file = "./Figures/Cell_Density/Stats/TSPAN_membrane_score_stats.csv", row.names = F)

### Cytoplasm
for(i in levels(data_TSPAN$Parameter)){
  print(i)
  df <- droplevels(subset(data_TSPAN, Parameter == i))
  p <- ggplot(df, aes(x = TSpan6.POS.vs.Neg..Cyt., y = Average.cell.Density)) +
    geom_boxplot(alpha = 0.5, width = 0.2) +
    geom_violin(aes(TSpan6.POS.vs.Neg..Cyt., fill = TSpan6.POS.vs.Neg..Cyt.), scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) + labs(x = "TSpan6.POS.vs.Neg..Cyt.", y = i) +
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
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Cell_Density/TSPAN_cyto/pos_neg/",
         height = 5, width = 5, units = "in")
}


TSPAN <- list()
for(i in levels(data_TSPAN$Parameter)){
  print(i)
  work <- droplevels(subset(data_TSPAN, Parameter == i))
  x <- compare_means(Average.cell.Density ~ TSpan6.POS.vs.Neg..Cyt., data = work, method = "wilcox.test")
  TSPAN[[i]] <- x
}
TSPAN1 = do.call(rbind, TSPAN)
TSPAN1$id <- rep(names(TSPAN), sapply(TSPAN, nrow))

write.csv(TSPAN1, file = "./Figures/Cell_Density/Stats/TSPAN_cytoplasm_PN_stats.csv", row.names = F)

#### Score
data_TSPAN$Tetraspanin.6..cytoplasm. <- as.factor(data_TSPAN$Tetraspanin.6..cytoplasm.)

for(i in levels(data_TSPAN$Parameter)){
  print(i)
  df <- droplevels(subset(data_TSPAN, Parameter == i))
  p <- ggplot(df, aes(x = Tetraspanin.6..cytoplasm., y = Average.cell.Density)) +
    geom_boxplot(alpha = 0.5, width = 0.2) +
    geom_violin(aes(Tetraspanin.6..cytoplasm., fill = Tetraspanin.6..cytoplasm.), scale = "width", alpha = 0.8) +
    scale_fill_manual(values = cbcols) + labs(x = "Tetraspanin.6..cytoplasm.", y = i) +
    geom_dotplot(binaxis = "y",
                 method = "histodot",
                 stackdir = "center",
                 binwidth = 20,
                 position = position_jitter(0.1),
                 alpha = 0,
                 dotsize = 0.4)+
    theme_bw()+
    theme(axis.text = element_text(size = 16)) +
    stat_compare_means(comparisons = list(c("0", "1"), c("0", "2"), c("0", "3"),
                                          c("1", "2"), c("1", "3"), c("2", "3")), label = "p.signif")
  filen <- paste0(i,".pdf")
  ggsave(filen, plot = p, device = "pdf",
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Cell_Density/TSPAN_cyto/score/",
         height = 5, width = 5, units = "in")
}

TSPAN <- list()
for(i in levels(data_TSPAN$Parameter)){
  print(i)
  work <- droplevels(subset(data_TSPAN, Parameter == i))
  x <- compare_means(Average.cell.Density ~ Tetraspanin.6..cytoplasm., data = work, method = "wilcox.test")
  TSPAN[[i]] <- x
}
TSPAN1 = do.call(rbind, TSPAN)
TSPAN1$id <- rep(names(TSPAN), sapply(TSPAN, nrow))

write.csv(TSPAN1, file = "./Figures/Cell_Density/Stats/TSPAN_cytoplasm_score_stats.csv", row.names = F)









