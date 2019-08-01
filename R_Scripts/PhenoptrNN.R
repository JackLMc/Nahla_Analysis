# Source Functions + Packages
library(UsefulFunctions)
library(tidyverse)
library(data.table)
library(phenoptr)
# Search for NN9 if already previously ran

# NEAREST NEIGHBOUR
## Convert old DF to a tsv file
### Batch 1
csvtotsv("./Output/Batch1_single.csv",
         "./Output/NN/Batch1_single.tsv")

#### Calculate distances
NearNeigh <- compute_all_nearest_distance("./Output/NN/Batch1_single.tsv")
write.csv(file = "./Output/Batch1_NN.csv", x = NearNeigh, row.names = F)

### Batch 2
csvtotsv("./Output/Batch2_single.csv",
         "./Output/NN/Batch2_single.tsv")

#### Calculate distances
NearNeigh1 <- compute_all_nearest_distance("./Output/NN/Batch2_single.tsv")
write.csv(file = "./Output/Batch2_NN.csv", x = NearNeigh1, row.names = F)

# Write a TSV or a CSV?
# Actually have to read this in, else it's as a tibble
Batch1_NN <- read.csv("Output/Batch1_NN.csv")
Batch2_NN <- read.csv("Output/Batch2_NN.csv")

rm(NearNeigh)
rm(NearNeigh1)

# Remove stuff that doesn't matter anymore
Batch1_NN$tag <- as.character(Batch1_NN$tag)
Batch2_NN$tag <- as.character(Batch2_NN$tag)
NN <- rbind(Batch1_NN, Batch2_NN)

write.csv(file = "./Output/Nearest_Neighbour_Start.csv", x = NN, row.names = F)
drop <- c("tag",
          "Cell.ID",
          "Cell.X.Position",
          "Cell.Y.Position",
          "Category.Region.ID",
          "Confidence",
          "Batch",
          "Distance.from.Tissue.Category.Edge..microns")
NN1 <- droplevels(NN[,!(names(NN) %in% drop)])
NN1$Slide.ID <- as.factor(NN1$Slide.ID)
NN1$Sample.Name <- as.factor(NN1$Sample.Name)
nlevels(NN1$Slide.ID)

# Gather Distance variables into column to
NN2 <- NN1 %>% gather(contains("Distance.to"), key = "PhenoTo", value = "Distance")

# Make a unique column 
NN2$Comb <- as.factor(paste(NN2$Tissue.Category, NN2$Phenotype, NN2$PhenoTo, NN2$Sample.Name, sep = "/"))

# Make two columns of pheno from and pheno to
## Drop some stuff to allow splitting of Comb
drop <- c("Sample.Name",
          "Tissue.Category",
          "Phenotype",
          "PhenoTo")
NN3 <- droplevels(NN2[,!(names(NN2) %in% drop)])

## Separate 
head(NN3)
NN4 <- NN3 %>%
  separate(Comb, c("Tissue.Category", "PhenoFrom", "PhenoTo", "Sample.Name"), "/")

## Make the PhenoTo and PhenoFrom the same
NN4$PhenoTo <- gsub("Distance.to.", "", NN4$PhenoTo)
NN4$PhenoTo <- gsub(" ", ".", NN4$PhenoTo)
NN4$PhenoFrom <- gsub(" ", ".", NN4$PhenoFrom)

## Remove rows where either PhenoTo or PhenoFrom is DAPI
NN4$PhenoTo <- as.factor(NN4$PhenoTo)
NN4$PhenoFrom <- as.factor(NN4$PhenoFrom)
NN6 <- NN4

# # Paste back together
# NN5$Comb <- as.factor(paste(NN5$Tissue.Category, NN5$PhenoFrom, ".Distance.To.", NN5$PhenoTo, NN5$Sample.Name, sep = "/"))
# 
# # Factor them all
# NN5$Comb <- as.factor(NN5$Comb)
# NN5$Slide.ID <- as.factor(NN5$Slide.ID)
# 
# ## Loop and find the mean distance (for the replicates to get a biological replicate)
# output <- data.frame(Comb = character(),
#                      MeanDist = double(),
#                      Slide.ID = character(),
#                      stringsAsFactors = FALSE)
# 
# c <- 1
# for(i in levels(NN5$Comb)){
#   name <- basename(i)
#   cat("Processing", i, "\n")
#   df <- droplevels(subset(NN5, Comb == i))
#   y <- mean(df$Distance, na.rm = T)
#   output[c, "Comb"] <- i
#   output[c, "MeanDist"] <- y
#   output[c, "Slide.ID"] <- as.character(levels(df$Slide.ID))
#   c <- c + 1}
# 
# # Separate back into original columns
# NN6 <- output %>%
#   separate(Comb, c("Tissue.Category", "PhenoFrom","Distance.To", "PhenoTo", "Sample.Name"), "/")

# Combine and reloop for mean nearest distance on a slide basis
head(NN6)
NN6$Parameter <- as.factor(paste(NN6$Tissue.Category, NN6$PhenoFrom, ".Distance.To.", NN6$PhenoTo, NN6$Slide.ID, sep = "/"))
NN6a <- data.frame(Parameter = character(),
                     Distance = double(),
                     stringsAsFactors = F)

c <- 1
for(i in levels(NN6$Parameter)){
  name <- basename(i)
  cat("Processing", i, "\n")
  df <- droplevels(subset(NN6, Parameter == i))
  y <- mean(df$Distance, na.rm = T)
  NN6a[c, "Parameter"] <- i
  NN6a[c, "Distance"] <- y
  c <- c + 1}

# Remove Slide.ID from combined Parameter column
NN7 <- NN6a %>%
  separate(Parameter, c("Tissue.Category", "PhenoFrom","Distance.To", "PhenoTo", "Slide.ID"), "/")
# cell_seg_path <- c("/Users/JackMcMurray/Desktop/")
# Phenotypes <- c("CD4", "CD8", "CD20", "CD68", "FOXP3", "DAPI")
# pairs <- find_pheno_comb(Phenotypes)
# 
# out_path <- path.expand("~/spatial_distribution_report.html")
# 
# spatial_distribution_report(cell_seg_path, pairs, output_path=out_path)
# base_path <- "/Users/JackMcMurray/Desktop/8-B 10-5866/"
# setwd("/Users/JackMcMurray/Desktop/8-B 10-5866/")
# colors <- c("CD4"="#999999",
#             "CD8"="#E69F00",
#             "CD20" = "#56B4E9",
#             "CD68" = "#009E73",
#             "FOXP3" = "#F0E442",
#             "DAPI" = "#0072B2")
# 
# 
# paths <- list_cell_seg_files(base_path)
# for (path in paths)
#   spatial_distribution_report(path, pairs, colors)

NN7$Parameter <- as.factor(paste(NN7$Tissue.Category, NN7$PhenoFrom, NN7$Distance.To, NN7$PhenoTo, sep = "_"))

# Drop what you don't need
drop <- c("Tissue.Category",
          "Distance.To",
          "PhenoFrom",
          "PhenoTo")
NN8 <- droplevels(NN7[,!(names(NN7) %in% drop)])
NN8$Slide.ID <- as.factor(NN8$Slide.ID)

# Spread
NN9 <- spread(NN8, key = "Parameter", value = "Distance")
NN9$Slide.ID <- gsub("  5PLEX|5plex Fedor| 5PLEX| breast tumour CD4 CD8 CD20 CD68 FOXP3", "", NN9$Slide.ID) 
NN9$Slide.ID <-trim.trailing(NN9$Slide.ID)

write.csv("./Output/Clean_Nearest_Neighbour.csv", x = NN9, row.names = F)

####### START!
library(UsefulFunctions)
library(tidyverse)
library(data.table)
library(phenoptr)
library(ggpubr)

# Read in data for comparisons
NN9 <- read.csv("Output/Clean_Nearest_Neighbour.csv")
clin <- read.csv("Data/Nahla_clinical.csv")
colnames(clin)[colnames(clin) == "Lablels.as.in.inform.cases"] <- "Slide.ID"
clin$Slide.ID <- trim.trailing(clin$Slide.ID)
clin <- droplevels(subset(clin, Slide.ID != ""))
clin$ER[clin$ER == "Negative"] <- "negative"
clin$HER2[clin$HER2 == "Positive"] <- "positive"
clin$Path.Response[clin$Path.Response == "pCR"] <- "PCR"

# Merge with Clinical data
NN10 <- droplevels(merge(NN9, clin, by = "Slide.ID"))
length(NN10$Slide.ID)
length(NN9$Slide.ID)
length(clin$Slide.ID)


#
library(tidyverse)
NN11 <- NN10 %>% gather(contains("Distance.To"), key = "Parameter", value = "Distance")
NN11$Parameter <- as.factor(NN11$Parameter)
NN12 <- na.omit(NN11)

cbcols <- c("#999999", "#56B4E9", "#E69F00", "#009E73")


for(i in levels(NN12$Parameter)){
  print(i)
  df <- droplevels(subset(NN12, Parameter == i))
  p <- ggplot(df, aes(x = ER, y = Distance)) +
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
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Nearest_Neighbour/ER",
         height = 5, width = 5, units = "in")
}

ER <- list()
for(i in levels(NN12$Parameter)){
  print(i)
  work <- droplevels(subset(NN12, Parameter == i))
  x <- compare_means(Distance ~ ER, data = work, method = "wilcox.test")
  ER[[i]] <- x
}
ER1 = do.call(rbind, ER)
ER1$id <- rep(names(ER), sapply(ER, nrow))

write.csv(ER1, file = "./Figures/Nearest_Neighbour/Stats/ER_stats.csv", row.names = F)


## HER2
for(i in levels(NN12$Parameter)){
  print(i)
  df <- droplevels(subset(NN12, Parameter == i))
  p <- ggplot(df, aes(x = HER2, y = Distance)) +
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
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Nearest_Neighbour/HER2",
         height = 5, width = 5, units = "in")
}

HER2 <- list()
for(i in levels(NN12$Parameter)){
  print(i)
  work <- droplevels(subset(NN12, Parameter == i))
  x <- compare_means(Distance ~ HER2, data = work, method = "wilcox.test")
  HER2[[i]] <- x
}
HER2a = do.call(rbind, HER2)

HER2a$id <- rep(names(HER2), sapply(HER2, nrow))
write.csv(HER2a, file = "./Figures/Nearest_Neighbour/Stats/HER2_stats.csv", row.names = F)


## grade
NN12$Grade <- as.factor(NN12$Grade)

for(i in levels(NN12$Parameter)){
  print(i)
  df <- droplevels(subset(NN12, Parameter == i))
  p <- ggplot(df, aes(x = Grade, y = Distance)) +
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
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Nearest_Neighbour/Grade",
         height = 5, width = 5, units = "in")
}

grade <- list()
for(i in levels(NN12$Parameter)){
  print(i)
  work <- droplevels(subset(NN12, Parameter == i))
  x <- compare_means(Distance ~ Grade, data = work, method = "wilcox.test")
  grade[[i]] <- x
}
grade1 = do.call(rbind, grade)
grade1$id <- rep(names(grade), sapply(grade, nrow))

write.csv(grade1, file = "./Figures/Nearest_Neighbour/Stats/grade_stats.csv", row.names = F)


## Path Response
levels(NN12$Path.Response)

for(i in levels(NN12$Parameter)){
  print(i)
  df <- droplevels(subset(NN12, Parameter == i))
  p <- ggplot(df, aes(x = Path.Response, y = Distance)) +
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
    theme_bw()+
    theme(axis.text = element_text(size = 16)) +
    stat_compare_means(comparisons = list(c("no", "PCR"),
                                          c("PCR", "PPR"),
                                          c("no", "PPR")), label = "p.signif")
  filen <- paste0(i,".pdf")
  ggsave(filen, plot = p, device = "pdf",
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Nearest_Neighbour/Path.Response",
         height = 5, width = 5, units = "in")
}

PResp <- list()
for(i in levels(NN12$Parameter)){
  print(i)
  work <- droplevels(subset(NN12, Parameter == i))
  x <- compare_means(Distance ~ Path.Response, data = work, method = "wilcox.test")
  PResp[[i]] <- x
}
PResp1 = do.call(rbind, PResp)
PResp1$id <- rep(names(PResp), sapply(PResp, nrow))

write.csv(PResp1, file = "./Figures/Nearest_Neighbour/Stats/PResp_stats.csv", row.names = F)

# ------------------------------------------------------------------------------------------

# Survival analysis
library(survival)
library(survminer)

NN12$OS_STATUS <- ifelse((NN12$Dead == "alive"), T, F)

for(i in levels(NN12$Parameter)){
  print(i)
  setwd("/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Nearest_Neighbour/Survival/")
  working <- droplevels(subset(NN12, Parameter == i))
  med <- median(working$Distance)
  working$Category <- ifelse((working$Distance >= med), "High", "Low")
  os.surv.fit <- survfit(Surv(survival.in.months, OS_STATUS) ~ Category,
                         data = working)
  filen <- paste0(i, ".pdf")
  pdf(filen, height = 6, width = 6)
  print(ggsurvplot(os.surv.fit, data = working, risk.table = F, pval = T))
  dev.off()
}

# ------------------------------------------------------------------------------------------


# TSPAN analysis
TSPAN <- read.csv("Data/TSPAN6_July.csv")
colnames(TSPAN)[colnames(TSPAN) == "Number"] <- "Slide.ID"
TSPAN$Slide.ID <- trim.trailing(TSPAN$Slide.ID)

TSPAN <- TSPAN[, c("Slide.ID", "Tetraspanin.6..membrane.",
                   "Tetraspanin.6..cytoplasm.",
                   "TSpan6.Pos.vs.Neg..memb.", "TSpan6.POS.vs.Neg..Cyt.")]

NN_Ts <- droplevels(merge(NN9, TSPAN, by = "Slide.ID"))
length(NN_Ts$Slide.ID)
length(NN9$Slide.ID)

#
library(tidyverse)
NN_Ts1 <- NN_Ts %>% gather(contains("Distance.To"), key = "Parameter", value = "Distance")
head(NN_Ts1)

NN_Ts1$Parameter <- as.factor(NN_Ts1$Parameter)

head(NN_Ts1)
# NN_impute <- rfImpute(ER + HER2 + Path.Response ~ ., NN_Ts1)
# View(NN_Ts1)

NN_Ts2 <- na.omit(NN_Ts1)

NN_Ts2 <- factorthese(NN_Ts2, c("Tetraspanin.6..membrane.",
                                "Tetraspanin.6..cytoplasm.",
                                "TSpan6.Pos.vs.Neg..memb.", "TSpan6.POS.vs.Neg..Cyt."))


NN_Ts2$TSpan6.Pos.vs.Neg..memb. <- tolower(NN_Ts2$TSpan6.Pos.vs.Neg..memb.)
NN_Ts2$TSpan6.POS.vs.Neg..Cyt. <- tolower(NN_Ts2$TSpan6.POS.vs.Neg..Cyt.)


levels(NN_Ts2$Tetraspanin.6..membrane.)
levels(NN_Ts2$Tetraspanin.6..cytoplasm.)
levels(NN_Ts2$TSpan6.Pos.vs.Neg..memb.)
levels(NN_Ts2$TSpan6.POS.vs.Neg..Cyt.)


## TSPAN
### Membrane
for(i in levels(NN_Ts2$Parameter)){
  print(i)
  df <- droplevels(subset(NN_Ts2, Parameter == i))
  p <- ggplot(df, aes(x = TSpan6.Pos.vs.Neg..memb., y = Distance)) +
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
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Nearest_Neighbour/TSPAN_mem/pos_neg/",
         height = 5, width = 5, units = "in")
}

TSPAN <- list()
for(i in levels(NN_Ts2$Parameter)){
  print(i)
  work <- droplevels(subset(NN_Ts2, Parameter == i))
  x <- compare_means(Distance ~ TSpan6.Pos.vs.Neg..memb., data = work, method = "wilcox.test")
  TSPAN[[i]] <- x
}
TSPAN1 = do.call(rbind, TSPAN)
TSPAN1$id <- rep(names(TSPAN), sapply(TSPAN, nrow))

write.csv(TSPAN1, file = "./Figures/Nearest_Neighbour/Stats/TSPAN_membrane_PN_stats.csv", row.names = F)

#### Score
my_comparisons <- list(c("0", "1"),
                       c("0", "2"),
                       c("0", "3"),
                       c("1", "2"),
                       c("1", "3"),
                       c("2", "3"))

NN_Ts2$Tetraspanin.6..membrane. <- as.factor(NN_Ts2$Tetraspanin.6..membrane.)

for(i in levels(NN_Ts2$Parameter)){
  print(i)
  df <- droplevels(subset(NN_Ts2, Parameter == i))
  p <- ggplot(df, aes(x = Tetraspanin.6..membrane., y = Distance)) +
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
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Nearest_Neighbour/TSPAN_mem/score/",
         height = 5, width = 5, units = "in")
}

TSPAN <- list()
for(i in levels(NN_Ts2$Parameter)){
  print(i)
  work <- droplevels(subset(NN_Ts2, Parameter == i))
  x <- compare_means(Distance ~ Tetraspanin.6..membrane., data = work, method = "wilcox.test")
  TSPAN[[i]] <- x
}
TSPAN1 = do.call(rbind, TSPAN)
TSPAN1$id <- rep(names(TSPAN), sapply(TSPAN, nrow))

write.csv(TSPAN1, file = "./Figures/Nearest_Neighbour/Stats/TSPAN_membrane_score_stats.csv", row.names = F)

### Cytoplasm
for(i in levels(NN_Ts2$Parameter)){
  print(i)
  df <- droplevels(subset(NN_Ts2, Parameter == i))
  p <- ggplot(df, aes(x = TSpan6.POS.vs.Neg..Cyt., y = Distance)) +
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
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Nearest_Neighbour/TSPAN_cyto/pos_neg/",
         height = 5, width = 5, units = "in")
}


TSPAN <- list()
for(i in levels(NN_Ts2$Parameter)){
  print(i)
  work <- droplevels(subset(NN_Ts2, Parameter == i))
  x <- compare_means(Distance ~ TSpan6.POS.vs.Neg..Cyt., data = work, method = "wilcox.test")
  TSPAN[[i]] <- x
}
TSPAN1 = do.call(rbind, TSPAN)
TSPAN1$id <- rep(names(TSPAN), sapply(TSPAN, nrow))

write.csv(TSPAN1, file = "./Figures/Nearest_Neighbour/Stats/TSPAN_cytoplasm_PN_stats.csv", row.names = F)

#### Score
NN_Ts2$Tetraspanin.6..cytoplasm. <- as.factor(NN_Ts2$Tetraspanin.6..cytoplasm.)

for(i in levels(NN_Ts2$Parameter)){
  print(i)
  df <- droplevels(subset(NN_Ts2, Parameter == i))
  p <- ggplot(df, aes(x = Tetraspanin.6..cytoplasm., y = Distance)) +
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
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Nearest_Neighbour/TSPAN_cyto/score/",
         height = 5, width = 5, units = "in")
}

TSPAN <- list()
for(i in levels(NN_Ts2$Parameter)){
  print(i)
  work <- droplevels(subset(NN_Ts2, Parameter == i))
  x <- compare_means(Distance ~ Tetraspanin.6..cytoplasm., data = work, method = "wilcox.test")
  TSPAN[[i]] <- x
}
TSPAN1 = do.call(rbind, TSPAN)
TSPAN1$id <- rep(names(TSPAN), sapply(TSPAN, nrow))

write.csv(TSPAN1, file = "./Figures/Nearest_Neighbour/Stats/TSPAN_cytoplasm_score_stats.csv", row.names = F)











#
NN9a <- data.frame(NN10[, names(NN10) != "Slide.ID"], row.names = NN10[, names(NN10) == "Slide.ID"])
library(dplyr)
NN11 <- NN9a %>% dplyr:: select(-contains("DAPI")) %>% dplyr:: select(-contains("Slide.ID")) 
head(NN11)

NN_impute <- na.omit(NN10)
# randomForest imputing
library(randomForest)
NN_impute <- rfImpute(Subtype ~ ., NN10)

# writeCsvO(NN_impute)

# Logistic Regression
corrplot(cor(NN_impute[, names(NN_impute) != "Subtype"]), type = "lower")

NN_impute$Subtype <- as.factor(NN_impute$Subtype)
library(arm)
model <- bayesglm(Subtype ~., family = binomial (link = "logit"), data = NN_impute)
display(model)





# PCA
### Complete PCA
head(NN11)

prin_comp <- prcomp(NN11, scale. = T)
plot(prin_comp)

Subtype <- NN_impute[, "Subtype"]
g <- ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, ellipse = T,
              circle = F,
              var.axes = F
)
g <- g + scale_color_manual(values = cbcols)
g <- g + theme_bw()
g <- g + theme(legend.direction = "horizontal", 
               legend.position = "top")
g <- g + ggtitle("PCA of the Mean Nearest Distance Between Immune Subsets")
g
ggsave("PCA of Nearest Neighbour.png" , plot = g, device = "png",
       path = "/Users/JackMcMurray/OneDrive/University_of_Birmingham/PhD/Extra/Nahla_Analysis/Figures/NearestNeighbour",
       height = 6, width = 6, units = "in", dpi = 600)


# Store the contribution to each PC
PC_NN <- as.data.frame(prin_comp$x)
PC_NN1 <- tibble:: rownames_to_column(PC_NN, "Slide.ID")
PC_NN1$Slide.ID <- as.factor(PC_NN1$Slide.ID)

# Change the Types to proper types
PC_NN1[grep("dMMR$", PC_NN1$Slide.ID), "Subtype"] <- "MSI"
PC_NN1[grep("pMMR$", PC_NN1$Slide.ID), "Subtype"] <- "MSS"
PC_NN1[grep(".*N[ ]S[[:digit:]]{6}$", PC_NN1$Slide.ID), "Subtype"] <- "hiCIRC"
PC_NN1$Subtype <- as.factor(PC_NN1$Subtype)

PC_NN3 <- PC_NN1 %>% gather(contains("PC"), key = Component, value = "ComponentScore")
PC_NN3$Component <- as.factor(PC_NN3$Component)
PC_NN3$Rank <- rank(PC_NN3$ComponentScore)

# Plot
for(i in levels(PC_NN3$Component)){
  name <- basename(i)
  cat("Processing", i, "\n")
  Chosen <- droplevels(subset(PC_NN3, Component == i))
  YTitle <- paste0("Rank of ", i)
  MainTitle <- paste0("Violin plot of ", i)
  temp_plot <- ggSubtype(Chosen, YTitle, MainTitle)
  filen <- paste0(i,".png")
  ggsave(filen, plot = temp_plot, device = "png",
         path = "/Users/JackMcMurray/OneDrive/University_of_Birmingham/PhD/Extra/Nahla_Analysis/Figures/NearestNeighbour",
         height=5, width=5, units="in", dpi=600)
}

#
## Find the Eigenvalues
fviz_screeplot(prin_comp, ncp = 10, choice = "eigenvalue")
eig <- (prin_comp$sdev)^2

## Variances in percentage
variance <- eig*100/sum(eig)

## Cumulative variances
cumvar <- cumsum(variance)

Eigenvalues.All <- data.frame(eig = eig, variance = variance,
                              cumvariance = cumvar)
# writeCsvO(Eigenvalues.All)

# Store the variances
var <- get_pca_var(prin_comp)

## Find the correlation between variables and principal components
loadings <- prin_comp$rotation
loads <- as.data.frame(loadings)
load1 <- tibble:: rownames_to_column(loads, "Parameter")
load1$Parameter <- as.factor(load1$Parameter)

# Biggest contributions to a given PC
# load1[order(load1$PC1, decreasing = T)[1:4],]
# droplevels(subset(load1, Parameter == "CD8.10.Class.2.Epithelium"))
# levels(load1$Parameter)

sdev <- prin_comp$sdev

### Find the correlataions
var.coord <- t(apply(loadings, 1, var_cor_func, sdev))

## Calculate the Cos2 (square of the coordinates)
var.cos2 <- var.coord^2

## Contribution of variables to the Components ((cos2) * 100) / (total cos2 of the component)
comp.cos2 <- apply(var.cos2, 2, sum)
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
contrib.var <- as.data.frame(var.contrib)

## Find the most contributing variable
contrib.var[order(contrib.var$PC1, decreasing = T)[1:10],]

# LDA
prop.pca <- prin_comp$sdev^2/sum(prin_comp$sdev^2)

r <- lda(Subtype ~ ., 
         NN_impute
         ,prior = c(0.3296089, 0.3351955, 0.3351955)
         )

prop.lda = r$svd^2/sum(r$svd^2)

plda <- predict(object = r,
                newdata = NN6)

dataset <- data.frame(Type = NN6[,"Type"],
                      pca = prin_comp$x, lda = plda$x)


p1 <- ggplot(dataset) + geom_point(aes(lda.LD1, lda.LD2, colour = Type), size = 2.5) + 
  labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
       y = paste("LD2 (", percent(prop.lda[2]), ")", sep="")) + theme_bw()


# tSNE
tsne <- Rtsne(as.matrix(NN_impute[, names(NN_impute) != "Subtype"]), dims = 2, perplexity = 8, verbose = T, max_iter = 10000)
plot(tsne$Y, main = "tsne", col = NN_impute$Subtype, pch = 16,
     axes = F)



# Look at interactions
head(NN_impute)

# Remake slide.ID column
NNI1 <- tibble:: rownames_to_column(NN_impute, var = "Slide.ID")
head(NNI1)
NNI1[grep("dMMR$", NNI1$Slide.ID), "Subtype"] <- "MSI"
NNI1[grep("pMMR$", NNI1$Slide.ID), "Subtype"] <- "MSS"
NNI1[grep(".*N[ ]S[[:digit:]]{6}$", NNI1$Slide.ID), "Subtype"] <- "hiCIRC"

NNI2 <- NNI1 %>% gather(contains("_.Distance.To._"), key = "param", value = "Min_Dist")

NNI2$Parameter <- as.factor(paste(NNI2$Subtype, NNI2$param, sep = "_"))

head(NNI2)

levels(NNI2$Parameter)

attempt <- droplevels(subset(NNI2, Parameter == "MSS_Stroma_CD4...T.bet_.Distance.To._CD4"))

mean(attempt$Min_Dist)

output <- data.frame(Parameter = character(),
                     Min_Dist = double(),
                     stringsAsFactors = F)

c <- 1

for(i in levels(NNI2$Parameter)){
  work <- droplevels(subset(NNI2, Parameter == i))
  meanDist <- mean(work$Min_Dist)
  output[c, "Parameter"] <- i
  output[c, "Min_Dist"] <- meanDist
  c <- c + 1
}

NNI3 <- output %>% separate(Parameter, c("Subtype", "Tissue.Category", "PhenoFrom", "Distance.To", "PhenoTo"), sep = "_")
NNI4 <- NNI3

NNI4$Graph <- as.factor(paste(NNI3$Subtype, NNI3$Tissue.Category, NNI3$PhenoTo, sep = "_"))
drop <- c("Subtype",
          "Tissue.Category",
          "PhenoTo",
          "Distance.To")
NNI5 <- droplevels(NNI4[,!(names(NNI4) %in% drop)])
NNI6 <- as.tibble(spread(NNI5, key = "Graph", value = "Min_Dist")) 
ggplot(NNI6, aes(`hiCIRC_Epithelium_CD4`, color=PhenoFrom)) + 
  geom_density(size = 1) + 
  theme_minimal()






head(NN12)
NN12a <- NN12[, c("Slide.ID", "Parameter", "Distance", "grade", "ER", "HER2", "Path.Response")]
NN13 <- spread(NN12a, key = "Parameter", value = "Distance")
pca1 <- column_to_rownames(NN13, var = "Slide.ID")





prin_comp <- prcomp(pca1[, names(pca1) != "grade" & names(pca1) != "HER2" & names(pca1) != "ER" & names(pca1) != "Path.Response"], scale. = T)
Subtype <- pca1[, "grade"]

library(ggbiplot)
g <- ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
              groups = Subtype, ellipse = T,
              circle = T,
              var.axes = F,
              choices = c(1,2)
)


g <- g + scale_color_manual(values = cbcols)
g <- g + theme_bw()
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')

