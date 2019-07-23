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

nlevels(clinical$Slide.ID)
nlevels(CS1d$Slide.ID)


try <- merge(CS1d, clinical, by = "Slide.ID") %>% droplevels()
try$Slide.ID <- as.factor(try$Slide.ID)


# levels(clinical$Slide.ID)[!('%in%'(levels(clinical$Slide.ID), levels(CS1d$Slide.ID)))]
# levels(CS1d$Slide.ID)[grep("S073561", levels(CS1d$Slide.ID))]


try1 <- try[, c("Slide.ID", "Average.cell.Density", "Parameter", "grade", "HER2", "ER", "Path.Response")]
try2 <- spread(try1, key = "Parameter", value = "Average.cell.Density") %>% column_to_rownames(., var = "Slide.ID")

head(try2)

pca1 <- try2 %>% dplyr:: select(-contains("All")) %>% dplyr:: select(-contains("DAPI"))
pca1[is.na(pca1)] <- 0

View(try2)


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

g <- g + ggtitle("PCA of Immune Cell Abundance (Density per Megapixel)")
ggsave("PCA of Immune Cell Abundance (Density per Mpix).png" , plot = g, device = "png",
       path = "/Users/JackMcMurray/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC/Figures",
       height = 6, width = 6, units = 'in', dpi = 600)





# Plot
for(i in levels(CS1d$Parameter)){
  print(paste("Working on:", i))
  Chosen <- droplevels(subset(CS1d, Parameter == i))
  YTitle <- paste0("Rank of the ", i)
  MainTitle <- paste0("Violin plot of the ", i)
  temp_plot <- ggSubtype(Chosen, YTitle, MainTitle)
  filen <- paste0(i,".png")
  ggsave(filen, plot = temp_plot, device = "png",
         path = "/Users/JackMcMurray/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Cell_Density",
         height = 5, width = 5, units = 'in', dpi = 600)
}

# Combine stats into a list
stat_list <- list()
c <- 1
for(i in levels(CS1d$Parameter)) {
  name <- basename(i)
  cat('Processing', i, '\n')
  workingon <- droplevels(subset(CS1d, Parameter == i))
  # assign your ggplot call to the i'th position in the list
  x <- compare_means(Average.cell.Density ~ Subtype, data = workingon, method = "wilcox.test")
  y <- as.data.frame(x)
  stat_list[[i]]  <- cbind(i, y)
c <- c + 1}

# Bind and remove row names
z <- do.call(rbind, stat_list)
rownames(z) <- c()

# Separate
CellDensity_Megapixel_Stat <- z 

# Write out the statistics
writeCsvO(CellDensity_Megapixel_Stat)

# Logistic Regression
CS1e <- CS1d %>% dplyr:: select(-contains("Phenotype")) %>% dplyr:: select(-contains("Tissue.Category")) %>% dplyr:: select(-contains("Rank"))
LR <- spread(CS1e, "Parameter", "Average.cell.Density")
LR[is.na(LR)] <- 0
LR2 <- LR %>% dplyr:: select(-contains("Slide.ID"))

# Using model to predict where hiCIRC patients should lie?
LR2 <- droplevels(subset(LR2, Subtype != "hiCIRC"))

### Find what is 0 and what is 1
contrasts(LR2$Subtype)

# Check for collinearity
corrplot(cor(LR2[,-1]), type = 'lower')

# Remove variables that are likely to be collinear
LR3 <- LR2 %>% dplyr:: select(-contains("All")) %>% dplyr:: select(-contains("Other"))

# Recheck collinearity
corrplot(cor(LR3[,-1]), type = 'lower')

# Do model
model <- glm(Subtype ~., family = binomial (link = 'logit'), data = LR3)

# Investigate model
summary(model)
anova(model, test = "Chisq")

## Find R squared
pR2(model)

# Predict where hiCIRC lie (rename hiCIRC - MSS [what the MSI test says])
LR2a <- droplevels(subset(LR, Subtype == "hiCIRC"))
LR2b <- LR2a %>% dplyr:: select(-contains("Subtype"))
LR2b$Subtype <- 1 # MSS binary

# See how accurately the model assigns hiCIRC to MSS
fitted.results <- predict(model, newdata = LR2b, type = 'response')
fitted.results <- ifelse((fitted.results > 0.5), 1, 0)
misClasificError <- mean(fitted.results != LR2b$Subtype)
print(paste("Accuracy", 1-misClasificError))

# PCA
## Spread so Phenotypes are columns
CS1e <- CS1d %>% dplyr:: select(-contains("Phenotype")) %>% dplyr:: select(-contains("Tissue.Category")) %>% dplyr:: select(-contains("Rank"))
CS2 <- spread(CS1e, "Parameter", "Average.cell.Density")

### Remove identifiers
CS3 <- data.frame(CS2[, names(CS2) != "Slide.ID"], row.names = CS2[, names(CS2) == "Slide.ID"])
pca1 <- CS3 %>% dplyr:: select(-contains("All")) %>% dplyr:: select(-contains("Other"))
pca1[is.na(pca1)] <- 0

# Measure Collinear variables
corrplot(cor(pca1), type = 'lower')

### Complete PCA for all parameters
prin_comp <- prcomp(pca1[, names(pca1) != "Subtype"], scale. = T)
Subtype <- pca1[, "Subtype"]
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

g <- g + ggtitle("PCA of Immune Cell Abundance (Density per Megapixel)")
ggsave("PCA of Immune Cell Abundance (Density per Mpix).png" , plot = g, device = "png",
       path = "/Users/JackMcMurray/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC/Figures",
       height = 6, width = 6, units = 'in', dpi = 600)

## Plot Vectors of PCA
theta <- seq(0,2*pi, length.out = 100)
circle <- data.frame(x = cos(theta), y = sin(theta))
p <- ggplot(circle,aes(x,y)) + geom_path()

loadings <- data.frame(prin_comp$rotation, 
                       Parameters = row.names(prin_comp$rotation))
p + geom_text(data = loadings, 
              mapping = aes(x = PC1, y = PC2, label = Parameters, colour = Parameters)) +
  coord_fixed(ratio=1) +
  labs(x = "PC1", y = "PC2") + theme_bw()


# How many PC's are important?
screeplot(prin_comp)
fviz_screeplot(prin_comp, ncp = 10, choice = "eigenvalue")

# Find the Eigenvalues
eig <- (prin_comp$sdev)^2

## Variances in percentage
variance <- eig*100/sum(eig)

## Cumulative variances
cumvar <- cumsum(variance)

## Store the variances as a dataframe
### Write as a DF
eigenvalues_Cell_Density <- data.frame(eig = eig, variance = variance,
                      cumvariance = cumvar)
writeCsvO(eigenvalues_Cell_Density)

# Store the variances
var <- get_pca_var(prin_comp)

## Find the coordinates of variables
var$coord[, 1:5]

## Find the correlation between variables and principal components
loadings <- prin_comp$rotation
sdev <- prin_comp$sdev

load <- as.data.frame(loadings)

load[order(load$PC14, decreasing = T)[1:10],]
### Find the correlataions
var.coord <- t(apply(loadings, 1, var_cor_func, sdev))

## Calculate the Cos2 (square of the coordinates)
var.cos2 <- var.coord^2

## Contribution of variables to the Components ((cos2) * 100) / (total cos2 of the component)
comp.cos2 <- apply(var.cos2, 2, sum)
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
contrib.var <- as.data.frame(var.contrib)

## Find the most contributing variable
contrib.var[order(contrib.var$PC14, decreasing = T)[1:10],]

# Use the analysis to compare the values for the groups across PC1-3
## Store the contribution to each PC
PC_CD <- as.data.frame(prin_comp$x)
PC_CD1 <- tibble:: rownames_to_column(PC_CD, "Slide.ID")
PC_CD1$Slide.ID <- as.factor(PC_CD1$Slide.ID)

## Gather the PCs 
PC_CD2 <- PC_CD1 %>% gather(contains("PC"), key = Component, value = "ComponentScore")
PC_CD2$Component <- as.factor(PC_CD2$Component)

# Plot
## Test for normality
PC_CD2[grep("dMMR$", PC_CD2$Slide.ID), "Subtype"] <- "MSI"
PC_CD2[grep("pMMR$", PC_CD2$Slide.ID), "Subtype"] <- "MSS"
PC_CD2[grep(".*N[ ]S[[:digit:]]{6}$", PC_CD2$Slide.ID), "Subtype"] <- "hiCIRC"
PC_CD2$Subtype <- as.factor(PC_CD2$Subtype)

Normality <- PC_CD2
Normality$uniq <- as.factor(paste(Normality$Subtype, Normality$Component, sep = ","))

## Shapiro test
normal_list <- list()
c <- 1
for(i in levels(Normality$uniq)){
  name <- basename(i)
  cat('Processing', i, '\n')
  working <- droplevels(subset(Normality, uniq == i))
  stat <- shapiro.test(working[,"ComponentScore"])
  normal <- as.data.frame(stat$p.value)
  normal_list[[i]] <- cbind(i, normal)
  c <- c + 1
}
z <- do.call(rbind, normal_list)
rownames(z) <- c()

## QQ plots
for(i in levels(Normality$uniq)){
  working <- droplevels(subset(Normality, uniq == i))
  temp_plot <- ggplot(working) +
    stat_qq(aes(sample = ComponentScore))
  filen <- paste0(i,".png")
  ggsave(filen, plot = temp_plot, device = "png",
         path = "/Users/JackMcMurray/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC/Figures/Cell_Density/PCA/Normality",
         height=5, width=5, units='in', dpi=600)
  }

# Store the components
PC_CD2$Rank <- rank(PC_CD2$ComponentScore)

for(i in levels(PC_CD2$Component)){
  name <- basename(i)
  cat('Processing', i, '\n')
  Chosen <- droplevels(subset(PC_CD2, Component == i))
  YTitle <- paste0("Rank of ", i)
  MainTitle <- paste0("Violin plot of ", i)
  temp_plot <- ggSubtype(Chosen, YTitle, MainTitle)
  filen <- paste0(i,".png")
  ggsave(filen, plot = temp_plot, device = "png",
         path = "/Users/JackMcMurray/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC/Figures/Cell_Density/Per_Megapixel/PCA/PC_Comp",
         height = 5, width = 5, units = 'in', dpi = 600)
}

PC_CD2[grep("dMMR$", PC_CD2$Slide.ID), "Subtype"] <- "MSI"
PC_CD2[grep("pMMR$", PC_CD2$Slide.ID), "Subtype"] <- "MSS"
PC_CD2[grep(".*N[ ]S[[:digit:]]{6}$", PC_CD2$Slide.ID), "Subtype"] <- "hiCIRC"
PC_CD2$Subtype <- as.factor(PC_CD2$Subtype)

# Combine stats into a list
stat_list <- list()
c <- 1
for(i in levels(PC_CD2$Component)) {
  name <- basename(i)
  cat('Processing', i, '\n')
  workingon <- droplevels(subset(PC_CD2, Component == i))
  # assign your ggplot call to the i'th position in the list
  x <- compare_means(ComponentScore ~ Subtype, data = workingon, method = "wilcox.test")
  y <- as.data.frame(x)
  stat_list[[i]]  <- cbind(i, y)
  c <- c + 1}

# Bind and remove row names
z <- do.call(rbind, stat_list)
rownames(z) <- c()

# Separate
CellDensity_Megapixel_PCAStats <- z

# Write out the statistics
writeCsvO(CellDensity_Megapixel_PCAStats)

# Completing individual, tissue category PCAs
pcaEpi <- pca1 %>% dplyr:: select(-contains("Stroma"))
pcaEpi[is.na(pcaEpi)] <- 0

### Complete PCA for Epithelium
prin_comp <- prcomp(pcaEpi[, names(pcaEpi) != "Subtype"], scale. = T)
Subtype <- pcaEpi[, "Subtype"]
g <- ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
              groups = Subtype, ellipse = TRUE,
              circle = T,
              var.axes = F,
              choices = c(1,2)
)
g <- g + scale_color_manual(values = cbcols)
g <- g + theme_bw()
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')

g <- g + ggtitle("PCA of Immune Cell Abundance in the Stroma (Density per Megapixel)")
g
head(pcaEpi)





# LDA - are they normal?
head(pca1)
fac <- pca1
att1 <- fac %>% dplyr:: select(-contains("Slide.ID")) %>% dplyr:: select(-contains("Other")) #%>% dplyr:: select(-contains("All"))

## Run a PCA
att1[is.na(att1)] <- 0
pca <- prcomp(att1[,-1],
              center = TRUE,
              scale. = TRUE) 

prop.pca <- pca$sdev^2/sum(pca$sdev^2)

## Use PCA to help make LDA
r <- lda(Subtype ~ ., 
         att1, 
         prior = c(1,1,1)/3)

prop.lda <- r$svd^2/sum(r$svd^2)

plda <- predict(object = r,
                newdata = att1)

dataset <- data.frame(Type = att1[,"Subtype"],
                      pca = pca$x, lda = plda$x)


p1 <- ggplot(dataset) + geom_point(aes(lda.LD1, lda.LD2, colour = Type), size = 2.5) + 
  labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
       y = paste("LD2 (", percent(prop.lda[2]), ")", sep="")) + theme_bw()


# tsne
library(Rtsne)
# tSNE
head(pca1)
tsne <- Rtsne(as.matrix(pca1[, names(pca1) != "Subtype"]), dims = 10, perplexity = 8, verbose=TRUE, max_iter = 500)
plot(tsne$Y, main = "tsne", col = pca1$Subtype, pch = 16,
     axes = F)
