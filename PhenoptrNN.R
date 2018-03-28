# Source Functions.R for packages and Functions
source("Functions.R")

# Search for NN9 if already previously ran

# NEAREST NEIGHBOUR
## Convert old DF to a tsv file
csvtotsv("/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Extra/Nahla_Analysis/Output/df2a.csv",
         "/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Extra/Nahla_Analysis/Output/NN/DF1_tab.tsv")

## Calculate distances
NearNeigh <- compute_all_nearest_distance("/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Extra/Nahla_Analysis/Output/NN/DF1_tab.tsv")

# Write a TSV or a CSV?
# write.table(NearNeigh, file = "/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Extra/Nahla_Analysis/Output/NN/NearNeigh.tsv", row.names = F, sep = "\t")
# writeCsvO(NearNeigh)
# Actually have to read this in, else it's as a tibble
NearNeigh <- read.csv("Output/NearNeigh.csv")

# Remove stuff that doesn't matter anymore
NN <- NearNeigh %>% 
  dplyr:: select(-contains("Opal")) %>% 
  dplyr:: select(-contains("Total.Weighting")) %>% 
  dplyr:: select(-contains("pixels"))
drop <- c("tag",
          "Cell.ID",
          "Cell.X.Position",
          "Cell.Y.Position",
          "Category.Region.ID",
          "Confidence",
          "Batch",
          "Distance.from.Tissue.Category.Edge..microns")
NN1 <- droplevels(NN[,!(names(NN) %in% drop)])
NearestNeighbour_all_Cells <- NN1
# writeCsvO(NearestNeighbour_all_Cells)

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
NN4 <- NN3 %>%
  separate(Comb, c("Tissue.Category", "PhenoFrom", "PhenoTo", "Sample.Name"), "/")

## Make the PhenoTo and PhenoFrom the same
NN4$PhenoTo <- gsub("Distance.to.", "", NN4$PhenoTo)
NN4$PhenoTo <- gsub(" ", ".", NN4$PhenoTo)
NN4$PhenoFrom <- gsub(" ", ".", NN4$PhenoFrom)

## Remove rows where either PhenoTo or PhenoFrom is DAPI
NN4$PhenoTo <- as.factor(NN4$PhenoTo)
NN4$PhenoFrom <- as.factor(NN4$PhenoFrom)

NN4a <- droplevels(subset(NN4, PhenoTo != "DAPI"))
NN5 <- droplevels(subset(NN4a, PhenoFrom != "DAPI"))

NearestNeighbour_no_DAPI_relationships <- NN5
# writeCsvO(NearestNeighbour_no_DAPI_relationships)

# Paste back together
NN5$Comb <- as.factor(paste(NN5$Tissue.Category, NN5$PhenoFrom, ".Distance.To.", NN5$PhenoTo, NN5$Sample.Name, sep = "/"))

# Factor them all
NN5$Comb <- as.factor(NN5$Comb)

## Loop and find the mean distance (for the replicates to get a biological replicate)
output <- data.frame(Comb = character(),
                     MeanDist = double(),
                     Type = character(),
                     Slide.ID = character(),
                     stringsAsFactors = FALSE)

c <- 1
for(i in levels(NN5$Comb)){
  name <- basename(i)
  cat('Processing', i, '\n')
  df <- droplevels(subset(NN5, Comb == i))
  y <- mean(df$Distance, na.rm = T)
  output[c, "Comb"] <- i
  output[c, "MeanDist"] <- y
  output[c, "Type"] <- as.character(levels(df$Type))
  output[c, "Slide.ID"] <- as.character(levels(df$Slide.ID))
  c <- c + 1}

# Separate back into original columns
NN6 <- output %>%
  separate(Comb, c("Tissue.Category", "PhenoFrom","Distance.To", "PhenoTo", "Sample.Name"), "/")

# Combine and reloop for mean nearest distance on a slide basis
NN6$Parameter <- as.factor(paste(NN6$Tissue.Category, NN6$PhenoFrom, NN6$Distance.To, NN6$PhenoTo, NN6$Slide.ID, sep = "/"))
NN6$Type <- as.factor(NN6$Type)
NN6a <- data.frame(Parameter = character(),
                     Distance = double(),
                     Type = character(),
                     stringsAsFactors = FALSE)


c <- 1
for(i in levels(NN6$Parameter)){
  name <- basename(i)
  cat('Processing', i, '\n')
  df <- droplevels(subset(NN6, Parameter == i))
  y <- mean(df$MeanDist, na.rm = T)
  NN6a[c, "Parameter"] <- i
  NN6a[c, "Distance"] <- y
  NN6a[c, "Type"] <- as.character(levels(df$Type))
  c <- c + 1}

# Remove Slide.ID from combined Parameter column
NN7 <- NN6a %>%
  separate(Parameter, c("Tissue.Category", "PhenoFrom","Distance.To", "PhenoTo", "Slide.ID"), "/")
NN7$Parameter <- as.factor(paste(NN7$Tissue.Category, NN7$PhenoFrom, NN7$Distance.To, NN7$PhenoTo, sep = "_"))

# Drop what you don't need
drop <- c("Tissue.Category",
          "Distance.To",
          "PhenoFrom",
          "PhenoTo")
NN8 <- droplevels(NN7[,!(names(NN7) %in% drop)])
NN8$Slide.ID <- as.factor(NN8$Slide.ID)
NN8$Type <- as.factor(NN8$Type)

# Spread
NN9 <- spread(NN8, key = "Parameter", value = "Distance")

# Change the Types to proper types
NN9[grep("dMMR$", NN9$Slide.ID), "Subtype"] <- "MSI"
NN9[grep("pMMR$", NN9$Slide.ID), "Subtype"] <- "MSS"
NN9[grep(".*N[ ]S[[:digit:]]{6}$", NN9$Slide.ID), "Subtype"] <- "hiCIRC"
NN9$Subtype <- as.factor(NN9$Subtype)

# writeCsvO(NN9)
NN9 <- read.csv("Output/NN9.csv")

## Make the Slide.ID the row.names
NN9a <- data.frame(NN9[, names(NN9) != "Slide.ID"], row.names = NN9[, names(NN9) == "Slide.ID"])
NN10 <- NN9a %>% dplyr:: select(-contains("Other")) %>% dplyr:: select(-contains("Slide.ID")) %>% dplyr:: select(-contains("Type", ignore.case = F))

# randomForest imputing
library(randomForest)
NN_impute <- rfImpute(Subtype ~ ., NN10)

# writeCsvO(NN_impute)

# Logistic Regression
corrplot(cor(NN_impute[, names(NN_impute) != "Subtype"]), type = 'lower')

NN_impute$Subtype <- as.factor(NN_impute$Subtype)
library(arm)
model <- bayesglm(Subtype ~., family = binomial (link = 'logit'), data = NN_impute)
display(model)


# PCA
### Complete PCA
prin_comp <- prcomp(NN_impute[, names(NN_impute) != "Subtype"], scale. = T)
plot(prin_comp)

Subtype <- NN_impute[, "Subtype"]
g <- ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
              groups = Subtype, ellipse = T,
              circle = F,
              var.axes = F
)
g <- g + scale_color_manual(values = cbcols)
g <- g + theme_bw()
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
g <- g + ggtitle("PCA of the Mean Nearest Distance Between Immune Subsets")
g
ggsave("PCA of Nearest Neighbour.png" , plot = g, device = "png",
       path = "/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC/Figures/NearestNeighbour",
       height = 6, width = 6, units = 'in', dpi = 600)


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
  cat('Processing', i, '\n')
  Chosen <- droplevels(subset(PC_NN3, Component == i))
  YTitle <- paste0("Rank of ", i)
  MainTitle <- paste0("Violin plot of ", i)
  temp_plot <- ggSubtype(Chosen, YTitle, MainTitle)
  filen <- paste0(i,".png")
  ggsave(filen, plot = temp_plot, device = "png",
         path = "/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC/Figures/NearestNeighbour",
         height=5, width=5, units='in', dpi=600)
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
