# Source Functions.R for packages and Functions
source("Functions.R")

# Search for ACW5
base_path <- "/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Export_Fixed/Cell_seg_only"

# create a dataframe with your phenotypes in
Phenotypes <- c("CD4", "CD8", "CD20", "CD68", "FOXP3", "DAPI")
pairs <- find_pheno_comb(Phenotypes)
categories <- c("Tumour", "Stroma")
radii <- c(10, 30)

# Perform these to get AveCountsWithin dataframe
AveCountsWithin <- count_within_batch(base_path, pairs, radii, category = categories, verbose = T)

read_cell_seg_data("/Volumes/ResearchData/Vectra/Vectra3/Fedor/nahla/Export_Fixed/Cell_seg_only/S0228449 25B 5PLEX_[45507,13296].txt")
## writeCsvO(AveCountsWithin)
ACW <- read.csv("Output/AveCountsWithin.csv")
ACW <- droplevels(subset(ACW, from != "DAPI"))
ACW <- droplevels(subset(ACW, to != "DAPI"))

# Factor
Factor <- c("slide_id",
            "source",
            "category",
            "from",
            "to")
ACW1 <- factorthese(ACW, Factor)

# Create parameter column
ACW1$parameter <- as.factor(paste(ACW1$from, ACW1$radius, ACW1$to, ACW1$category, ACW1$slide_id, sep = "/"))

# Loop for average number per slide
ACW2 <- ACW1
output <- data.frame(Parameter = character(),
                     MeanNum = double(),
                     stringsAsFactors = F)
c <- 1
for(i in levels(ACW2$parameter)){
  name <- basename(i)
  cat('Processing', i, '\n')
  working <- droplevels(subset(ACW2, parameter == i))
  aver <- mean(working$within_mean, na.rm = T)
  output[c, "Parameter"] <- i
  output[c, "MeanNum"] <- aver
  c <- c + 1
}

ACW3 <- output %>%
  separate(Parameter, c("PhenoFrom", "Radius", "PhenoTo", "Tissue.Category", "Slide.ID"), sep = "/")

# Get rid of [+] and spaces
ACW3$PhenoFrom <- gsub(" [+] ", "...", ACW3$PhenoFrom)
ACW3$PhenoFrom <- gsub(" ", ".", ACW3$PhenoFrom)
ACW3$PhenoTo <- gsub(" [+] ", "...", ACW3$PhenoTo)
ACW3$PhenoTo <- gsub(" ", ".", ACW3$PhenoTo)

somecolumns <- c("PhenoFrom",
                 "PhenoTo",
                 "Radius",
                 "Tissue.Category",
                 "Subtype",
                 "Slide.ID")

ACW4 <- factorthese(ACW3, somecolumns)

# Paste a parameter back together
ACW4$Parameter <- as.factor(paste(ACW4$PhenoFrom, ACW4$Radius, ACW4$PhenoTo, ACW4$Tissue.Category, sep = "/"))

# Spread
# Drop what you don't need
drop <- c("Tissue.Category",
          "Radius",
          "PhenoFrom",
          "PhenoTo")
ACW5 <- droplevels(ACW4[,!(names(ACW4) %in% drop)])
writeCsvO(ACW5)

ACW5 <- read.csv("Output/ACW5.csv")

# Make a ggplot for all
ACW5a <- ACW5
ACW5a$Rank <- rank(ACW5a$MeanNum)
ACW5a$Parameter <- as.factor(gsub("/", ",", ACW5a$Parameter))
ACW5a$Parameter <- as.factor(gsub("[...]", "-", ACW5a$Parameter))

# Plot
for(i in levels(ACW5a$Parameter)){
  name <- basename(i)
  cat('Processing', i, '\n')
  Chosen <- droplevels(subset(ACW5a, Parameter == i))
  YTitle <- paste0("Rank of ", i)
  MainTitle <- paste0("Violin plot of ", i)
  temp_plot <- ggSubtype(Chosen, YTitle, MainTitle)
  filen <- paste0(i,".png")
  ggsave(filen, plot = temp_plot, device = "png",
         path = "/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC/Figures/AverageCountsWithin/Raw",
         height = 6, width = 6, units = 'in', dpi = 600)
}

# Spread
ACW6 <- spread(ACW5, key = "Parameter", value = "MeanNum")

## Make Slide.ID the Row.Names
ACW7 <- data.frame(ACW6[, names(ACW6) != "Slide.ID"], row.names = ACW6[, names(ACW6) == "Slide.ID"])
ACW8 <- ACW7 %>% dplyr:: select(-contains("Other"))

# Check Collinearity
corrplot(cor(ACW8[, names(ACW8) != "Subtype"]), type = "lower")

# Linear Regression
model <- bayesglm(Subtype ~., family = binomial (link = 'logit'), data = ACW8)
display(model)

# PCA
prin_comp <- prcomp(ACW8[, names(ACW8) != "Subtype"], scale. = T)
plot(prin_comp)

Subtype <- ACW8[, "Subtype"]
g <- ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
              groups = Subtype, ellipse = T,
              circle = F,
              var.axes = F
)
g <- g + scale_color_manual(values = cbcols)
g <- g + theme_bw()
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
g <- g + ggtitle("PCA of the Average Number of Immune Cells 
                 within a 30 micron Radius")

ggsave("PCA of Radius.png" , plot = g, device = "png",
       path = "/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC/Figures/AverageCountsWithin/PCA",
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
Eigenvalues.ACW8 <- data.frame(eig = eig, variance = variance,
                      cumvariance = cumvar)
writeCsvO(Eigenvalues.ACW8)

# Store the variances
var <- get_pca_var(prin_comp)

## Find the coordinates of variables
var$coord[, 1:4]

## Find the correlation between variables and principal components
loadings <- prin_comp$rotation
sdev <- prin_comp$sdev


load <- as.data.frame(loadings)
load[order(load$PC1, decreasing = T)[1:5],]


### Find the correlataions
var.coord <- t(apply(loadings, 1, var_cor_func, sdev))

## Calculate the Cos2 (square of the coordinates)
var.cos2 <- var.coord^2

## Contribution of variables to the Components ((cos2) * 100) / (total cos2 of the component)
comp.cos2 <- apply(var.cos2, 2, sum)
var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
contrib.var <- as.data.frame(var.contrib)

## Find the most contributing variable
contrib.var[order(contrib.var$PC1, decreasing = T)[1:5],]

# Store the contribution to each PC
PC_ACW <- as.data.frame(prin_comp$x)
PC_ACW1 <- tibble:: rownames_to_column(PC_ACW, "Slide.ID")
PC_ACW1$Slide.ID <- as.factor(PC_ACW1$Slide.ID)

## Gather the PCs 
PC_ACW2 <- PC_ACW1 %>% gather(contains("PC"), key = Component, value = "ComponentScore")
PC_ACW2$Component <- as.factor(PC_ACW2$Component)
PC_ACW2$Rank <- rank(PC_ACW2$ComponentScore)

# Plot
for(i in levels(PC_ACW2$Component)){
  name <- basename(i)
  cat('Processing', i, '\n')
  Chosen <- droplevels(subset(PC_ACW2, Component == i))
  YTitle <- paste0("Rank of ", i)
  MainTitle <- paste0("Violin plot of ", i)
  temp_plot <- ggSubtype(Chosen, YTitle, MainTitle)
  filen <- paste0(i,".png")
  ggsave(filen, plot = temp_plot, device = "png",
         path = "/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC/Figures/AverageCountsWithin/PCA/PC_Comp",
         height = 6, width = 6, units = 'in', dpi = 600)
}

PC_ACW2[grep("dMMR$", PC_ACW2$Slide.ID), "Subtype"] <- "MSI"
PC_ACW2[grep("pMMR$", PC_ACW2$Slide.ID), "Subtype"] <- "MSS"
PC_ACW2[grep(".*N[ ]S[[:digit:]]{6}$", PC_ACW2$Slide.ID), "Subtype"] <- "hiCIRC"
PC_ACW2$Subtype <- as.factor(PC_ACW2$Subtype)

# Combine stats into a list
stat_list <- list()
c <- 1
for(i in levels(PC_ACW2$Component)) {
  name <- basename(i)
  cat('Processing', i, '\n')
  workingon <- droplevels(subset(PC_ACW2, Component == i))
  # assign your ggplot call to the i'th position in the list
  x <- compare_means(Rank ~ Subtype, data = workingon, method = "wilcox.test")
  y <- as.data.frame(x)
  stat_list[[i]]  <- cbind(i, y)
  c <- c + 1}

# Bind and remove row names
z <- do.call(rbind, stat_list)
rownames(z) <- c()

# Separate
AverageCounts_PCAStats <- z

# Write out the statistics
writeCsvO(AverageCounts_PCAStats)



# Plot
## Test for normality
PC_ACW2[grep("dMMR$", PC_ACW2$Slide.ID), "Subtype"] <- "MSI"
PC_ACW2[grep("pMMR$", PC_ACW2$Slide.ID), "Subtype"] <- "MSS"
PC_ACW2[grep(".*N[ ]S[[:digit:]]{6}$", PC_ACW2$Slide.ID), "Subtype"] <- "hiCIRC"
PC_ACW2$Subtype <- as.factor(PC_ACW2$Subtype)

Normality <- PC_ACW2
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
         path = "/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC/Figures/AverageCountsWithin/PCA/Normality",
         height = 6, width = 6, units = 'in', dpi = 600)
}