# Source Functions.R for packages and Functions
source("Functions.R")

# Read CSV
CS1d <- read.csv("Output/Cell_Density_DF.csv")
ACW5 <- read.csv("Output/ACW5.csv")


CellDense <- CS1d %>% dplyr:: select(-contains("Rank")) %>% dplyr:: select(-contains("Tissue.Category")) %>% dplyr:: select(-contains("Phenotype"))
AverCount <- spread(ACW5, key = "Parameter", value = "MeanNum")
CD1 <- spread(CellDense, key = "Parameter", value = "Average.cell.Density")

# Merge and replace NA for 0 (they're actually 0)
Comb <- merge(CD1, AverCount, by = c("Slide.ID", "Subtype"))
Comb[is.na(Comb)] <- 0
DF1 <- Comb %>% dplyr:: select(-contains("Other")) %>% dplyr:: select(-contains("50")) %>% dplyr:: select(-contains("All"))

## Start PCA for Counts and Radius
DF2 <- data.frame(DF1[, names(DF1) != "Slide.ID"], row.names = DF1[, names(DF1) == "Slide.ID"])
prin_comp <- prcomp(DF2[, names(DF2) != "Subtype"], scale. = T)

Subtype <- DF2[, "Subtype"]
g <- ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
              groups = Subtype, ellipse = T,
              circle = F,
              var.axes = F,
              choices = c(1,2)
)
g <- g + scale_color_manual(values = cbcols)
g <- g + theme_bw()
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
g <- g + ggtitle("PCA of Immune Cell Density and Location")
g
g + geom_vline(xintercept = 0) 
g + geom_hline(yintercept = 0)

ggsave("PCA of Immune Cell Density and Interactions.png" , plot = g, device = "png",
       path = "/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC/Figures",
       height = 6, width = 6, units = 'in', dpi = 600)

# Read in Nearest Neighbour Stuff if choosing to impute
NN9 <- read.csv("Output/NN9.csv")
NN10 <- NN9 %>% dplyr:: select(-contains("Other")) %>% dplyr:: select(-contains("Type", ignore.case = F))
Combine <- merge(Comb, NN10, by = c("Slide.ID", "Subtype"))
DF1 <- Combine %>% dplyr:: select(-contains("Other")) %>% dplyr:: select(-contains("50")) %>% dplyr:: select(-contains("All"))

NN10
# randomForest imputing
head(DF1)

DF1a <- rfImpute(Subtype ~ ., DF1)

# Use RandomForest imputed with NN (Top) Or just ACW and CD (bottom)
DF2a <- data.frame(DF1a[, names(DF1a) != "Slide.ID"], row.names = DF1a[, names(DF1a) == "Slide.ID"])
prin_comp <- prcomp(DF2a[, names(DF2a) != "Subtype"], scale. = T)

Subtype <- DF2a[, "Subtype"]
g <- ggbiplot(prin_comp, obs.scale = 1, var.scale = 1, 
              groups = Subtype, ellipse = T,
              circle = F,
              var.axes = F
)
g <- g + scale_color_manual(values = cbcols)
g <- g + theme_bw()
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
g <- g + ggtitle("PCA of Immune Cell Density and Location")
# g + geom_vline(xintercept = 0) 
# g + geom_hline(yintercept = 0)

ggsave("PCA of Immune Cell Density, Location and Interactions.png" , plot = g, device = "png",
       path = "/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC/Figures",
       height = 6, width = 6, units = 'in', dpi = 600)

# PCA interrogation (Remember which set you're investigating with, both inclusion of NN and exclusion are prin_comp)
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
contrib.var[order(contrib.var$PC2, decreasing = T)[1:10],]


# Plot two biggest predictors against one another
DF2a <- tibble:: rownames_to_column(DF2, "Slide.ID")

# Try gathering and then ranking to plot? (Or is that not ok?)
DF3 <- DF2a %>% gather(-"Subtype", - "Slide.ID", key = "Parameter", value = "NormalisedVal")
DF3$Rank <- rank(DF3$NormalisedVal)
DF4 <- DF3 %>% dplyr:: select(-contains("Normal"))
DF5 <- spread(DF4, key = "Parameter", value = "Rank")


g <- ggplot(DF5, aes(y = CD4.10.CD8...T.bet.Epithelium, x = Epithelium.CD8...T.bet, color = Subtype))+
  geom_point(
    alpha = 0.8,
    size = 4)+
  labs(x = "Tc1 % in Epithelium",
       y = "Tc1 within 10 microns of CD4 in Epithelium"
  )+
  theme_bw()+ 
  geom_text(aes(x = 3000, y = 500, label = lm_eqn(lm(Epithelium.CD8...T.bet ~ CD4.10.CD8...T.bet.Epithelium, DF5))), parse = TRUE) +
  
  #stat_ellipse() +
  scale_color_manual(values = cbcols) 


ggsave("Predictors_of_PC1.png" , plot = g, device = "png",
       path = "/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC/Figures",
       height = 6, width = 6, units = 'in', dpi = 600)

# Store the contribution to each PC
PC_ALL <- as.data.frame(prin_comp$x)
PC_ALL1 <- tibble:: rownames_to_column(PC_ALL, "Slide.ID")
PC_ALL1$Slide.ID <- as.factor(PC_ALL1$Slide.ID)

## Gather the PCs 
PC_ALL2 <- PC_ALL1 %>% gather(contains("PC"), key = Component, value = "ComponentScore")
PC_ALL2$Component <- as.factor(PC_ALL2$Component)
PC_ALL2$Rank <- rank(PC_ALL2$ComponentScore)

# Plot
for(i in levels(PC_ALL2$Component)){
  name <- basename(i)
  cat('Processing', i, '\n')
  Chosen <- droplevels(subset(PC_ALL2, Component == i))
  YTitle <- paste0("Rank of ", i)
  MainTitle <- paste0("Violin plot of ", i)
  temp_plot <- ggSubtype(Chosen, YTitle, MainTitle)
  filen <- paste0(i,".png")
  ggsave(filen, plot = temp_plot, device = "png",
         path = "/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC/Figures/ALL",
         height = 5, width = 5, units = 'in', dpi = 600)
}

PC_ACW2[grep("dMMR$", PC_ACW2$Slide.ID), "Subtype"] <- "MSI"
PC_ACW2[grep("pMMR$", PC_ACW2$Slide.ID), "Subtype"] <- "MSS"
PC_ACW2[grep(".*N[ ]S[[:digit:]]{6}$", PC_ACW2$Slide.ID), "Subtype"] <- "hiCIRC"
PC_ACW2$Subtype <- as.factor(PC_ACW2$Subtype)
