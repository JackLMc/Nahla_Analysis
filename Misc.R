# Source Functions.R for packages and Functions
source("Functions.R")


# Gain the dataframe
## Source
### Source("ToCSV.R")

## Import
DF1 <- read.csv("Output/df2a.csv")

# Assuming that the cells are circular
diameter <- find_diameter(DF1)
DF <- cbind(DF1, diameter)

# Correct the Column (numeric)
DF$Distance.from.Tissue.Category.Edge..pixels. <- as.character(DF$Distance.from.Tissue.Category.Edge..pixels.)
DF$Distance.from.Tissue.Category.Edge..pixels.[DF$Distance.from.Tissue.Category.Edge..pixels. == "#N/A"]  <- NA 
DF[DF=="<NA>"]=NA # Takes a while, but seems to do the job
DF$Distance.from.Tissue.Category.Edge..pixels. <- as.numeric(DF$Distance.from.Tissue.Category.Edge..pixels.)

# Combine into a single, loopable fashion
DF$Combined <- as.factor(paste(DF$Slide.ID, DF$Tissue.Category, DF$Phenotype, sep = "/"))
DFa <- DF

# Loop around to find Min, Max, Median and Mean
output <- data.frame(Combine = character(),
                     Min = double(),
                     Max = double(),
                     Median = double(),
                     Mean = double(),
                     Subtype = character(),
                     stringsAsFactors = FALSE
)

c <- 1
for(i in levels(DFa$Combined)){
  name <- basename(i)
  cat('Processing', i, '\n')
  df <- droplevels(subset(DFa, Combined == i))
  output[c, "Combine"] <- i
  output[c, "Mean"] <- mean(df$Distance.from.Tissue.Category.Edge..pixels., na.rm = T)
  output[c, "Median"] <- median(df$Distance.from.Tissue.Category.Edge..pixels., na.rm = T)
  output[c, "Max"] <- max(df$Distance.from.Tissue.Category.Edge..pixels., na.rm = T)
  output[c, "Min"] <- min(df$Distance.from.Tissue.Category.Edge..pixels., na.rm = T)
  output[c, "Subtype"] <- as.character(levels(df$Type))
  c <- c + 1
}

# Separate back out
Dist_Next <- output %>%
  separate(Combine, c("Sample.Name", "Tissue.Category", "Phenotype"), "/")

# Get rid of NaN, Inf and -Inf
Dist_Next$Min <- as.character(Dist_Next$Min)
Dist_Next$Max <- as.character(Dist_Next$Max)
Dist_Next$Median <- as.character(Dist_Next$Median)
Dist_Next$Mean <- as.character(Dist_Next$Mean)

Dist_Next$Min[Dist_Next$Min == "NaN" |Dist_Next$Min == "Inf"| Dist_Next$Min == "-Inf"]  <- NA
Dist_Next$Max[Dist_Next$Max == "NaN" |Dist_Next$Max == "Inf"| Dist_Next$Max == "-Inf"]  <- NA
Dist_Next$Median[Dist_Next$Median == "NaN" |Dist_Next$Median == "Inf"| Dist_Next$Median == "-Inf"]  <- NA
Dist_Next$Mean[Dist_Next$Mean == "NaN" |Dist_Next$Mean == "Inf"| Dist_Next$Mean == "-Inf"]  <- NA

Dist_Next$Min <- as.numeric(Dist_Next$Min)
Dist_Next$Max <- as.numeric(Dist_Next$Max)
Dist_Next$Median <- as.numeric(Dist_Next$Median)
Dist_Next$Mean <- as.numeric(Dist_Next$Mean)

# Factor these
Factor <- c("Sample.Name",
            "Tissue.Category",
            "Phenotype",
            "Subtype")

Dist_Next_Cat <- factorthese(Dist_Next, Factor)

# Factor the combination + Remove individual columns
Dist_Next_Cat$PhenoSub <- as.factor(paste(Dist_Next_Cat$Phenotype, Dist_Next_Cat$Subtype, sep = "_"))
DNC <- Dist_Next_Cat %>% dplyr:: select(-contains("Subtype")) %>% dplyr:: select(-contains("Phenotype"))

## Drop to one tissue category for impution
Epi <- droplevels(subset(DNC, Tissue.Category == "Epithelium"))
Stro <- droplevels(subset(DNC, Tissue.Category == "Stroma"))

### Impute values based on Cell Phenotype, Subtype of slide and their tissue localisation
RFEpi <- rfImpute(PhenoSub~., Epi)
RFStro <- rfImpute(PhenoSub~., Stro)

# Bind back together, separate and rebind for graph making
RF <- rbind(RFEpi, RFStro)
IDNC <- RF %>%
  separate(PhenoSub, c("Phenotype", "Subtype"), "_")
Imputed_DistanceToNextCat <- IDNC %>% gather(key = "Measure", value = "value", Mean, Min, Median, Max)

Imputed_DistanceToNextCat$Category.To <- ifelse((Imputed_DistanceToNextCat$Tissue.Category == "Epithelium"), "Stroma", "Epithelium")
Imputed_DistanceToNextCat$PhenoCatMeasure <- as.factor(paste(Imputed_DistanceToNextCat$Measure, "distance of", Imputed_DistanceToNextCat$Phenotype, "to", Imputed_DistanceToNextCat$Category.To, sep = " "))
# writeCsvO(Imputed_DistanceToNextCat)
Imputed_DistanceToNextCat$Rank <- rank(Imputed_DistanceToNextCat$value)

# Loop to produce graphs of the Average distance to next category 
for(i in levels(Imputed_DistanceToNextCat$PhenoCatMeasure)){
  print(paste("Working on:", i))
  Working <- droplevels(subset(Imputed_DistanceToNextCat, PhenoCatMeasure == i))
  YTitle <- paste0("Rank of ", i)
  MainTitle <- paste0("Violin plot of ", i)
  temp_plot <- ggSubtype(Working, YTitle, MainTitle)
  filen <- paste0(i,".png")
  ggsave(filen, plot = temp_plot, device = "png",
         path = "/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC/Figures/Dist_To_Next",
         height = 5, width = 5, units = 'in', dpi = 600)
}



# Average Tumour Burden
keep <- c("Entire.Cell.Area..pixels.",
          "Sample.Name",
          "Tissue.Category",
          "Phenotype",
          "Cell.ID",
          "Type",
          "Slide.ID")
           
DFb <- droplevels(DF[,(names(DF) %in% keep)])
DFb$uniqID <- as.factor(paste(DFb$Slide.ID, DFb$Tissue.Category, sep = "/"))

Output <- data.frame(Combine = character(),
                     Total.Cell.Area = double(),
                     stringsAsFactors = F)
c <- 1
for(i in levels(DFb$uniqID)){
  print(paste("Working on:", i))
  step1 <- droplevels(subset(DFb, uniqID == i))
  step2 <- sum(step1$Entire.Cell.Area..pixels.)
  Output[c, "Combine"] <- i
  Output[c, "Total.Cell.Area"] <- step2
  c <- c + 1
}

DFc <- Output %>%
  separate(Combine, c("Slide.ID", "Tissue.Category"), "/")

# Load this in, gains total.category.burdens 
CS2c <- read.csv("Output/Cat_Burden.csv")
remove <- c("Phenotype",
            "Rank")
CS2d <- droplevels(CS2c[,!(names(CS2c) %in% remove)])
CS2e <- droplevels(subset(CS2d, Tissue.Category != "All"))

## Merge the dataframes
DS <- merge(DFc, CS2e, by = c("Slide.ID", "Tissue.Category"))

## Calculate the Cell.Free Area
DS$Total.Free.Area <- DS$Total.Cat.Burden - DS$Total.Cell.Area

## Gather Parameters
DS1 <- DS %>% gather(contains("Total."), key = "Parameter", value = "Area")
DS1$Parameter <- as.factor(paste(DS1$Parameter, DS1$Tissue.Category, sep = ","))
DS1$Rank <- rank(DS1$Area)

# Plot
for(i in levels(DS1$Parameter)){
  print(paste("Working on:", i))
  Chosen <- droplevels(subset(DS1, Parameter == i))
  YTitle <- paste0("Rank of ", i)
  MainTitle <- paste0("Violin plot of ", i)
  temp_plot <- ggSubtype(Chosen, YTitle, MainTitle)
  filen <- paste0(i,".png")
  ggsave(filen, plot = temp_plot, device = "png",
         path = "/Users/jlm650/OneDrive/University_of_Birmingham/PhD/Vectra_MSI_MSS_hiCIRC/Figures/Checking_burdens",
         height = 5, width = 5, units = 'in', dpi = 600)
}
