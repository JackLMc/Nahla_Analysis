library(UsefulFunctions)


## Create Pies for Nahla stuff
clin <- read.csv("Data/New_Clin.csv")
levels(clin$grade)[levels(clin$grade)=="1"] <- "I"
levels(clin$grade)[levels(clin$grade)=="2"] <- "II"
levels(clin$grade)[levels(clin$grade)=="3"] <- "III"

# Grades
grade_percents <- data.frame(stringsAsFactors = F)
c <- 1
clin$grade <- as.factor(clin$grade)
for(i in levels(clin$grade)){
  work <- droplevels(subset(clin, grade == i))
  number <- nrow(work)
  grade_percents[c, "grade"] <- i
  grade_percents[c, "number_of_pats"] <- number
  c <- c + 1
}

library(tidyverse)
bp <- ggplot(grade_percents, aes(x="", y = number_of_pats, fill = grade))+
  geom_bar(width = 1, stat = "identity")

bp <- bp + coord_polar("y", start = 0)
bp <- bp + theme_blank() + theme(axis.text.x = element_blank()) 

cols <- c("I" = "#009E73", "II" = "#F0E442", "III" = "#D55E00")
bp <- bp + scale_fill_manual(values = cols)

ggsave("Grade_pie.pdf" , plot = bp, device = "pdf",
       path = "/Users/jlm650/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Paper/")

# HER2/ER2 status
clin$HER2_ER_Status <- ifelse((clin$ER == "positive" & clin$HER2 == "positive"), "doublepos", 
                              ifelse((clin$ER == "negative" & clin$HER2 == "negative"), "doubleneg",
                                     ifelse((clin$ER == "positive" & clin$HER2 == "negative"), "ER_pos", 
                                            ifelse((clin$ER == "negative" & clin$HER2 == "positive"), "HER2_pos", "this"))))

clin$HER2_ER_Status <- as.factor(clin$HER2_ER_Status)
HER2_ER_percents <- data.frame(stringsAsFactors = F)
c <- 1

for(i in levels(clin$HER2_ER_Status)){
  work <- droplevels(subset(clin, HER2_ER_Status == i))
  number <- nrow(work)
  HER2_ER_percents[c, "HER2_ER_Status"] <- i
  HER2_ER_percents[c, "number_of_pats"] <- number
  c <- c + 1
}

library(tidyverse)
bp <- ggplot(HER2_ER_percents, aes(x="", y = number_of_pats, fill = HER2_ER_Status))+
  geom_bar(width = 1, stat = "identity")

bp <- bp + coord_polar("y", start = 0)
bp <- bp + theme_blank() + theme(axis.text.x = element_blank()) 

cols <- c("doubleneg" = "#009E73", "doublepos" = "#F0E442", "ER_pos" = "#D55E00", "HER2_pos" = "#0072B2")
levels(clin$HER2_ER_Status)
bp <- bp + scale_fill_manual(values = cols)

ggsave("HER2_ER_Status_pie.pdf", plot = bp, device = "pdf",
       path = "/Users/jlm650/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Paper/")

# Overall Survival
dead_percents <- data.frame(stringsAsFactors = F)
c <- 1
clin$Dead <- as.factor(clin$Dead)
for(i in levels(clin$Dead)){
  work <- droplevels(subset(clin, Dead == i))
  number <- nrow(work)
  dead_percents[c, "Dead"] <- i
  dead_percents[c, "number_of_pats"] <- number
  c <- c + 1
}

library(tidyverse)
bp <- ggplot(dead_percents, aes(x="", y = number_of_pats, fill = Dead))+
  geom_bar(width = 1, stat = "identity")

bp <- bp + coord_polar("y", start = 0)
bp <- bp + theme_blank() + theme(axis.text.x = element_blank()) 

cols <- c("alive" = "#009E73", "dead" = "#F0E442")
bp <- bp + scale_fill_manual(values = cols)

ggsave("OS_pie.pdf", plot = bp, device = "pdf",
       path = "/Users/jlm650/OneDrive/UoB/PhD/Projects/5_Extra/Nahla_Analysis/Figures/Paper/")
