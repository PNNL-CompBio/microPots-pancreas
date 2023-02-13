library(readxl)
library(dplyr)

#1.3.1
# Create fractions table
datasets <- unique(intersect(msnid$Dataset, masic_data$Dataset))
fractions <- data.frame(Dataset = datasets) %>% 
  mutate(PlexID = "S1")

#1.3.2
# Create samples table
samples <- read.csv("Hubmap/Analysis/0/data/plexKey.csv") %>% 
  rename(ReporterName = Channel) %>% 
  mutate(PlexID = "S1", 
         QuantBlock = 1, 
         ReporterAlias = sprintf("S_%s_%s", X, Y), 
         MeasurementName = ReporterAlias) %>%
  select(-c(X,Y))

#1.3.3
# Create references table
# Create references table
references <- data.frame(PlexID = "S1", 
                         QuantBlock = 1,
                         Reference = 1)

#1.4.2
crosstab <- 2^crosstab
rowMedian <- apply(crosstab, 1, median, na.rm = T)
crosstab <- sweep(crosstab, 1, rowMedian, FUN = "/")
crosstab <- log2(crosstab)

#Norm
normcoef <- apply(crosstab, 2, median, na.rm = T)
crosstab <- sweep(crosstab, 2, normcoef, FUN = "-") 

#V2
load("Hubmap/Analysis/0/data/0_processed_msnid_and_masic.RData")
##Make crosstab
rowMedian <- apply(crosstab, 1, median, na.rm = T)
norm_table <- crosstab-rowMedian
colMedian <- apply(norm_table, 2, median, na.rm = T)
final_norm_table <- t(crosstab)-colMedian
crosstab <- t(final_norm_table)

write.table(crosstab, file = "Hubmap/Analysis/0/data/0_global_crosstab_V2.txt", sep = "\t", quote = FALSE, row.names = TRUE)

boxplot(crosstab)



