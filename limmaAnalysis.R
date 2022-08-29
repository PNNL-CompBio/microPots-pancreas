library(dplyr)
library(purrr)
library(MSnSet.utils)
library(ggplot2)

islets <- c(0:4, 7, 10)
files <- sprintf("%g/data/%g_global_crosstab.txt", islets, islets)

crosstabList <- lapply(files, function(file_i) 
  {
    islet_i <- sub("/.*", "", file_i)
    
    crosstab_i <- read.delim(file_i, row.names = 1)
    
    colnames(crosstab_i) <- paste0(islet_i, "_", colnames(crosstab_i))
    crosstab_i <- crosstab_i %>% tibble::rownames_to_column("feature")
    return(crosstab_i)
  })

crosstab <- crosstabList %>% purrr::reduce(full_join, by = "feature") %>% 
  tibble::column_to_rownames("feature") %>% 
  as.matrix()

isletMeta <- read.csv("~/Desktop/2022/Paul/Hubmap/IsletMeta_New_correct.csv") %>% 
  mutate(X_Y = paste(Islet_Number, "S", Xcoord, Ycoord, sep = "_")) %>% 
  tibble::column_to_rownames("X_Y") %>% 
  mutate(Islet_Number = as.character(Islet_Number))
  
isletMeta <- isletMeta[colnames(crosstab), ]

m1 <- MSnSet(exprs = crosstab, pData = isletMeta)
m1 <- m1[rowMeans(!is.na(exprs(m1))) >= .5, ]

plot_pca(m1, phenotype = "Islet_Number") #no batch effect
plot_pca(m1, phenotype = "Islet_status") #check biplot = T?

save(m1, file = "combinedIslets_Results_New/combinedMSnSet.RData", compress = T)

#Run limma
res <- limma_contrasts(eset = m1, model.str = "~ 0 + Islet_status", 
                       coef.str = "Islet_status", contrasts = "Islet_statusProximal - Islet_statusDistal", trend = T, robust = T) #plot = T?

write.table(res, "combinedIslets_Results_New/Proximal_Distal/limma/limmaResults.txt", quote = F, row.names = F, sep = "\t")

#now take a look at the p-values
hist(res$P.Value, breaks = seq(0,1,.025))
prop.table(table(res$P.Value < .025))
hist(res$adj.P.Val, breaks = seq(0,1,.025))
prop.table(table(res$adj.P.Val < .05))

#scatterplot
numProt <-colSums(!is.na(exprs(m1)))
names(numProt) <- sampleNames(m1)

ggplot(mapping = aes(x = names(numProt), y = numProt, color = m1$Islet_Number)) + 
  geom_point() + 
  geom_text(label = ifelse(numProt < 4000, names(numProt), NA), vjust = 1) +
  theme(axis.text.x = element_text(angle = 90))



