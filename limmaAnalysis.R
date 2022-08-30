library(dplyr)
library(purrr)
library(MSnbase)
library(ggplot2)

islets <- c(0:4, 7, 10)

#files <- sprintf("%g/data/%g_global_crosstab.txt", islets, islets)

source("synapseUtil.R")
syn<-synapseLogin()

crosstabs<-list(`0`='syn30389393',`10`='syn30385467',`1`='syn30385349',`2`='syn30385407',
                `3`='syn30385441',`4`='syn30385450',`7`='syn30385456')

#there are also v2 crosstabs?
v2_crosstabs<-list(`0`='syn31972408',`10`='syn31972599',`1`='syn31972430',`2`='syn31972462',
                   `3`='syn31972520',`4`='syn31972565',`7`='syn31972584')

##load in crosstabs into single file
crosstabList <- lapply(names(crosstabs), function(file_i) 
  {
    islet_i <- file_i#sub("/.*", "", file_i)
    
    crosstab_i <- read.delim(syn$get(crosstabs[[file_i]])$path, row.names = 1)
    
    colnames(crosstab_i) <- paste0(islet_i, "_", colnames(crosstab_i))
    crosstab_i <- crosstab_i %>% tibble::rownames_to_column("feature")
    return(crosstab_i)
  })


##create single matrix
crosstab <- crosstabList %>% purrr::reduce(full_join, by = "feature") %>% 
  tibble::column_to_rownames("feature") %>% 
  as.matrix()


##read in metadata from synaspe
library(readxl)
isletMeta <- readxl::read_xlsx(syn$get('syn30070383')$path,sheet='Metadata table')%>%#"~/Desktop/2022/Paul/Hubmap/IsletMeta_New_correct.csv") %>% 
  mutate(X_Y = paste(`Islet Number`, "S", `X-coord`, `Y-coord`, sep = "_")) %>% 
  tibble::column_to_rownames("X_Y") %>% 
  mutate(`IsletNumber` = as.character(`Islet Number`))%>%
  dplyr::rename(IsletStatus='Islet status')
  
isletMeta <- isletMeta[colnames(crosstab), ]


#m1 <- m1[

library(ggplot2)
library(ggfortify)


crosstab <- crosstab[rowMeans(!is.na(crosstab)) >= .5, ]
crosstab[which(is.na(crosstab),arr.ind=TRUE)]<-0.0
crosstab <- crosstab[rowMeans(!is.na(crosstab)) >= .5, ]
pc<-prcomp(t(crosstab))

autoplot(pc,data=isletMeta,colour='IsletNumber')
autoplot(pc,data=isletMeta,colour='IsletStatus')

##Need to find the plot_pca funciton, is it in msnbase.utils?
#plot_pca(m1, phenotype = "Islet_Number") #no batch effect
#plot_pca(m1, phenotype = "Islet_status") #check biplot = T?

#save(m1, file = "combinedIslets_Results_New/combinedMSnSet.RData", compress = T)

m1 <- MSnSet(exprs = crosstab, pData = isletMeta)
#Run limma
res <- limma_contrasts(eset = m1, model.str = "~ 0 + Islet_status", 
                       coef.str = "IsletStatus", contrasts = "IsletStatusProximal - IsletStatusDistal", trend = T, robust = T) #plot = T?

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



