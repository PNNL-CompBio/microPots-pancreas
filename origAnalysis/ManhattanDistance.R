library(dplyr)
library(plyr)
library(tidyr)
library(tibble)
library(data.table)

load("combinedIslets_Results_New/combinedMSnSet.RData")
View(exprs(m1))
View(pData(m1))


m1$NewCol <- case_when(
  m1$Plex == "127C" ~ 1,
  m1$Plex == "129C" ~ 2,
  m1$Plex == "126" ~ 3,
  m1$Plex == "128N" ~ 0,
  m1$Plex == "129N" ~ 1,
  m1$Plex == "131C" ~ 2,
  m1$Plex == "127N" ~ 1,
  m1$Plex == "128C" ~ 2,
  m1$Plex == "130C" ~ 3)




manhattan_dist <- function(a, b){
  dist <- abs(a-b)
  dist <- sum(dist)
  return(dist)
}

a <- c(1,5,7)
b <- c(3,2,1)
c <- c(3,2,1)

mat <- data.frame("0" = c(1,NA,NA), "1" = c(3,2,4), "2" = c(3,2,1), "3" = c(c(4,2,NA)))
mat <- t(mat)
mat2 <- dist(mat, method = "manhattan")


View(exprs(m1_sub))
View(pData(m1))

m1_sub <- m1["sp|A0AV96|RBM47_HUMAN",m1$Islet_Number == "0"]

#######

x <- pData(m1) %>% select(Islet_Number, NewCol, Plex) %>% 
  cbind(., t(exprs(m1))) %>% 
  arrange(Islet_Number, NewCol, Plex) %>% 
  group_by(Islet_Number, NewCol) %>% 
  dplyr::mutate(idx = paste0("c", 1:n())) %>% 
  pivot_longer(cols = -c(Islet_Number, NewCol, Plex, idx), names_to = "Protein", values_to = "exprs") %>% 
  select(-Plex) %>% 
  pivot_wider(names_from = idx, values_from = exprs) %>% arrange(Islet_Number, Protein, NewCol) %>% 
  filter(!is.na(c1)) %>% 
  ungroup()

x2 <- split.data.frame(select(x, -Islet_Number), x$Islet_Number) %>% 
  lapply(function(xi)
    {
      split.data.frame(select(xi, -Protein), xi$Protein) %>% 
      lapply(function(pi)
        {
          tmp <- column_to_rownames(pi, "NewCol") %>% 
          as.matrix() %>% 
          dist(method = "manhattan") %>% as.matrix()
          
          tmp[upper.tri(tmp, diag = T)] <- NA
          
          tmp <- tmp %>% as.data.frame() %>% rownames_to_column("NewCol1") #%>% pivot_longer(cols = -NewCol1, names_to = "NewCol2", values_to = "mhDist")
          return(tmp)
        }) %>% rbindlist(idcol = "Protein", fill = T)
    }) %>% rbindlist(idcol = "Islet_Number")

x3 <- x2 %>% pivot_longer(cols = -c(0,1,2,3), names_to = "NewCol2", values_to = "mhDist") %>% filter(!is.na(mhDist))

###
x <- pData(m1) %>%
  split.data.frame(.$Islet_Number) %>% 
  lapply(function(xi)
  {
    #xi$ix_coord <- xi[xi$Islet_status == "Islet", "Xcoord"]
    #xi$iy_coord <- xi[xi$Islet_status == "Islet", "Ycoord"]
    xi$mhD <- abs(xi$Xcoord - xi[xi$Islet_status == "Islet", "Xcoord"]) + abs(xi$Ycoord - xi[xi$Islet_status == "Islet", "Ycoord"]) 
    
    
    return(xi)
  }) %>% rbindlist(idcol = "Islet_Number")


tmp <- cbind(pData(m1)[,c("NewCol", "Islet_Number"), drop = FALSE], t(exprs(m1))) %>% 
  pivot_longer(cols = -c("NewCol","Islet_Number"), names_to = "Protein", values_to = "exprs") %>% #filter(!is.na(exprs)) %>% 
  #mutate(NewCol = factor(NewCol, levels = as.character(0:3))) %>% 
  as.data.frame()

tmp2 <- tmp %>% 
  #dplyr::mutate(NewCol = factor(NewCol, levels = as.character(0:3))) %>%
  dplyr::group_by(Islet_Number, Protein) %>% 
  dplyr::summarise(corVal = cor(NewCol, exprs, method = "pearson", use = "pairwise.complete.obs"), std.dev = sd(exprs, na.rm = T))

tmp2 <- tmp %>% 
  #dplyr::mutate(NewCol = factor(NewCol, levels = as.character(0:3))) %>%
  dplyr::group_by(Protein) %>% 
  dplyr::summarise(corVal = cor(NewCol, exprs, method = "pearson", use = "pairwise.complete.obs"), std.dev = sd(exprs, na.rm = T))

#don't worry about for now
tmp3 <- tmp2 %>% 
  group_by(Islet_Number) %>% 
  pivot_wider(id_cols = Protein) %>% as.matrix() %>% 
  dplyr::summarise(cor2 = cor(., method = "pearson"))




