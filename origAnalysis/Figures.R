library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)

#BP for Biological Process, MF for Molecular Function, and CC for Cellular Component

### Barplots of top +/- 10
BP_path <- "combinedIslets_Results_New/Proximal_Islet/gsea/gsea_results_BP.csv"
MF_path <- "combinedIslets_Results_New/Proximal_Islet/gsea/gsea_results_MF.csv"
CC_path <- "combinedIslets_Results_New/Proximal_Islet/gsea/gsea_results_CC.csv"

gsea <- read.csv(CC_path)

gsea_pos <- gsea %>% filter(NES > 0) %>% arrange(p.adjust) %>% slice(1:10)
gsea_neg <- gsea %>% filter(NES < 0) %>% arrange(p.adjust) %>% slice(1:10)

gsea_sub <- gsea_pos %>% full_join(gsea_neg)

pdf(file = "combinedIslets_Results_New/Proximal_Islet/gsea/CC_top20_pvalue.pdf", width = 11, height = 8.5)
pdf(file = "combinedIslets_Results_New/Proximal_Islet/gsea/CC_top20.pdf", width = 11, height = 8.5)

ggplot(data = gsea_sub, aes(x = reorder(Description, +NES), y = NES, fill = NES < 0)) +
  geom_col() +
  coord_flip() +
  theme(legend.position="none") +
  xlab("Cellular Component")

ggplot(data = gsea_sub, aes(x = reorder(Description, +NES), y = NES, fill = log10(pvalue))) +
  geom_col() +
  coord_flip() +
  xlab("Cellular Component")

dev.off()

#Combine multiple databases
gseaBP <- read.csv(BP_path)
gsea_posBP <- gseaBP %>% filter(NES > 0) %>% arrange(p.adjust) %>% slice(1:10)
gsea_negBP <- gseaBP %>% filter(NES < 0) %>% arrange(p.adjust) %>% slice(1:10)
gsea_subBP <- gsea_posBP %>% full_join(gsea_negBP) %>% mutate(Group = "BP")

gseaMF <- read.csv(MF_path)
gsea_posMF <- gseaMF %>% filter(NES > 0) %>% arrange(p.adjust) %>% slice(1:10)
gsea_negMF <- gseaMF %>% filter(NES < 0) %>% arrange(p.adjust) %>% slice(1:10)
gsea_subMF <- gsea_posMF %>% full_join(gsea_negMF) %>% mutate(Group = "MF")

gseaCC <- read.csv(CC_path)
gsea_posCC <- gseaCC %>% filter(NES > 0) %>% arrange(p.adjust) %>% slice(1:10)
gsea_negCC <- gseaCC %>% filter(NES < 0) %>% arrange(p.adjust) %>% slice(1:10)
gsea_subCC <- gsea_posCC %>% full_join(gsea_negCC) %>% mutate(Group = "CC")

gsea_sub <- gsea_subBP %>% full_join(gsea_subMF) %>% full_join(gsea_subCC)

pdf(file = "combinedIslets_Results_New/Proximal_Islet/gsea/GSEA_plot.pdf", width = 11, height = 8.5)

ggplot(data = gsea_sub, aes(x = interaction(reorder(Description, +NES), Group), y = NES, fill = Group)) +
  scale_x_discrete("",breaks=interaction(gsea_sub$Description,gsea_sub$Group),labels=gsea_sub$Description) +
  geom_col() +
  coord_flip()

dev.off()

###Add approach to combine multiple gsea results
library(dplyr)
library(ggplot2)
gsea <- read.csv("combinedIslets_Results_New/correlationRank/islet_gsea.csv") %>% 
  filter(ONTOLOGY == "MF") %>% 
  split.data.frame(.$Islet_Number) %>% 
  lapply(function(xi)
    {
      xi_pos <- xi %>% filter(NES > 0) %>% arrange(p.adjust) %>% slice(1:20)
      xi_neg <- xi %>% filter(NES < 0) %>% arrange(p.adjust) %>% slice(1:20)
      xi_top20 <- xi_pos %>% full_join(xi_neg)
  })

gsea2 <- gsea %>% tibble::enframe(name = "Islet_Number2") %>% 
  tidyr::unnest(cols = value) %>% 
  select(-Islet_Number2) %>% 
  group_by(ID) %>% 
  filter(n() >= 4) %>% 
  arrange(NES)

ggplot(gsea2, aes(x = as.factor(Islet_Number), y = factor(Description, levels = unique(gsea2$Description)))) +
  xlab("Image Number") +
  ylab("Description") +
  geom_tile(color = "black", aes(fill = NES)) +
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  coord_fixed() + 
  facet_wrap(vars(ONTOLOGY)) +
  theme_dark()

### Heatmaps
set.seed(23)
library(ComplexHeatmap)
load("combinedIslets_Results_New/combinedMSnSet.RData")

m1_sub <- m1[sample(nrow(m1), size=50), ]


limma <- read.delim("combinedIslets_Results_New/Proximal_Islet/limma/limmaResults.txt") %>% 
  arrange(adj.P.Val) %>% slice(1:50)
m1_sub <- m1[rownames(m1) %in% limma$feature, ]

pdf(file = "combinedIslets_Results_New/Proximal_Islet/limma/top50DI_heatmap.pdf", width = 11, height = 8.5)

column_ha = HeatmapAnnotation(Image = m1_sub$Islet_Number, Status = m1_sub$Islet_status, Plex = m1_sub$Plex)
Heatmap(exprs(m1_sub), top_annotation = column_ha, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8))

dev.off()

###Volcano plot
limma <- read.delim("combinedIslets_Results_New/Proximal_Islet/limma/limmaResults.txt") %>% 
  mutate(sig = ifelse(.$P.Value<=0.05 & .$logFC<=-1, "Down", ifelse(.$P.Value<=0.05 & .$logFC>=1, "Up", "Normal")))

pdf(file = "combinedIslets_Results_New/Proximal_Islet/limma/volcanoplot.pdf", width = 11, height = 8.5)

volc = ggplot(limma, aes(logFC, -log10(P.Value))) + #log2Foldchange vs FDR
  geom_point(aes(col=sig)) + # the points are colored by significance
  scale_color_manual(values=c("green", "black", "red")) + 
  ggtitle("Proximal vs Islet - Volcano Plot") #e.g. 'Title description'
volc+
  geom_text_repel(data= limma[(limma$sig == "Down" & -log10(limma$P.Value) >= 23) | ((limma$sig == "Up" & -log10(limma$P.Value) >= 16)),], aes(label=feature),
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) +
  geom_hline(yintercept= -log10(.05),linetype=2, alpha = .4) + 
  geom_vline(xintercept= c(-1,1),linetype=2, alpha = .4)

dev.off()
