#------------------------------------------------------------------------------------
# Cell volume regulation during T lymphoblast controls T cell activation and function
#------------------------------------------------------------------------------------
#---I. Prepare work----
workdir <- "H:/CTD/COOPERATION/WangYuMan/20230919-GitHubProgram"
setwd(workdir)
library(DESeq2)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(forcats)
library(pheatmap)
library(forcats)
library(homologene)
library(KEGGREST)
library(readxl)
library(ggpubr)
library(ggpattern)
R.utils::setOption("clusterProfiler.download.method",'auto')
#---II. Input data----
cts1 <- read.delim("./Data/expression_data/gene_count_matrix.txt")
cts1 <- cts1[!duplicated(cts1$gene_id),]
cts <- cts1[,-1]
rownames(cts) <- cts1$gene_id
coldata <- data.frame(genetype = gsub("^(WT|KO).*", "\\1", names(cts)),
                      condition = gsub("^(WT|KO)[0-9]+([A-Z]+)", "\\2", names(cts)))
coldata$group <- paste(coldata$genetype, coldata$condition, sep = "_")
datExp1 <- read_xlsx("./data/expression_data/expression.xlsx")
datExp1 <- as.data.frame(datExp1[!duplicated(datExp1$SYMBOL),])
datExp <- datExp1[,-1]
rownames(datExp) <- datExp1$SYMBOL
colnames(datExp) <- gsub("^([A-Z]{2})([A-Z]{1,2}).*_[0-9]([0-9])_.*$", "\\1\\3\\2", toupper(colnames(datExp)))
#---III. DESeq2 Analysis----
## 1. DESeq2 Differential Expression Genes filter
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ group)
dds <- DESeq(dds)
## 2. Comparison of different groups
CT_compare_within_gr <- list(c("group", "WT_NH", "WT_CT"),
                             c("group", "WT_NL", "WT_CT"),
                             c("group", "WT_G", "WT_CT"),
                             c("group", "KO_NH", "KO_CT"),
                             c("group", "KO_NL", "KO_CT"),
                             c("group", "KO_G", "KO_CT"))
G_compare_within_gr <- list(c("group", "WT_G", "WT_NH"),
                            c("group", "WT_G", "WT_NL"),
                            c("group", "KO_G", "KO_NH"),
                            c("group", "KO_G", "KO_NL"))
KO_compare_between_gr <- list(c("group","KO_CT", "WT_CT"),
                              c("group","KO_G","WT_G"),
                              c("group","KO_NH", "WT_NH"),
                              c("group","KO_NL", "WT_NL"))
contrast.name <- list(CT_compare_within_gr, G_compare_within_gr, KO_compare_between_gr)
# 2.1 Comparison of different stimuli within groups
gene_counts <- list(CT_wit=c(),G_wit=c(),KO_bet=c())
res <- list(CT_wit=list(),G_wit=list(),KO_bet=list())
for (j in 1:length(contrast.name)) {
  for (i in 1:length(contrast.name[[j]])) {
    res[[j]][[i]] <- results(dds, contrast = contrast.name[[j]][[i]])
    filenames <- paste("./Data/DEGs_data/DEGs",contrast.name[[j]][[i]][2],contrast.name[[j]][[i]][3],"comparison.csv", sep = "_")
    write.csv(res[[j]][[i]], file = filenames)
    k <- (i-1) * 2 + 1
    gene_counts[[j]][k] <- res[[j]][[i]] %>% as.data.frame() %>%
      filter(padj < 0.05 & log2FoldChange > 1) %>%
      count() %>%
      as.numeric()
    gene_counts[[j]][k+1] <- res[[j]][[i]] %>% as.data.frame() %>%
      filter(padj < 0.05 & log2FoldChange < -1) %>%
      count() %>%
      as.numeric()
  }
}
#----Figure 4.a----
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd)
percentVar <- c(35.7, 15.5)
pca_plot <- pcaData$data
pca_plot$type <- c(rep("KO",12),rep("WT",12))
ggplot(pca_plot, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(pcaData$labels$x) +
  ylab(pcaData$labels$y) +
  scale_color_manual(values = c("#627f96", "#ceba75","#d34e2c", "#529054"))
ggsave("./Result/Figure4a/Figure4a.pdf", width=8, height=6)
#----Figure 4.b----
volume <- read.csv("./Data/CellVolume/cell_volume.csv")
total_counts <- apply(datExp, 2, sum)
meltp <- data.frame(cellcounts = as.numeric(total_counts),
                    group = factor(rep(c("KO","WT"), each = 12)),
                    color = factor(rep(rep(c("unstim.", "G4", "N4_High", "N4_low"),each = 3), 2)))
for (i in c(1,4)) {
  meltp$volume <- as.numeric(volume[i,-1])
  meltp %>% ggplot(mapping = aes(x = cellcounts, y = volume)) +
    geom_point(mapping = aes(x = cellcounts, y = volume, shape = group,
                             color = color), size = 5)+
    scale_color_manual(values = c(unstim. = "#4582b4",
                                  G4 = "#e7b815",
                                  N4_High = "#fc4e07",
                                  N4_low = "#3aa44e")) +
    #color = c("#4582b4", "#e7b815", "#fc4e07", "#3aa44e")), size = 5)+
    geom_smooth(method = "lm", color = "darkred")+
    stat_cor(method = "spearman")+
    ylab(volume$X[i])+
    theme_classic()
  # scale_x_break(c(0,922609), scales = "free", expand = T,space = 1)+
  # scale_x_continuous(limits = range(meltp$cellcounts))
  filename = paste("./Result/Figure4b/",volume$X[i], "_lm_plot.pdf", sep="")
  ggsave(filename = filename, width=10, height=8)
}
#----Figure 4.c----
### KO&WT inner group compare
contrast.name <- list(c("group", "WT_NH", "WT_CT"),
                      c("group", "WT_NL", "WT_CT"),
                      c("group", "WT_G", "WT_CT"),
                      c("group", "KO_NH", "KO_CT"),
                      c("group", "KO_NL", "KO_CT"),
                      c("group", "KO_G", "KO_CT"))
nCountKO <- c()
resKO <- list()
for (i in 1:length(contrast.name)) {
  resKO[[i]] <- results(dds, contrast = contrast.name[[i]])
  #filenames <- paste(contrast.name[[i]][2],contrast.name[[i]][3],"KO_DEGs.csv", sep = "_")
  #write.csv(resKO[[i]], file = filenames)
  nCountKO[i] <- resKO[[i]] %>% as.data.frame() %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
    count() %>%
    as.numeric()
}
meltp <- data.frame(Name = c("WT_NH_CT","WT_NL_CT","WT_G_CT","KO_NH_CT","KO_NL_CT","KO_G_CT"),
                    Count = nCountKO)
write.csv(meltp, file = "KO and WT nCount plot.csv")
meltp$Name <- factor(meltp$Name, levels = c("WT_NL_CT","WT_G_CT","WT_NH_CT",
                                            "KO_NL_CT","KO_G_CT","KO_NH_CT"))
meltp %>% 
  ggplot(mapping = aes(x=Name, y=Count, fill = Name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Set3")[-2])+
  theme_classic() +
  theme(axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        legend.text = element_text(size = 14))
ggsave("./Result/Figure4c/Figure4c.pdf")
#----Figure 4d----
#dir.create("./Result/Figure4d")
### WT&KO group
idx.name <- list(c("group","KO_CT", "WT_CT"),
                 c("group","KO_G","WT_G"),
                 c("group","KO_NH", "WT_NH"),
                 c("group","KO_NL", "WT_NL"))
nCount_KO_WT <- c()
res <- list()
for (i in 1:length(idx.name)) {
  res[[i]] <- results(dds, idx.name[[i]])
  filenames <- paste(idx.name[[i]][1], idx.name[[i]][2], "DEGs.csv",sep = "_")
  write.csv(res[[i]], file = filenames)
  nCount_KO_WT[i] <- res[[i]] %>%
    as.data.frame() %>%
    filter(padj < 0.01 & abs(log2FoldChange) > 1) %>%
    count() %>%
    as.numeric()
}
meltp <- data.frame(Name = c("CT","G","NH","NL"),
                    Count = nCount_KO_WT)
write.csv(meltp, file = "KO&WT DEGs nCount plot only.csv")
meltp %>% mutate(Name = fct_reorder(Name, Count)) %>%
  ggplot(meltp, mapping = aes(x=Name, y=Count, fill = Name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Set2"))+
  theme_classic() +
  theme(axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        legend.text = element_text(size = 14))
ggsave("./Result/Figure4d/Figure4d.pdf")
#----Figure 4e----
#dir.create("./Result/Figure4e")
### G compared with NL & NH
contrast.name <- list(c("group", "WT_G", "WT_NH"),
                      c("group", "WT_G", "WT_NL"),
                      c("group", "KO_G", "KO_NH"),
                      c("group", "KO_G", "KO_NL"))
nCountWT <- c()
resWT <- list()
for (i in 1:length(contrast.name)) {
  resWT[[i]] <- results(dds, contrast = contrast.name[[i]])
  filenames <- paste(contrast.name[[i]][2],contrast.name[[i]][3],"WT_DEGs.csv", sep = "_")
  write.csv(resWT[[i]], file = filenames)
  nCountWT[i] <- resWT[[i]] %>% as.data.frame() %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
    count() %>%
    as.numeric()
}
meltp <- data.frame(Name = c("WT_G_NH","WT_G_NL","KO_G_NH", "KO_G_NL"),
                    Count = nCountWT)
write.csv(meltp, file = "G compared with N nCount plot only.csv")
meltp$Name <- factor(meltp$Name, levels = c("WT_G_NL","WT_G_NH",
                                            "KO_G_NL", "KO_G_NH"))
meltp %>%
  ggplot(mapping = aes(x=Name, y=Count, fill = Name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Set3")[5:8]) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
        legend.text = element_text(size = 14))
ggsave("./Result/Figure4e/Figure4e.pdf")
#----Figure 4f----
#dir.create("./Result/Figure4f")
## 1. prepare work
idx1 <- c("WT", "KO")
idx2 <- c(1,2,3)
idx3 <- c("CT", "NL", "G", "NH")
idx.name <- c()
for (i in 1:length(idx3)) {
  for (j in 1:length(idx1)) {
    for (m in 1:length(idx2)) {
      idx.name[length(idx.name)+1] <- paste(idx1[j], idx2[m], idx3[i], sep = "")
    }
  }
}
## 2. heat map input
file_names <- c("./data/DEGs_data/DEGs_KO_CT_WT_CT_comparison.csv",
                "./data/DEGs_data/DEGs_KO_G_WT_G_comparison.csv",
                "./data/DEGs_data/DEGs_KO_NH_WT_NH_comparison.csv",
                "./data/DEGs_data/DEGs_KO_NL_WT_NL_comparison.csv")
total_DEGs_list <- c()
for (i in 1:length(file_names)) {
  k <- read.csv(file = file_names[i], header = T)
  k <- k %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1)
  total_DEGs_list <- c(total_DEGs_list, k$X)
}
total_DEGs_list <- unique(total_DEGs_list)
break_bar <- c(seq(-2.75,-0.1,by = 0.01), seq(0,3,by = 0.01))
#pdf(file = "Total_heatmap_DEGs.pdf")
(pf <- pheatmap(cts[total_DEGs_list,idx.name], scale = "row", cutree_rows = 3,
               cluster_cols = F, show_rownames = F, cluster_rows = T,
               color = c(colorRampPalette(colors = c("blue", "white"))(round(length(break_bar)/2)),
                         colorRampPalette(colors = c("white", "red"))(round(length(break_bar)/2))),
               breaks = break_bar,
               filename = "./Result/Figure4f/Figure4f.pdf",
               height = 10))
#----Figure 4g----
#dir.create("Result/Figure4g")
#filenames <- c("group_KO_CT_DEGs.csv", "group_KO_G_DEGs.csv", "group_KO_NH_DEGs.csv", "group_KO_NL_DEGs.csv")
filenames <- c("group_KO_G_DEGs.csv", "group_KO_NH_DEGs.csv", "group_KO_NL_DEGs.csv")
## 9.1 KEGG ORA Analysis
kk <- list()
for (i in 1:length(filenames)) {
  k <- read.csv(file = filenames[i], header = T)
  k <- k %>% arrange(desc(log2FoldChange)) %>% filter(pvalue < 0.05 & abs(log2FoldChange) > 1)
  genelist <- k$log2FoldChange
  names(genelist) <- AnnotationDbi::select(org.Mm.eg.db, k$X, keytype = "SYMBOL", columns = "ENTREZID")[,2]
  kk[[i]] <- enrichKEGG(gene = names(genelist), organism = "mmu", pvalueCutoff = 0.05)
}
## 9.2 Visualization the result of ORA KEGG
### 9.2.1 merge the KEGG plot
merge_list<-list()
for (i in 1:length(kk)) {
  merge_list[[i]] <- kk[[i]]@result$ID
}
pathways <- list(common = Reduce(intersect, merge_list),
                 #CT = Reduce(setdiff, merge_list),
                 G = Reduce(setdiff, list(merge_list[[1]],merge_list[[2]],merge_list[[3]])),
                 NH = Reduce(setdiff, list(merge_list[[2]],merge_list[[1]],merge_list[[3]])),
                 NL = Reduce(setdiff, list(merge_list[[3]],merge_list[[2]],merge_list[[1]])))
merge_meltp <- list()
merge_column <- c("Description", "GeneRatio", "p.adjust")
Condition <- c("G","NH","NL")
for (i in 1:length(kk)) {
  merge_meltp[[i]]<-kk[[i]]@result[,merge_column]
  merge_meltp[[i]]$Condition <- rep(Condition[i], nrow(kk[[i]]@result))
}
merge_meltp <- bind_rows(merge_meltp)
#Gene Ratio character convert into number
merge_meltp$GeneRatio <- apply(merge_meltp, 1, function(x){numerator=as.numeric(strsplit(x[2],"/")[[1]][1])
denominator = as.numeric(strsplit(x[2],"/")[[1]][2])
GeneRatio = numerator/denominator})
#Filter and arrange the data
merge_meltp %<>% filter(p.adjust < 0.05) %>% arrange(desc(GeneRatio))
#Visualization
merge_meltp$Condition <- factor(merge_meltp$Condition, levels = c("G","NH","NL"))
merge_meltp %>% mutate(Description = fct_reorder(Description, GeneRatio, .desc = T)) %>%
  ggplot(mapping = aes(x=Description, y=GeneRatio, size = -log(p.adjust), color = Condition)) +
  geom_point() +
  scale_size_continuous(range = c(4,10)) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10))
# dir.create("Ameliorate KEGG")
ggsave("./Result/Figure4g/Figure4g.pdf", width = 16, height = 9)
#----Figure 4h----
#dir.create("Result/Figure4h")
idx1 <- c("WT", "KO")
idx2 <- c(1,2,3)
idx3 <- c("CT", "G", "NH", "NL")
idx.name <- c()
for (i in 1:length(idx3)) {
  for (j in 1:length(idx1)) {
    for (m in 1:length(idx2)) {
      idx.name[length(idx.name)+1] <- paste(idx1[j], idx2[m], idx3[i], sep = "")
    }
  }
}
total_DEGs_list <- c("Nr4a1", "Nr4a2", "Nr4a3", "Cd69", "Tnf", "Itgb2", "Cd2", "Cd5",
                     "Irf4", "Irf8", "Egr2", "Il2", "Il2ra", "Cd40lg", "Pdcd1", "Tbx21", "Tnfrsf4",
                     "Ifng", "Gzmb", "Il12rb2", "Tnf",
                     "Ctla4", "Icos", "Maf",
                     "Il10", "Nfil3", "Lag3", "Tigit", "Tnfrsf4", "Ctla4", "Havcr2", "Cd244a",
                     "Cd3g", "Cd3e", "Lck", "Rasgrp1", "Gata3")
pheatmap(cts[total_DEGs_list,idx.name], scale = "row",
              cluster_cols = F, show_rownames = F, cluster_rows = T,
              color = c(colorRampPalette(colors = c("blue", "white"))(round(length(break_bar)/2)),
                        colorRampPalette(colors = c("white", "red"))(round(length(break_bar)/2))),
              breaks = break_bar,
              filename = "./Result/Figure4h/Figure4h.pdf")
#----Figure 4i----
#dir.create("Result/Figure4i")
TF_genes <- read.csv("./data/TF_analysis/TFlist.csv", header = F)
TF_symbol <- AnnotationDbi::select(org.Mm.eg.db, keys = TF_genes$V1,keytype = "ENSEMBL",columns = "SYMBOL")
TF_cts <- cts[TF_symbol$SYMBOL,]
TF_cts <- na.omit(TF_cts)
TF_exp <- datExp[TF_symbol$SYMBOL,]
TF_exp <- na.omit(TF_exp)
## 2.differential TF genes analysis
dds <- DESeqDataSetFromMatrix(countData = TF_cts,colData = coldata,design = ~ group)
dds <- DESeq(dds)
contrast.name <- list(c("group","KO_CT", "WT_CT"),
                      c("group","KO_G","WT_G"),
                      c("group","KO_NH", "WT_NH"),
                      c("group","KO_NL", "WT_NL"))
res <- list()
for (i in 1:length(contrast.name)) {
  res[[i]] <- results(dds, contrast.name[[i]])
}
## 3. heat map of the TFs
### 3.1 prepare
break_bar <- c(seq(-2.75,-0.1,by = 0.01), seq(0,3,by = 0.01))
idx1 <- c("WT", "KO")
idx2 <- c(1,2,3)
idx3 <- c("CT", "NL", "G", "NH")
idx.name <- c()
for (i in 1:length(idx3)) {
  for (j in 1:length(idx1)) {
    for (m in 1:length(idx2)) {
      idx.name[length(idx.name)+1] <- paste(idx1[j], idx2[m], idx3[i], sep = "")
    }
  }
}
### 3.2 select and filter the DEG
DEGs_TF <- c()
for (i in 1:4){
  DEGs_G_TF <- res[[i]] %>% as.data.frame() %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
  DEGs_TF <- c(DEGs_TF, rownames(DEGs_G_TF))
}
DEGs_TF <- unique(DEGs_TF)
DEGs_TF <- na.omit(DEGs_TF)
### 3.3 plot the heat map
break_bar <- c(seq(-2,-0.1,by = 0.01), seq(0,2,by = 0.01))
(p <- pheatmap(na.omit(datExp[DEGs_TF,idx.name]),scale = "row", cluster_cols = F,cutree_rows = 4,
               #cellheight = 10, cellwidth = 10,
               color = c(colorRampPalette(colors = c("blue", "white"))(round(length(break_bar)/2)),
                         colorRampPalette(colors = c("white", "red"))(round(length(break_bar)/2))),
               breaks = break_bar,
               filename = "./Result/Figure4i/Figure4i.pdf"))
#----Supplementary Figure 6.b----
#dir.create("Result/SFigure6b")
ego<-list()
kk<-list()
### 3.2 heat map idx number
cluster_genes <- c("Klra8","Gm33877","Gm32525","Spats2l")
idx <- list()
for (i in 1:length(cluster_genes)) {
  idx[[i]] <- which(pf$tree_row$labels[pf$tree_row$order] %in% cluster_genes[i])
}
genes <- list()
for (i in 1:3) {
  genes[[i]] <- pf$tree_row$labels[pf$tree_row$order][idx[[i]]:idx[[(i+1)]]]
  write.csv(genes[[i]], file = paste("./Result/SFigure6b/Q",i,".csv", sep = ""))
  genelist <- AnnotationDbi::select(org.Mm.eg.db, keys = genes[[i]], keytype = "SYMBOL", columns = "ENTREZID")[,2]
  kk[[i]]<-enrichKEGG(gene = genelist, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05)
  ego[[i]]<-enrichGO(gene = genelist, OrgDb = "org.Mm.eg.db", pvalueCutoff = 0.05)
}
egoq3 <- enrichGO(gene = genes[[3]], OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", pvalueCutoff = 0.05)
write.csv(egoq3@result,file = "./Result/SFigure6b/GOQ3-Genes.csv")
ego_p <- list()
for (i in 1:length(ego)) {
  ego_p[[i]] <- dotplot(ego[[i]])
  ego_names <- paste("./Result/SFigure6b/GO enrichment Q",i," plot.pdf", sep = "")
  pdf(ego_names)
  ego_p[[i]]
  print(ego_p[[i]])
  dev.off()
}
#----Supplementary Figure 6.c----
#dir.create("Result/SFigure6c")
idx1 <- c("WT", "KO")
idx2 <- c(1,2,3)
idx3 <- c("CT", "G", "NH", "NL")
idx.name <- c()
for (i in 1:length(idx3)) {
  for (j in 1:length(idx1)) {
    for (m in 1:length(idx2)) {
      idx.name[length(idx.name)+1] <- paste(idx1[j], idx2[m], idx3[i], sep = "")
    }
  }
}
genes <- c("Il3","Ifng","Ccl1","Lif","Ccl9","Ccl3","Xcl1","Ccl4","Il23a","Ccl6","Il2",
           "Cxcl9","Csf2","Cxcr3","Il12rb1","Il2ra","Cxcr5","Osm","Tnfsf4","Tnfsf8",
           "Cd40lg", "Ccr8","Il21","Ccr1","Tnfrsf8","Ccr5","Il1r1","Tnfsf10","Inhba",
           "Cth", "Cbs", "Asns","Gpt2","Psph","Bcat1","Psat1","Phgdh","Shmt2","Kyat3",
           "Srm","Maoa","Aoc3")
pheatmap(cts[genes,idx.name], scale = "row",
         cluster_cols = F, show_rownames = T, cluster_rows = F,
         color = c(colorRampPalette(colors = c("blue", "white"))(round(length(break_bar)/2)),
                   colorRampPalette(colors = c("white", "red"))(round(length(break_bar)/2))),
         breaks = break_bar,
         filename = "./Result/SFigure6c/SFigure6c.pdf")

