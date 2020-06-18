install.packages("remotes")
remotes::install_github("constantAmateur/SoupX")
library(remotes)
library(SoupX)
library(Seurat)

#Load CellRanger output
sc = load10X('LUNGCELL1/outs')

#Set contamination fraction (10 or 20%)
sc = setContaminationFraction(sc, 0.2)

#Adjust count matrix (should automatically use clustering data from CellRanger)
out = adjustCounts(sc, roundToInt=T, nCores = 20)


#Generate Seurat object for merge
set.seed(1)
Cell1twentypct = CreateSeuratObject(out)
Cell1zeropct = CreateSeuratObject(scCELLS1wIntrons$toc)
Nuc1NoSoup[["percent.mt"]] <- PercentageFeatureSet(Nuc1NoSoup, pattern = "^mt-")
Nuc1NoSoup
mean(Cell1zeropct@meta.data$nFeature_RNA)
mean(Cell1zeropct@meta.data$nCount_RNA)
Nuc1NoSoup = subset(Nuc1NoSoup, subset = nFeature_RNA > 500 & nFeature_RNA < 2500)
Nuc1NoSoup = NormalizeData(Nuc1NoSoup)
Nuc1NoSoup = FindVariableFeatures(Nuc1NoSoup)
Nuc1NoSoup = ScaleData(Nuc1NoSoup)
Nuc1NoSoup = RunPCA(Nuc1NoSoup)
ElbowPlot(Nuc1NoSoup, ndims = 50)
Nuc1NoSoup = RunUMAP(Nuc1NoSoup, dims = 1:35)
Nuc1NoSoup = FindNeighbors(Nuc1NoSoup,dims = 1:35)
Nuc1NoSoup = FindClusters(Nuc1NoSoup, resolution = 0.8, dims = 1:35)
DimPlot(SNCombzCl4, label = T, reduction = "umap") + NoLegend()

#Seurat2: sratFrame$Cluster = factor(srat@active.ident[rownames(sratFrame),'res.1'])
sratFrame$Cluster<-srat@active.ident

#Save object
saveRDS(Nuc1NoSoup, "Nuc1noSoupSeurat.rds")

#Soup Object Harmony Merge
library(harmony)
library(dplyr)
library(magrittr)
library(harmony)

srat@meta.data$stim <- "0%"
Cell1TenPctSeurat@meta.data$stim <- "10%"
Cell1twentypct@meta.data$stim <- "20%"

Nuc1NoSoup@meta.data$stim <- "0%"
Nuc1TenPctSeurat@meta.data$stim <- "10%"
Nuc1TwentyPctSeurat@meta.data$stim <- "20%"

#Merger of nucleus soup objects (same procedure for cell objects)
NucSoupMerge <- merge(x = Nuc1NoSoup, y = list(Nuc1TenPctSeurat, Nuc1TwentyPctSeurat), add.cell.ids = c("0%", "10%", "20%"))
SNCombzCl2

NucSoupMerge %<>% Seurat::NormalizeData()

NucSoupMerge %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 
NucSoupMerge %<>% ScaleData()
NucSoupMerge %<>% RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE)
NucSoupMerge <- NucSoupMerge %>% 
  RunHarmony("stim", plot_convergence = TRUE)
ElbowPlot(NucSoupMerge)
NucSoupMerge <- NucSoupMerge %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

DimPlot(NucSoupMerge, group.by = "stim", label = T, label.size = 3) + NoLegend()
FeaturePlot(NucSoupMerge, c("Rgs6"), split.by = "stim", order =T)
FeaturePlot(NucSoupMerge, c("Myh6"), order = T)

saveRDS(NucSoupMerge, file = "NucSoupMerge.rds")
write.csv(NucSoupMergeMarkers, file = "NucSoupMergeMarkers.csv")
NucSoupMergeMarkers <- FindAllMarkers(NucSoupMerge, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Remove suspected doublet clusters
nSoupMergeCl <- subset(NucSoupMerge, idents = c("0","1","2","3","4","5","6","8","9","10","11","12","13","14","15","16","17",
                                                "19","24","37","47","49","52","54"))


