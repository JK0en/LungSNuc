library(SeuratData)
library(Seurat)
library(dplyr)
library(harmony)

#Assign labels to Seurat objects

SNucLung1z@meta.data$stim3 <- "Nuc1"
SNucLung2zCl@meta.data$stim3 <- "Nuc2"
SCL1zCl4@meta.data$stim3 <- "Cell"

#Single nuc data merge

SNCombz <- merge(SNucLung1z, y = SNucLung2zCl, add.cell.ids = c("Nuc1","Nuc2"))
SNCombz <- NormalizeData(SNCombz) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
SNCombz <- RunHarmony(SNCombz, group.by.vars = "stim3")
ElbowPlot(SNCombz, ndims = 50)
SNCombz <- RunUMAP(SNCombz, reduction = "harmony", dims = 1:40)
SNCombz <- FindNeighbors(SNCombz, reduction = "harmony", dims = 1:40) %>% FindClusters()
DimPlot(SNCombz, group.by = c("stim3"), ncol = 3)

#Cell and nuc data merge

NCCombz <- merge(SCL1zCl4, y = c(SNucLung1z, SNucLung2zCl), add.cell.ids = c("Cell","Nuc1","Nuc2"), project = "NCCombz")
NCCombz <- NormalizeData(NCCombz) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
NCCombz <- RunHarmony(NCCombz, group.by.vars = "stim3")
ElbowPlot(NCCombz, ndims = 50)
NCCombz <- RunUMAP(NCCombz, reduction = "harmony", dims = 1:40)
NCCombz <- FindNeighbors(NCCombz, reduction = "harmony", dims = 1:40) %>% FindClusters()
DimPlot(NCCombz, group.by = c("stim3"), ncol = 3)

