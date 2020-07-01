library(Seurat)
library(dplyr)
library(Matrix)

#Start with Seurat object from intron+exon zUMIs pipeline (nuc1, nuc2, cell1)
SNucLung1z <- subset(SNucLung1z, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
SNucLung1z <- FindVariableFeatures(SNucLung1z, selection.method = "vst", nfeatures = 2000)
SNucLung1z <- RunPCA(SNucLung1z, features = VariableFeatures(object = SNucLung1z))
ElbowPlot(SNucLung1z)
SNucLung1z <- FindNeighbors(SNucLung1z, dims = 1:40)
SNucLung1z <- FindClusters(SNucLung1z, resolution = 1.0)
SNucLung1z <- RunUMAP(SNucLung1z, dims = 1:40)
DimPlot(SNucLung1z, reduction = "umap")

SNucLung2z <- subset(SNucLung2z, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
SNucLung2z <- FindVariableFeatures(SNucLung2z, selection.method = "vst", nfeatures = 2000)
SNucLung2z <- RunPCA(SNucLung2z, features = VariableFeatures(object = SNucLung2z))
ElbowPlot(SNucLung2z)
SNucLung2z <- FindNeighbors(SNucLung2z, dims = 1:40)
SNucLung2z <- FindClusters(SNucLung2z, resolution = 1.0)
SNucLung2z <- RunUMAP(SNucLung2z, dims = 1:40)
DimPlot(SNucLung2z, label = T, reduction = "umap")

SNucLung2z.markers <- FindAllMarkers(SNucLung2z, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.10)
write.csv(SNucLung2z.markers, file = "SNucLung2zMarkers.csv")

new.cluster.ids <- c("AT2","EC1","FB1","Club", "AT1", "soup","Cil","AM", "TC", "BC",
                     "Mes","cEC","EC2","Peri","EC3","IM","DC","Mono","AT2b", "Club2",
                     "Club3","Cil2","SMC","LEC","FB2","AT1b","AT1c","AT1d","Ccl22+DC",
                     "FB3","NEC")
#Remove doublet clusters...
SNucLung2zCl <- subset(SNucLung2z, idents = c("AT2","EC1","FB1","Club", "AT1","Cil","AM", "TC", "BC",
                                              "Mes","cEC","EC2","Peri","EC3","IM","DC","Mono",
                                              "Cil2","SMC","LEC","FB2","Ccl22+DC",
                                              "FB3","NEC"))

DimPlot(SNucLung2zCl, label = T, reduction = "umap")

#For cells include mitochondrial parameter...
SCellLung1z[["percent.mt"]] <- PercentageFeatureSet(SCellLung1z, pattern = "^mt-")
SCellLung1z <- subset(SCellLung1z, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
SCellLung1z <- FindVariableFeatures(SCellLung1z, selection.method = "vst", nfeatures = 2000)
SCellLung1z <- RunPCA(SCellLung1z, features = VariableFeatures(object = SCellLung1z))
ElbowPlot(SCellLung1z)
SCellLung1z <- FindNeighbors(SCellLung1z, dims = 1:40)
SCellLung1z <- FindClusters(SCellLung1z, resolution = 1.0)
SCellLung1z <- RunUMAP(SCellLung1z, dims = 1:40)
DimPlot(SCellLung1z, reduction = "umap")




#Single cell immune clustering for comparison to snuc
SCL1zImmune <- subset(SCL1zCl4, idents = c("AM","Mono1","Mono2","BC","DC","DC2","TC","TC2","NK","IM","Nphil","GB","pDC","Baso"))
SCL1zImmune <- FindVariableFeatures(SCL1zImmune, selection.method = "vst", nfeatures = 2000)
SCL1zImmune <- ScaleData(SCL1zImmune, features = rownames(SCL1zImmune))
DimPlot(SCL1zCl4, label = T, label.size = 3) + NoLegend()
SCL1zImmune <- RunPCA(SCL1zImmune, features = VariableFeatures(object = SCL1zImmune))
ElbowPlot(SCL1zImmune, ndims = 50)
SCL1zImmune <- FindNeighbors(SCL1zImmune, dims = 1:35)
SCL1zImmune <- FindClusters(SCL1zImmune, resolution = 0.7)
SCL1zImmune <- RunUMAP(SCL1zImmune, dims = 1:35)
DimPlot(SCL1zCl4, reduction = "umap", label = T, label.size = 3) + NoLegend()
FeaturePlot(SCL1zImmune, features = "Il12b", order = T)

SCL1zImmune.markers <- FindAllMarkers(SCL1zImmune, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.10)
write.csv(SCL1zImmuneCl2.markers, file = "SCL1zImmuneCl2Markers.csv")

new.cluster.ids <- c("AM","cMono","TC","ncMono", "BC", "GB","Cd209+ DC","Mmp12+ AM (M2)", "Cd103+ DC",
                     "NK","Nphil","IM","unk1","unk2","pDC","Ccl22+ DC","Baso/mast","Div",
                     "Il12b/Ccl22+","Baso")
                    
names(new.cluster.ids) <- levels(SCL1zImmuneCl2)
SCL1zImmuneCl2 <- RenameIdents(SCL1zImmuneCl2, new.cluster.ids)

SCL1zImmuneCl2 <- subset(SCL1zImmuneCl, idents = c("AM","TC","cMono","ncMono","Cd209+ DC","BC","GB","Cd103+ DC",
                                                "Mmp12+ AM","NK","Nphil","IM","?lymph","?lymph2","Il12a+?",
                                                "Il12b/Ccl22+","Baso"))
SCL1zImmuneCl2 <- FindVariableFeatures(SCL1zImmuneCl2, selection.method = "vst", nfeatures = 2000)
SCL1zImmuneCl2 <- ScaleData(SCL1zImmuneCl2, features = rownames(SCL1zImmuneCl2))
SCL1zImmuneCl2 <- RunPCA(SCL1zImmuneCl2, features = VariableFeatures(object = SCL1zImmuneCl2))
ElbowPlot(SCL1zImmuneCl2, ndims = 50)
SCL1zImmuneCl2 <- FindNeighbors(SCL1zImmuneCl2, dims = 1:30)
SCL1zImmuneCl2 <- FindClusters(SCL1zImmuneCl2, resolution = 0.7)
SCL1zImmuneCl2 <- RunUMAP(SCL1zImmuneCl2, dims = 1:30)
DimPlot(SCL1zImmuneCl2, reduction = "umap", label = F, label.size = 3) + NoLegend()
FeaturePlot(SCL1zImmuneCl2, features = "Ccl22", order = T)

SCL1zImmuneCl3 <- subset(SCL1zImmuneCl2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
DimPlot(SCL1zCl5, label = T, label.size = 3) + NoLegend()
SCL1zImmuneCl2
SCL1zImmuneCl2.markers <- FindAllMarkers(SCL1zImmuneCl2, only.pos = T, min.pct = 0.10, logfc.threshold = 0.1)

saveRDS(SCL1zImmuneCl3, file = "SCL1zImmuneCl3.rds")

