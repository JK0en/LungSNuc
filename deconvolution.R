#Deconvolution of bulk RNASeq data using sc/snRNASeq data

# devtools::install_github('shenorrlab/bseqsc')
library(bseqsc)
library(here)
library(BisqueRNA) # for Seurat to ExpressionSet
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape)
bseqsc_config('comparison_human_datasets/fan_PMID31578193/CIBERSORT.R')

#Create single cell/nuc expression set
NCCombzCL@meta.data$celltype <- NCCombzCL@active.ident
scrna_eset <- SeuratToExpressionSet(NCCombzCL, delimiter="", position=1, version="v3")
scrna_eset@phenoData$celltype <- NCCombzCL$celltype

# select top cluster specific snRNA markers for deconvolution. Vary 'slice' parameter to select different # of markers.
sc.markers <- FindAllMarkers(NCCombzCL, min.pct = 0.25, only.pos = TRUE) # only.pos=TRUE
list.markers <- sapply(levels(NCCombzCL), function(x) {
  df <- dplyr::filter(sc.markers, cluster == x)  %>%
    dplyr::filter(avg_logFC > 1 & p_val_adj < 0.01) %>%
    dplyr::arrange(desc(avg_logFC, p_val_adj)) %>%
    dplyr::slice(1:20)  %>%
    dplyr::select(gene)
})
names(list.markers) <- levels(NCCombzCL)

# build reference basis matrix of expression profiles for individual celltypes
# average counts computed within each cell type in each sample
plotCellTotals(scrna_eset, 'cellType', 'SubjectName')
B <- bseqsc_basis(scrna_eset, list.markers, clusters = 'celltype', samples = 'SubjectName', ct.scale = TRUE)
plotBasis(B, list.markers, Colv = NA, Rowv = NA, col = 'Blues')

# load the bulk eset
bulk_eset <- readRDS("comparison_mouse_datasets/eset_lung_bulk.rds")

# estimate celltype proportions with cibersort
fit <- bseqsc_proportions(eset_lung_bulk, B, verbose = TRUE)

# add estimated celltype proportions as phenoData
pData(eset_lung_bulk) <- cbind(pData(eset_lung_bulk), t(coef(fit)))

saveRDS(eset_lung_bulk, file = "eset_lung_bulkNCCombz3.rds")

#To plot cell type frequency in groups from bulk seq (young v old for mouse lung)
toplot <- eset_lung_bulk@phenoData@data %>%
  dplyr::select(group, AM)
p1 <- ggplot(toplot, aes(x=group, y=AM, fill=group)) + 
  geom_boxplot() + 
  xlab("Group") +
  ylab("Proportion of XXX Cells") +
  ggtitle("Proportion of XXX Cells in Murine IRI") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  NoLegend()
p1

##For all cell type frequencies in bulk data 

toplot4 <- eset_lung_bulk@phenoData@data
toplot4 <- toplot4[5:ncol(toplot4)]
toplot4 <- toplot4[-c(1,2,3), ] #remove data from old mice
toplot4 <- melt(toplot4)

p3 <- ggplot(toplot4, aes(x=variable, y=value, fill=variable)) + 
  geom_boxplot() + 
  ylab("Proportion of Cells") +
  ggtitle("Proportion of Cells in Bulk Data")
p3 + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(), axis.line = element_line(colour="black"),
                        axis.text = element_text(size=12, face="bold"), axis.text.x=element_text(angle=45,vjust=1,hjust=1))
