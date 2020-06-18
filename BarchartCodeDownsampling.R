library(ggplot2)
library(reshape2)
library(Seurat)

#Object downsampling...
Idents(NWCellNucCL) <- "stim2"
NWCellNucDown <- subset(NWCellNucCL, idents = c("Cell","Nucleus"), downsample = 10000)
table(NWCellNucDown@active.ident)
saveRDS(NWCellNucDown, file = "NWCellNucDown.rds")

Idents(NWCellNucCL) <- "stim"
NWCellNucDownReps <- subset(NWCellNucCL, idents = c("NWCt1","NWCt2","CTRL1","CTRL2"), downsample = 6000) #round to largest #1000s that all groups have
table(NWCellNucDownReps@active.ident)
NWCellNucDownReps@active.ident <- factor(NWCellNucDownReps@active.ident, levels = c("NWCt1","NWCt2","CTRL1","CTRL2"))
saveRDS(NWCellNucDownReps, file = "NWCellNucDown.rds")
levels(NWCellNucDownReps@active.ident)
NWCellNucDownReps@meta.data$Idents2 <- NWCellNucDownReps@active.ident

Idents(NWCellNucCLTypes) <- "stim2"
NWCellNucTypesDown <- subset(NWCellNucCLTypes, idents = c("Cell","Nucleus"), downsample = 10000)
table(NWCellNucTypesDown@active.ident)
NWCellNucTypesDown@meta.data$Idents <- factor(NWCellNucTypesDown@meta.data$Idents, levels = c("Epi","Imm","Mes","EC"))
saveRDS(NWCellNucTypesDown, file = "NWCellNucTypesDown.rds")

Idents(NCCombzCLTypes) <- "stim2"
NCCombzTypesDown <- subset(NCCombzCLTypes, idents = c("Cell","Nucleus"), downsample = 10000)
table(NCCombzTypesDown@active.ident)
NCCombzTypesDown@meta.data$Idents <- factor(NCCombzTypesDown@meta.data$Idents, levels = c("Epi","Imm","Mes","EC"))
saveRDS(NCCombzTypesDown, file = "NWCellNucTypesDown.rds")

NCCombzCL@meta.data$Idents <- NCCombzCL@active.ident
Idents(NCCombzCL) <- "stim2"
NCCombzCLDown <- subset(NCCombzCL, idents = c("Cell","Nucleus"), downsample = 10000)
saveRDS(NCCombzCLDown, file = "NCCombzCLDown.rds")

Idents(NCCombzCL) <- "stim"
table(Idents(NCCombzCL))
NCCombzCLDownReps <- subset(NCCombzCL, idents = c("Cell","Nuc1","Nuc2"), downsample = 7000)

SNCombzCl3@active.ident <- SNCombzCl3@meta.data$Idents
Idents(SNCombzCl3) <- "stim"
SNCombzCl3Down <- subset(SNCombzCl3, idents = c("CTRL1","CTRL2"), downsample = 6500)
table(Idents(SNCombzCl3Down))
SNCombzCl3Down@active.ident <- SNCombzCl3Down@meta.data$Idents


#basic plot of clusters by replicate
MyPlot <- ggplot(NCCombzCLDown@meta.data, aes(x=Idents, fill=stim)) + geom_bar(width = 0.5) + scale_y_continuous(expand = c(0, 0))
MyPlot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12, face="bold"),
                            axis.text.x = element_text(angle = 45, vjust=1,hjust=1))
?geom_bar

table(NCCombzCL@meta.data$stim2)
NCCombzCL$saved.idents <- Idents(NCCombzCL)
NCCombzCL@active.ident <- NCCombzCL@meta.data$stim

Idents(NCCombzCL) <- "stim2"

#set cluster orders to match...
levels(NWCellNucCL@active.ident)
my_levels <- c("AT2","BC","EC","TC","FB1","AM","AT1","cEC","Mono1","Club/Gob","NK", 
               "avEC", "Cil", "PMN", "IM", "Peri", "Ccl17DC", "Mes","Cd103DC",
               "SMCmyoFB","LEC", "Div")
factor(object@meta.data$res.1, levels= my_levels)
object@meta.data$res.1 <- factor(object@meta.data$res.1, levels= my_levels)
ggplot(object@meta.data, aes(x=res.1, fill=replicate)) + geom_bar()

saveRDS(NCCombzCLTypes, file = "NCCombzCLTypes.rds")

NWCellNucCL@meta.data$Idents <- factor(NWCellNucCL@meta.data$Idents, levels = my_levels)
NWCellNucCLTypes@meta.data$Idents <- factor(NWCellNucCLTypes@meta.data$Idents, levels = c("Epi","Imm","Mes","EC"))
NCCombzCLTypes@meta.data$Idents <- factor(NCCombzCLTypes@meta.data$Idents, levels = c("Epi","Imm","Mes","EC"))


#plot as proportion or percentage of cluster
MyPlot <- ggplot(NCCombzCLDown@meta.data, aes(x=Idents, fill=stim)) + geom_bar(position = "fill", width = 0.5) + scale_y_continuous(expand = c(0, 0))
MyPlot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                           axis.text=element_text(size=20, face="bold"), axis.text.x = (angle = 45), vjust=1,hjust=1)

MyPlot <- ggplot(NWCellNucTypesDown@meta.data, aes(x=Idents, fill=stim2)) + geom_bar(position = "fill", width = 0.5) + scale_y_continuous(expand = c(0, 0))
MyPlot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                            axis.text.x=element_text(size=18, face="bold",angle = 45, vjust=1,hjust=1),
                            axis.text.y=element_text(size=18, face="bold"))

pg <- ggplot_build(MyPlot)
pg2 <- layer_data(MyPlot, 1)
write.csv(pg2, file = "Fig1Dnwdata.csv")

#Pull cells/nuclei separately for bar plot

NCCombzCLb <- subset(NCCombzCL, 
                     idents = c("Cell")) 
NCCombzCLb@active.ident <- NCCombzCLb$Idents

NCCombzCLc <- subset(NCCombzCL, idents = c("Nucleus"))
NCCombzCLc@active.ident <- NCCombzCLc$Idents

df2 <- as.data.frame(prop.table(table(Idents(NCCombzCLc))))
write.csv(df2, file = "nucsubtypes.csv")
NucCellSubs2<-read.csv("NucCellSubs2.csv")
data.m <- melt(NucCellSubs2, id.vars='CellType')
data.m$CellType <- factor(data.m$CellType, levels=NucCellSubs$CellType)
p1 <- ggplot(data.m, aes(CellType, value)) + geom_bar(aes(fill = variable), 
                                              width = 0.5, position = position_dodge(width=0.5), stat="identity") +
                                              scale_y_continuous(expand = c(0, 0)) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12, face="bold"),
              axis.text.x = element_text(angle = 45, vjust=1,hjust=1))
p1

#Now same for NW merged data (Reyfman et al.)
Idents(NWCellNucCL) <- "stim2"
NWCellNucCLb <- subset(NWCellNucCL, idents = "Cell")
NWCellNucCLb@active.ident <- NWCellNucCLb$Idents
df3 <- as.data.frame(prop.table(table(Idents(NWCellNucCLb))))
write.csv(df3, file = "cellsubtypesNW.csv")

NWCellNucCLc <- subset(NWCellNucCL, idents = "Nucleus")
NWCellNucCLc@active.ident <- NWCellNucCLc$Idents
df4 <- as.data.frame(prop.table(table(Idents(NWCellNucCLc))))
write.csv(df4, file = "nucsubtypesNW.csv")
#combined in excel to one spreadsheet
NucCellSubsNW <- read.csv("NucCellSubsNW.csv")
data.g <- melt(NucCellSubsNW, id.vars='CellType')
data.g$CellType <- factor(data.g$CellType, levels=NucCellSubsNW$CellType)
p2 <- ggplot(data.g, aes(CellType, value)) + geom_bar(aes(fill = variable), 
                                                      width = 0.5, position = position_dodge(width=0.5), stat="identity") +
  scale_y_continuous(expand = c(0, 0)) + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text=element_text(size=12, face="bold"),
        axis.text.x = element_text(angle = 45, vjust=1,hjust=1))
p2




CellType <- NucCellSubs$CellType
CellType

my_levels <- c("AT2","BC","EC","TC","FB1","AM","AT1","cEC","Mono1","Club/Gob","NK", 
               "avEC", "Cil", "PMN", "IM", "Peri", "Ccl17DC", "Mes","Cd103DC",
               "SMCmyoFB","LEC", "Div")


new.cluster.ids <- c("EC","Epi","Imm","Mes","Epi","EC","Epi","Imm","Epi","Mes","Imm",
                     "EC","Imm","Imm","EC","Epi","Imm","Imm","Imm","Mes","Epi","Imm",
                     "Imm","Imm","Mes","EC","Imm","EC","J","J","Mes","Imm","J","EC","J",
                     "J")
names(new.cluster.ids) <- levels(NCCombz2)
NCCombz2 <- RenameIdents(NCCombz2, new.cluster.ids)
DimPlot(NCCombz2, label = T, label.size = 3) + NoLegend()
NCCombz3 <- subset(NCCombz2, idents = c("EC","Epi","Imm","Mes"))

NCCombz2@meta.data$Idents <- NCCombz2@active.ident
NCCombz3@meta.data$Idents <- NCCombz3@active.ident
saveRDS(NCCombz3, file = "NCCombzCellCat.rds")

table(Idents(NWCellNucCL))
df <- as.data.frame(prop.table(table(Idents(NCCombzCLb))))
p3 <- ggplot(df, aes(x=Var1, y=Freq, fill=Var1)) + 
  barplot() 
p3
head(df)
write.csv(df, file= "cellsubtypes.csv")


#Group cell types as endo, immune, mesench, or epithelial

DimPlot(NCCombz2, label = T, label.size = 3) + NoLegend()
readRDS("NCCombz2.rds")
new.cluster.ids <- c("EC","AT2","AM","Col13+FB","AT1","cEC","Club","TC1","Cil","Peri",
                     "cMono","vEC","BC","ncMono","aEC","Mes","Th17TC","Ccl22+DC","Cd103+DC",
                     "Col14+FB","at2/ecDB","IM","GB","Th2TC","MyoFB","cEC","PMN","LEC","imm/mesDB",
                     "Div","SMC","pDC","DB","EC","DB","DB") #DB are doublets, remove in clean object
names(new.cluster.ids) <- levels(NCCombz2)
NCCombz2 <- RenameIdents(NCCombz2, new.cluster.ids)
NCCombz3 <- RenameIdents(NCCombz3, new.cluster.ids)

NCCombz3@active.ident <- NCCombz3@meta.data$seurat_clusters
new.cluster.ids <- c("EC","Imm","Epi","Mes","Imm","Imm","Epi","EC","Imm","Epi",
                     "Imm","Epi","Mes","Imm","EC","Imm","EC","Imm","Epi","Imm",
                     "Imm","Imm","Epi","Mes","EC","EC","Mes","Imm","Mes","Epi")
names(new.cluster.ids) <- levels(NCCombz3)
NCCombz3 <- RenameIdents(NCCombz3, new.cluster.ids)
DimPlot(NCCombz2, label = T, label.size=3) + NoLegend()


DimPlot(NWCellNuc, label = T, label.size=3, group.by = "stim") + NoLegend()
NCCombzCL <- subset(NCCombz2, idents = c("EC","AT2","AM","Col13+FB","AT1","cEC","Club","TC1","Cil","Peri",
                                         "cMono","vEC","BC","ncMono","aEC","Mes","Th17TC","Ccl22+DC","Cd103+DC",
                                         "Col14+FB","IM","GB","Th2TC","MyoFB","PMN","LEC","Div","SMC","pDC")) #remove suspected doublet clusters
saveRDS(NCCombzCL, file = "NCCombzCL.rds")
NCCombzCL@meta.data$Idents <- NCCombzCL@active.ident
NCCombzCL

#Set idents -> general types for barchart

new.cluster.ids <- c("EC","Epi","Imm","Mes","Epi","EC","Epi","Imm","Epi","Mes","Imm","EC","Imm",
                     "Imm","EC","Epi","Imm","Imm","Imm","Mes","Imm","Imm","Imm","Mes","Imm","EC","NA","Mes","Imm")
names(new.cluster.ids) <- levels(NCCombzCL)
NCCombzCL <- RenameIdents(NCCombzCL, new.cluster.ids)
NCCombzCLTypes <- subset(NCCombzCL, idents = c("Epi","Imm","EC","Mes"))
NCCombzCLTypes@meta.data$Idents <- NCCombzCLTypes@active.ident

new.cluster.ids <- c("EC","AT2","AM","Col13+FB","AT1","cEC","Club","TC1","Cil","Peri",
  "cMono","vEC","BC","ncMono","aEC","Mes","Th17TC","Ccl22+DC","Cd103+DC",
  "Col14+FB","IM","GB","Th2TC","MyoFB","PMN","LEC","Div","SMC","pDC")
names(new.cluster.ids) <- levels(NCCombzCL)
NCCombzCL <- RenameIdents(NCCombzCL, new.cluster.ids)


#Switch uMAP colors to match NCCombz object

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=36)
write.csv(color_list, file = "color_list.csv")

DimPlot(SCL1zCl4, label = T)
FeaturePlot(SCL1zCl4, c("Mreg"), order = T)

DimPlot(SNCombzCl3, label = F, cols = my_color_palette) + NoLegend()

p <- DimPlot(NCCombzCL)
pbuild <- ggplot2::ggplot_build(p)
pdata <- pbuild$data[[1]]

pdata <-  pdata[order(pdata$group), ] # Order the plot data by group
ucols <- unique(pdata$colour) # Get a vector of unique colors
names(ucols) <- unique(pdata$group)


my_color_palette <- c("#EF7F48","#F8766D","#C69900","#D69100","#9EA700","#60B200","#B4A000",
                   "#E48800","#FF64B0","#2BB600","#19B700","#83AD00","#00C088","#00BECF",
                   "#00C0BA","#00B9E3","#00BB48","#00ABFD","#00B9E1","#FF67A6","#B086FF",
                   "#46A0FF","#FC61D4","#E46DF6","#CE79FF","#FE6E8B")

my_color_palette_cell <- c("#E48800","#F8766D","#B4A000","#19B700","#EF7F48","#00BB48",
                           "#D69100","#00C1A2","#83AD00","#00ABFD","#00C088","#00C0BA",
                           "#C69900","#00B3F1","#FF64B0","#B086FF","#46A0FF","#8794FF",
                           "#00B9E1","#FF699C","#FF62BE","#60B200","#E46DF6","#FE6E8B",
                           "#F365E6","#FF62BE","#FC61D4","#00BECF","#FF67A6","#FD6F86")

my_color_palette_n1 <- c("#F8766D","#EF7F48","#D69100","#C69900","#B4A000","#19B700",
                         "#83AD00","#00C0BA","#46A0FF","#8794FF","#00BE6B","#E48800",
                         "#60B200","#9EA700","#00C088","#FF67A6","#00BECF","#00ABFD",
                         "#60B200","#FF699C","#FC61D4","#00B9E1","#B086FF","#00BB48")

my_color_palette_n2 <- c("#EF7F48","#F8766D","#D69100","#9EA700","#C69900","#60B200",
                         "#E48800","#83AD00","#00C088","#00BECF","#B4A000","#00C0BA",
                         "#19B700","#00BE6B","#8794FF","#00ABFD","#00BB48","#FF67A6",
                         "#FC61D4","#46A0FF","#00B3F1","#FD6F86")


SNucLung2zCl2 <- subset(SNucLung2zCl, idents = c("AT2","EC1","FB1","Club","AT1","Cil","AM",
                                                 "TC","BC","Mes","Kdr+ EC","EC2","Peri","EC3",
                                                 "IM","DC","Mono","SMC","LEC","FB2","IMM?","NEC"))


