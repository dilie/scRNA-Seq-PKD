library(tidyverse)
library(Seurat)
library(harmony)
library(cowplot)

setwd("~/Desktop/local-scRNA-PKD/")

x = readRDS("1-integrated.rds")
#x$new <- ifelse(x$orig.ident == "Cont", "Ctrl", "PKD")
#x$orig.ident = NULL

# remove PKD2
x =subset(x, sample!="PKD2")


# Linear dimensional reduction (PCA):

# variable genes
# method 1
#VariableFeatures(x) <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)

# method 2
VariableFeatures(x) = split(colnames(x), x$sample) %>% lapply(function(cells_use) {
#VariableFeatures(x) = split(colnames(x), x$orig.ident) %>% lapply(function(cells_use) {
    x[,cells_use] %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    VariableFeatures()
}) %>% unlist %>% unique

# method 3: good
wt <- subset(x, idents = "Cont")
VariableFeatures(x) = FindVariableFeatures(wt, selection.method = "vst", nfeatures = 2000) %>% VariableFeatures()
rm(wt)

x <- ScaleData(x)
x = x  %>% RunPCA(npcs = 30)

#print(x[["pca"]], dims = 1:4, nfeatures = 5)
# first 4 principal components; for each PC, find top 5 genes contributing to the PC

#VizDimLoadings(x, dims = 1:2, reduction = "pca")

DimPlot(x, reduction = "pca", pt.size = 0.01)

ElbowPlot(x) # Determine the ‘dimensionality’ of the dataset for UMAP, heatmap, etc

#DimHeatmap(x, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(x, dims = 1:15, cells = 500, balanced = TRUE)


## UMAP before homony

#x <- x %>% FindNeighbors(reduction = "pca") %>% FindClusters(resolution = 0.5)# default resolution 0.8
#table(Idents(x))

x <- RunUMAP(x, dims = 1:13, verbose = F)

DimPlot(x, reduction = "umap", label = TRUE, pt.size = 0.01)# + NoLegend()
DimPlot(x, reduction = "umap", pt.size = 0.01)# + NoLegend()
#DimPlot(x, reduction = "umap", group.by = "orig.ident", pt.size = 0.01)


# run harmony

x <- x %>% RunHarmony("orig.ident")

#p1 <- DimPlot(object = x, reduction = "harmony", pt.size = .4, group.by = "orig.ident")
#p2 <- VlnPlot(object = x, features = "harmony_1", group.by = "orig.ident",  pt.size = .4)
#p1 + p2
#DimHeatmap(object = x, reduction = "harmony", cells = 500, dims = 1:3)

x <- x %>% RunUMAP(reduction = "harmony",  dims = 1:13)

saveRDS(x, file = "2-umap-noPKD2.rds", compress = TRUE)
x = readRDS("2-umap-noPKD2.rds")

x <- x %>% FindNeighbors(reduction = "harmony") %>% FindClusters(resolution=0.35)
table(Idents(x))

mycolor = c("#00bfc4", "#f8766d")
names(mycolor) = c("Cont","PKD")
#x$orig.ident = factor(x$orig.ident, levels =c("PKD", "Cont"))

p1 <- DimPlot(x,reduction="umap", label=TRUE, pt.size = 0.01)# + NoLegend()
p1
p2 <- DimPlot(x,reduction="umap", group.by="orig.ident", cols = mycolor, pt.size = 0.01)
p2
p1+p2

FeaturePlot(x, "PRKG1", pt.size = 0.01, raster=FALSE)


# label cell type
library(readxl)
type.marker = read_xlsx("~/Dropbox/result-NYU/scRNA-PKD/cell_type.xlsx", sheet="marker") %>% select(1:4)
type.marker = type.marker %>% pivot_longer(2:4, names_to = "x", values_to = "gene") %>% filter(!is.na(gene)) %>% select(-x)

cell.type = unique(type.marker$cell_type)
VlnPlot(x, raster=FALSE, pt.size=0.01, features=type.marker[which(type.marker$cell_type==cell.type[1]),]$gene)


#Assigning cell type identity to clusters

cellType = read_xlsx("~/Dropbox/result-NYU/scRNA-PKD/cell_type.xlsx", sheet="type_grp") %>% select(1:2) %>% mutate(grp = as.character(grp))

# method 1: replace idents (not good)

grp.cell <- cellType$cell_type
names(grp.cell) <- cellType$grp

y <- RenameIdents(x, grp.cell)

p1 <- DimPlot(y, reduction = "umap", label = TRUE, pt.size = 0.01)# + NoLegend()
p2 <- DimPlot(y, reduction = "umap", group.by = "orig.ident", pt.size = 0.01)
p1+p2

# method 2: add new slot

id_grp = Idents(x)
id.grp <- tibble(id = names(id_grp), grp = id_grp)

id.grp.cell = id.grp %>% left_join(cellType)

id.cell = id.grp.cell$cell_type
names(id.cell) = id.grp.cell$id

x[["cell.type"]] = id.cell
#table(x$cell.type)

p1 <- DimPlot(x, reduction="umap", label=TRUE, pt.size=0.01)
p1
p2 <- DimPlot(x, reduction="umap", label=TRUE, group.by="cell.type", pt.size=0.01) + NoLegend()
p2
p3 <- DimPlot(x, reduction="umap", group.by="orig.ident", pt.size=0.01, cols=mycolor)
p3
p1+p2+p3

saveRDS(x, file = "3-res0.35-celltype-noPKD2.rds", compress = TRUE)
x = readRDS("3-res0.35-celltype-noPKD2.rds")


# slim down object:

x.slim = x
x.slim[["RNA"]]$data.1 = NULL
x.slim[["RNA"]]$data.2 = NULL
x.slim[["RNA"]]$data.3 = NULL
x.slim[["RNA"]]$data.4 = NULL
x.slim[["RNA"]]$data.5 = NULL
x.slim[["RNA"]]$data.6 = NULL
x.slim[["RNA"]]$data.7 = NULL
x.slim[["RNA"]]$data.8 = NULL
x.slim[["RNA"]]$data.9 = NULL
x.slim[["RNA"]]$data.10 = NULL
x.slim[["RNA"]]$data.11 = NULL
x.slim[["RNA"]]$data.12 = NULL
x.slim[["RNA"]]$data.13 = NULL

x.slim[["RNA"]]$counts.1 = NULL
x.slim[["RNA"]]$counts.2 = NULL
x.slim[["RNA"]]$counts.3 = NULL
x.slim[["RNA"]]$counts.4 = NULL
x.slim[["RNA"]]$counts.5 = NULL
x.slim[["RNA"]]$counts.6 = NULL
x.slim[["RNA"]]$counts.7 = NULL
x.slim[["RNA"]]$counts.8 = NULL
x.slim[["RNA"]]$counts.9 = NULL
x.slim[["RNA"]]$counts.10 = NULL
x.slim[["RNA"]]$counts.11 = NULL
x.slim[["RNA"]]$counts.12 = NULL
x.slim[["RNA"]]$counts.13 = NULL

x.slim[["RNA"]]$scale.data.2 = NULL
x.slim[["RNA"]]$scale.data.3 = NULL
x.slim[["RNA"]]$scale.data.4 = NULL
x.slim[["RNA"]]$scale.data.5 = NULL
x.slim[["RNA"]]$scale.data.6 = NULL
x.slim[["RNA"]]$scale.data.7 = NULL
x.slim[["RNA"]]$scale.data.8 = NULL
x.slim[["RNA"]]$scale.data.9 = NULL
x.slim[["RNA"]]$scale.data.10 = NULL
x.slim[["RNA"]]$scale.data.11 = NULL
x.slim[["RNA"]]$scale.data.12 = NULL
x.slim[["RNA"]]$scale.data.13 = NULL

#x.slim@reductions$pca=NULL
#x.slim@reductions$harmony=NULL

p1 <- DimPlot(x.slim, reduction="umap", label=TRUE, group.by="cell.type", pt.size=0.01, raster=FALSE) + NoLegend()
p1
p2 <- DimPlot(x.slim, reduction="umap", group.by="orig.ident", pt.size=0.01, raster=FALSE)
p2
p1+p2

saveRDS(x.slim, file = "3-res0.35-celltype-noPKD2-slim.rds", compress = TRUE)
x.slim = readRDS("3-res0.35-celltype-noPKD2-slim.rds")

## random sample

ids <- colnames(x.slim)
#ids <- names(x.slim@active.ident)

sam <- ids |> map(\(x) str_remove(x, "_.+_")) |> unlist()
table(sam)
min = 200#min(table(pat)) #2492
sam.df <- tibble(id = ids, sam = sam)

id.rand <- sam.df |> split(sam) |> lapply(\(x) x |> pull(id) |> sample(min)) |> unlist()
x1 = x.slim[,id.rand]
table(x1$sample)
saveRDS(x1, file = "~/Dropbox/QiuDi/kidney1.rds", compress = TRUE)

ct = as.data.frame(x1@assays$integrated@scale.data)
ct = as.data.frame(round(t(ct),6))
ct = ct %>% mutate(id = rownames(ct))

cellType = tibble(id=colnames(x1), cell.type=x1$cell.type)

ct = ct %>% left_join(cellType)

n = ncol(ct)
y = ct %>% select(n-1, n, 1:(n-2))

write_tsv(y, "~/Dropbox/QiuDi/kidney-scRNA-for-class/kidney1.tsv")