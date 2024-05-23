#data source:
#1. Defining cellular complexity in human autosomal dominant polycystic kidney disease by multimodal single cell analysis. Nat Commun. 2022 (https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/36310237/)
#2. Single-cell RNA sequencing of human kidney. scientific data. 2020 (https://www.nature.com/articles/s41597-019-0351-8#Sec11)

library(tidyverse)
library(Seurat)
setwd("~/Desktop/local-scRNA-PKD/")

x = readRDS("0-GSE185948_count_RNA.rds")
x <- CreateSeuratObject(counts=x, project = "kidney", min.cells = 3, min.features = 200)

head(x@meta.data)
table(x@meta.data$orig.ident)

# subset random sample from each sample
ids <- names(x.full@active.ident)
sam <- ids |> map(\(x) str_remove(x, "_.+_")) |> unlist()
table(sam)
min = min(table(sam)) #2492
sam.df <- tibble(id = ids, sam = sam)
id.rand <- sam.df |> split(sam) |> lapply(\(x) x |> pull(id) |> sample(min)) |> unlist()

x = x.full[,id.rand]
table(x$orig.ident)
rm(x.full, sam, id.rand, ids, sam.df, min)
####

x[["sample"]] = str_remove(colnames(x), "_.+_")
table(x$sample)

#x <- subset(x, idents = "Cont")
#x <- subset(x, idents = "PKD")


# Integration

x.list <- SplitObject(x, split.by = "sample")
rm(x)

for (i in 1:length(x.list)) {
  x.list[[i]] <- NormalizeData(x.list[[i]], verbose = FALSE)
  x.list[[i]] <- FindVariableFeatures(x.list[[i]], selection.method="vst",nfeatures=2000,verbose =FALSE)
}
rm(i)

anchor <- FindIntegrationAnchors(object.list = x.list)
rm(x.list)

x <- IntegrateData(anchorset = anchor)
rm(anchor)

saveRDS(x, file = "1-integrated.rds", compress = TRUE)
