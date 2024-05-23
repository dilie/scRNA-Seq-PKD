library(tidyverse)
library(Seurat)
#library(harmony)
#library(cowplot)

library(pheatmap)
library(gplots)

setwd("~/Desktop/local-scRNA-PKD/")

x = readRDS("3-res0.35-celltype-noPKD2-slim.rds")
#x = readRDS("3-res0.35-celltype-noPKD2.rds")

#Method 1: t-test for sub-group
library(broom)

subGrp = function(grp,res){
  y = subset(x, idents = grp)
  y %>% FindNeighbors(reduction="harmony") %>% FindClusters(resolution=res)
}

obj2tb = function(y){
  ct = as.data.frame(y@assays$integrated@scale.data)
  ct = ct %>% mutate(gene=rownames(ct))
  ct.long = ct %>% pivot_longer(1:(ncol(ct)-1), names_to = "id", values_to = "count")
  
  tb = tibble(id=colnames(y), pkd=y$orig.ident, subGrp=Idents(y))
  ct.long %>% left_join(tb) %>% mutate(pkd = ifelse(pkd=="Cont","Ctrl",pkd)) %>% mutate(sam = str_remove(id, "_.+_"))
}


##########
thisGrp = "11"
y = subGrp(thisGrp, 0.1)
DimPlot(y, reduction="umap", label=TRUE, pt.size=0.01, raster=FALSE) + NoLegend()
table(Idents(y))

#mycolor = c("#00bfc4", "#f8766d")
#names(mycolor) = c("Cont","PKD")
#DimPlot(y, reduction="umap", group.by="orig.ident", pt.size=0.01, raster=FALSE, cols = mycolor)

#if not sub-group
#y = subset(x, idents = thisGrp)

ct = obj2tb(y)

#1) t-test by sample: better

# option 1: combine subGrps or no subGrp
ct.med = ct %>% group_by(gene,pkd,sam) %>% summarise(count=median(count))
t.out = ct.med %>% group_by(gene) %>% do(tidy(t.test(data=., count~pkd)))
t.sig = t.out %>% filter((estimate > 1 | estimate < -1) & p.value<1e-03) %>% select(gene, estimate, p.value) %>% rename(fc = estimate, pv = p.value)
ct.sig = ct.med %>% left_join(t.sig, by=c("gene")) %>%
  filter(!is.na(fc)) %>%
  select(-fc,-pv) %>% 
  mutate(grp=thisGrp, subGrp="")

# option 2: sep subGrps
ct.med = ct %>% group_by(subGrp,gene,pkd,sam) %>% summarise(count=median(count))
t.out = ct.med %>% group_by(subGrp,gene) %>% do(tidy(t.test(data=., count~pkd)))
t.sig = t.out %>% filter((estimate > 1 | estimate < -1) & p.value<1e-03) %>%
  select(gene, subGrp, estimate, p.value) %>%
  rename(fc = estimate, pv = p.value)
ct.sig = ct.med %>% left_join(t.sig, by=c("gene","subGrp")) %>%
  filter(!is.na(fc)) %>%
  select(-fc,-pv) %>% 
  mutate(grp=thisGrp)

ct.all = ct.all %>% bind_rows(ct.sig)
#ct.all = ct.sig


########
# install and plot enhanced volcana
#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

t0 <- t.out %>% filter(subGrp == '0')
EnhancedVolcano(t0, lab = t0$gene, x = "estimate", y = "p.value", title = cell, subtitle = paste(thisGrp, "-0", sep=""), pCutoff = 1e-04)

t1 <- t.out %>% filter(subGrp == '1')
EnhancedVolcano(t1, lab = t1$gene, x = "estimate", y = "p.value", title = cell, subtitle = paste(thisGrp, "-1", sep=""), pCutoff = 1e-04)

#if no sub-group
EnhancedVolcano(t.med, lab=t.med$gene, x="estimate", y="p.value", title = cell, subtitle=thisGrp)

#2) t-test by cell
#t.out = ct %>% group_by(grp, gene) %>% do(tidy(t.test(data = ., count ~ pkd)))

#t.out0 <- t.out %>% filter(grp == '0')
#EnhancedVolcano(t.out0, lab = t.out1$gene, x = "estimate", y = "p.value")

#t.out1 <- t.out %>% filter(grp == '1')
#EnhancedVolcano(t.out1, lab = t.out1$gene, x = "estimate", y = "p.value")

##boxplot

# if subGrp:
subG = "0"
thisGene = "WDFY4"
thisTitle = paste(cell, "-", subG, ": ", thisGene, sep = "")

z = ct %>% filter(gene==thisGene & subGrp==subG) %>% mutate(sam = str_remove(id, "_.+_"))

#if no subGrp:
thisTitle = paste(cell, ": ", thisGene, sep = "")
z = ct %>% filter(gene==thisGene) %>% mutate(sam = str_remove(id, "_.+_"))


z %>%
  ggplot(aes(x=sam, y=count, color=pkd)) +
  geom_boxplot(outlier.colour = NA, color="gray") +
  geom_jitter(position=position_jitter(0.2), shape=1) +
  theme_light() +
  scale_color_manual(values=c("#00bfc4", "#f8766d")) +
  labs(title = thisTitle)

#ct.all = z %>% mutate(tt=thisTitle)

z = z %>% mutate(tt=thisTitle)
ct.all = bind_rows(ct.all, z)
write_tsv(ct.all, "~/Dropbox/result-NYU/scRNA-PKD/count-sig.tsv")
ct.all = read_tsv("~/Dropbox/result-NYU/scRNA-PKD/count-sig.tsv")

# heatmap
library(readxl)
cellType = read_xlsx("~/Dropbox/result-NYU/scRNA-PKD/cell_type.xlsx", sheet="type_grp") %>% select(1:2) %>% mutate(grp = as.character(grp))

z = ct.all %>% left_join(cellType) %>% mutate(tt = paste(cell_type, ":", subGrp, ":", gene))
z = tibble(tt=z$tt, pkd=z$pkd, count=z$count, sam=z$sam)
z.wide = z %>% pivot_wider(names_from = "tt", values_from = "count")
z.df = as.data.frame(z.wide)[, -1]
rownames(z.df) = z.df$sam
z.df = z.df[,-1]
z.df = as.data.frame(t(as.matrix(z.df)))

ann.col = data.frame(samType=z.wide$pkd, row.names = z.wide$sam)

row = rownames(z.df)
ann.row = data.frame(cellType=str_remove(row, " :.+"), row.names = row)

pheatmap(z.df, color = colorpanel(500,"white","blue"), annotation_col = ann.col, annotation_row = ann.row, scale = "row")

#####

ct.stat = ct.all %>% group_by(grp, subGrp, gene) %>% summarise(number=length(sam))
tmp = ct.stat %>% group_by(gene) %>% summarise(nGene = length(grp))

tmp2 = ct.stat %>% select(-number) %>% left_join(cellType) %>% arrange(cell_type, subGrp, gene)
write_tsv(tmp2[,c(4,2,3)], "~/Downloads/sig-genes.tsv")

#########

ct.all %>%
  ggplot(aes(x = sam, y = count, color=pkd)) +
  geom_boxplot(outlier.colour = NA, color="gray") +
  geom_jitter(position=position_jitter(0.2), shape=".") +
  facet_wrap(~tt, ncol = 2) +
  theme_light() +
  scale_color_manual(values=c("#00bfc4", "#f8766d")) +
  theme(text = element_text(size = 14),
        axis.text = element_text(color = "#777777"),
        legend.title=element_blank(),
        legend.key.size = unit(0.75,"cm")) +
  xlab("") + ylab("")

#########

# sub-group DCT (3)
y %>% subset(idents = 1) %>% VlnPlot(features = "SLC12A3", group.by = "orig.ident", cols = mycolor)
y %>% subset(idents = 0) %>% VlnPlot(features = "TRPM6", group.by = "orig.ident", cols = mycolor)

#y %>% subset(idents = 1) %>% VlnPlot(features = "TMEM52B", group.by = "orig.ident", cols=mycolor)
ct.all = ct.all %>% filter(tt!="DCT-1: TMEM52B")


# sub-group FIB (5)
y %>% subset(idents = 1) %>% VlnPlot(features = "EBF2", group.by = "orig.ident", cols = mycolor)
#y %>% subset(idents = 1) %>% VlnPlot(features = "SEMA4A", group.by = "orig.ident", cols = mycolor)
y %>% subset(idents = 1) %>% VlnPlot(features = "HPSE2", group.by = "orig.ident", cols = mycolor)

y %>% subset(idents = 0) %>% VlnPlot(features = "DGKG", group.by = "orig.ident", cols = mycolor)

# sub-group LEUK (8)
y %>% subset(idents = 0) %>% VlnPlot(features = "ARHGAP15", group.by = "orig.ident", cols = mycolor)
y %>% subset(idents = 1) %>% VlnPlot(features = "WDFY4", group.by = "orig.ident", cols = mycolor)
#y %>% subset(idents = 0) %>% VlnPlot(features = "PRKCB", group.by = "orig.ident", cols = mycolor)
#y %>% subset(idents = 1) %>% VlnPlot(features = "ITGAX", group.by = "orig.ident", cols = mycolor)

# sub-group CNT_PC (2)
y %>% subset(idents = 0) %>% VlnPlot(features = "LINC01098", group.by = "orig.ident", cols = mycolor)
y %>% subset(idents = 0) %>% VlnPlot(features = "LINC01099", group.by = "orig.ident", cols = mycolor)


# sub-group PT2 (1)
#y %>% subset(idents = 1) %>% VlnPlot(features = "AC079298.3", group.by = "orig.ident", cols = mycolor)

# sub-group PT1 (4)

# sub-group TAL1 (6)
y <- y %>% FindNeighbors(reduction = "harmony") %>% FindClusters(resolution=0.15)

# sub-group TAL2 (0)

# sub-group ENDO (7)

# sub-group ICA-ICB (9)
y %>% subset(idents = 0) %>% VlnPlot(features = "SLC26A4", group.by = "orig.ident", cols = mycolor)

# sub-group PEC (10)
#no sub groups

# sub-group PODO (11)

#######


#Method 2: heatmap for each group
out = SplitObject(x, split.by ="ident") %>% lapply(function(obj) {
  tmp = FindMarkers(obj, ident.1="PKD", group.by = "orig.ident") %>% slice_head(n=30)
  tmp$gene = rownames(tmp)
  rownames(tmp) = NULL
  tmp$cell.type = unique(obj$cell.type)
  tmp$grp = unique(Idents(obj))
  tmp
})
markers = bind_rows(out)
write_tsv(markers, "~/Dropbox/result-NYU/scRNA-PKD/sigGene-top30.tsv")

heatmap.singleGrp = function(y.sig){
  y.mat = y.sig@assays$integrated@scale.data
  y.df = as.data.frame(y.mat)
  df <- data.frame(pkd = y.sig$orig.ident, row.names = colnames(y.sig))
  pheatmap(y.df, show_colnames=F, annotation_col = df, color = colorpanel(1000,"white","blue"), main = unique(y.sig$cell.type), annotation_legend = F)
}

p = list()
for (id in 0:12) {
  g = as.character(id)
  genes = markers[which(markers$grp==g),]$gene
  x.sig = subset(x, features = genes, idents = g)
  p[[id+1]] = heatmap.singleGrp(x.sig)[[4]]
}

library(gridExtra)
grid.arrange(p[[2]], p[[5]], p[[11]], nrow=1)
grid.arrange(p[[6]], p[[1]], p[[4]], nrow=1)
grid.arrange(p[[3]], p[[9]], p[[12]], nrow=1)
grid.arrange(p[[10]], p[[7]], p[[8]], nrow=1)
grid.arrange(p[[13]], nrow=1)





#x %>% VlnPlot(features = marker[1:15], group.by = "orig.ident")

x.sig = subset(x, features = marker[1:30], idents = as.character(0:13))
myheatmap(x.sig)

myheatmap = function(x.sig){
    id.grp = Idents(x.sig)
    id.rand <- id.grp |> split(id.grp) |> lapply(\(x) names(x) |> sample(100)) |> unlist()
    
    y.sig = x.sig[,id.rand]
    
    y.mat = y.sig@assays$integrated@scale.data
    y.df = as.data.frame(y.mat)
    
    tb = tibble(id=colnames(y.sig), pkd=y.sig$orig.ident, cell=y.sig$cell.type)
    cells <- sort(unique(tb$cell))
    
    for(i in 1:length(cells)) {
      c.id <- cells[i]
      tb <- tb %>% mutate(mu = if_else(cell == c.id, 1, 0))
      colnames(tb)[i+3] = c.id
    }
    df = as.data.frame(tb) %>% select(-cell)
    rownames(df) = df$id
    df = df[,-1]
    
    pheatmap(y.df, show_colnames=F, annotation_col = df, color = colorpanel(1000,"white","blue"), annotation_legend = F)
}

## single cell type

sub.heatmap = function(cell){
    xx = subset(x, cell.type==cell)
    table(xx$orig.ident)
    xx.gene = FindMarkers(xx, ident.1 = "PKD", group.by = "orig.ident", only.pos = TRUE)
#    xx %>% VlnPlot(features = rownames(xx.gene)[1:6], group.by = "orig.ident")
    myheatmap(subset(xx, features = rownames(dct.gene)[1:30]))
}

sub.heatmap("TAL1")

#############################


DimHeatmap(wt, dims = 1, cells = 500, reduction = "harmony")

#markers0 <- FindMarkers(x, ident.1=0, logfc.threshold=0.25, test.use="roc", only.pos=TRUE)
#markers1 <- FindMarkers(x, ident.1=1, logfc.threshold=0.25, test.use="roc", only.pos=TRUE)

markers.wt <- FindAllMarkers(wt, only.pos = TRUE)

top10 <- markers.wt %>%
  group_by(cluster) %>% 
  dplyr::filter(avg_log2FC > 1) %>% 
  slice_head(n = 10) %>% 
  ungroup()

write.table(top10, "top10-wt.csv", row.names = FALSE, sep = ",", quote = FALSE)
#marker = read_csv("top10-wt.csv", col_names = FALSE)
#VlnPlot(x, features = marker$X1)

wt.top = subset(wt, features = unique(top10$gene))

mat = wt.top@assays$integrated@scale.data

df = as.data.frame(mat)

id_grp = Idents(wt.top)
id.df <- data.frame(grp = id_grp, row.names = names(id_grp))

grp.col <- rainbow(21) # character vector of 22 rainbow colors
names(grp.col) <- levels(id_grp) # add names to the color vector (hash)

ann_color <- list(grp = grp.col)
pheatmap(df, labels_col = '', annotation_col = id.df, annotation_colors = ann_color)




x.sig = subset(x, features = unique(type.marker$gene), idents = as.character(0:8))

id.grp = Idents(x.sig)
id.rand <- id.grp |> split(id.grp) |> lapply(\(x) names(x) |> sample(50)) |> unlist()

y.sig = x.sig[,id.rand]

y.mat = y.sig@assays$integrated@scale.data

y.df = as.data.frame(y.mat)

#id.grp = Idents(y.sig)
id.df <- data.frame(grp = id.grp, row.names = names(id.grp))

grp.col <- rainbow(9) # character vector of 22 rainbow colors
names(grp.col) <- levels(id.grp) # add names to the color vector (hash)

ann_color <- list(grp = grp.col)
pheatmap(y.df, labels_col = '', annotation_col = id.df, annotation_colors = ann_color)




sig.gene = c("SLC34A1", "SLC5A2","CDH6","VCAM1","CCL2","PRICKLE1","TGFB2","DCC","CRYAB","CFH","ALDH1A2","UMOD","CLC12A1","CLDN16","PAPPA2","TMPRSS4","COL7A1","SLC12A3","TRPM6","KCNH7","SLC8A1-AS1","CALB1","AQP2","SLC14A2","SLC26A7","CFTR","SLC26A4","DGKI","NPHS2","PLA2R1","EMCN","FLT1","PDGFRA","PDGFRB","COL1A1","FBLN1","PTPRC","SLC11A1","UPK3A","PSCA")
#sig.gene = read.table("sig-gene.txt")

sig.gene[sig.gene %in% top10$gene]

# You can use the LabelClusters function to help label individual clusters


markers.all <- FindAllMarkers(x, only.pos = TRUE)
top10 <- markers.all %>%
  group_by(cluster) %>% 
  dplyr::filter(avg_log2FC > 1) %>% 
  slice_head(n = 10) %>% 
  ungroup()

DimHeatmap(x, dims = 1, cells = 500)