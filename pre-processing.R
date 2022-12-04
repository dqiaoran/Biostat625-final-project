# load data
sparse_RNA = readRDS('/home/asmauger/biostat625final/sparse_RNA.rds')
sparse_prot = readRDS('/home/asmauger/biostat625final/sparse_prot.rds')
sparse_prot_norm = readRDS('/home/asmauger/biostat625final/sparse_prot.rds')

all.equal(dim(sparse_RNA), c(22085, 70988))
all.equal(dim(sparse_prot), c(140, 70988))
all.equal(dim(sparse_prot_norm), c(140, 70988))

# extract gene symbol from rownames
# rownames look like this : "ENSG00000081913-PHLPP1" (ensembl id - gene symbol)
library(stringr)
names = rownames(sparse_RNA)
hyphen_symbols = str_extract(names, "_.*")
gene_sym = str_sub(hyphen_symbols, 2)

# change rownames to just gene symbols
rownames(sparse_RNA) = gene_sym


# create seurat object
library(Seurat)
metadata = read.csv('/home/asmauger/biostat625final/metadata.csv', header=T)
rownames(metadata) = metadata[,1]
# add RNA and metadata
seuratobj = CreateSeuratObject(counts = sparse_RNA, meta.data=metadata)
# add raw protein
prot_assay = CreateAssayObject(counts=sparse_prot)
seuratobj[['prot']] = prot_assay
# add normalized protein
prot_norm_assay = CreateAssayObject(data = sparse_prot_norm)
seuratobj[['protnorm']] = prot_norm_assay

# remove original data
rm(sparse_RNA) # access with seuratobj@assays$RNA@counts
rm(sparse_prot) # access with seuratobj@assays$prot@counts
rm(sparse_prot_norm) # access with seuratobj@assays$protnorm@counts

## QC ??? not sure what we want to do here
VlnPlot(seuratobj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
# nFeature_RNA is number of unique genes expressed in a given cell
# nCount_RNA is number of transcripts in a given cell
rna_counts = seuratobj@assays$RNA@counts
all.equal(colSums(rna_counts), seuratobj$nCount_RNA)
rowSums(rna_counts!=0) # number of cells for which a gene has nonzero expression

gene_sym[grep('MT', gene_sym, perl=T)]


## Normalization
seuratobj = NormalizeData(seuratobj, normalization.method = "LogNormalize", scale.factor = 10000)

DefaultAssay(seuratobj) = 'prot'
# not sure if this is best way to normalize protein
seuratobj = NormalizeData(seuratobj, normalization.method = 'CLR', margin = 2)

prot = seuratobj@assays$prot@counts
#install.packages('dsb') - this is done in 'protnorm'
#library(dsb) # https://github.com/niaid/dsb
#prot_norm = DSBNormalizeProtein( # we need the empty droplets to do this normalization
#  cell_protein_matrix = prot,
#  empty_drop_matrix = NULL,
#  denoise.counts = TRUE,
#  use.isotype.control = TRUE,
#  isotype.control.name.vec = rownames(prot)[c(32:35, 85:86)]
#)
#32,33,34,35, 85, 86 are mouse/rat antibodies
rownames(prot)[c(32:35, 85:86)]


## Variable features
DefaultAssay(seuratobj) = 'RNA'
seuratobj = FindVariableFeatures(seuratobj)

# here they do not include isotype controls in analysis
# https://cran.r-project.org/web/packages/dsb/vignettes/end_to_end_workflow.html
# set 'variable features' of proteins to be the non-isotype control proteins
DefaultAssay(seuratobj) = 'prot'
VariableFeatures(seuratobj) = rownames(prot)[-c(32:35, 85:86)]

DefaultAssay(seuratobj) = 'protnorm'
VariableFeatures(seuratobj) = rownames(prot)[-c(32:35, 85:86)]

# Forward from this point, only variable features will be used by default

## Scale data
seuratobj =ScaleData(seuratobj, assay='RNA')
seuratobj = ScaleData(seuratobj, assay='prot')
seuratobj = ScaleData(seuratobj, assay='protnorm')

## PCA
seuratobj = RunPCA(seuratobj, assay='RNA')
seuratobj = RunPCA(seuratobj, assay='prot', reduction.name='protpca')
seuratobj = RunPCA(seuratobj, assay='protnorm', reduction.name='protnormpca')

## UMAP
ElbowPlot(seuratobj)
ElbowPlot(seuratobj, reduction = 'protpca')
ElbowPlot(seuratobj, reduction = 'protnormpca')
seuratobj = RunUMAP(seuratobj, dims=1:20)
seuratobj = RunUMAP(seuratobj, dims=1:10, reduction='protpca', reduction.name='protumap')
seuratobj = RunUMAP(seuratobj, dims=1:10, reduction='protnormpca', reduction.name='protnormumap')

p1 = DimPlot(seuratobj, reduction='umap', group.by='day', shuffle = T)
p2 = DimPlot(seuratobj, reduction='umap', group.by='cell_type', shuffle=T)
p3 = DimPlot(seuratobj, reduction='protumap', group.by='day', shuffle=T)
p4 = DimPlot(seuratobj, reduction='protumap', group.by='cell_type', shuffle=T)
p5 = DimPlot(seuratobj, reduction='protnormumap', group.by='day', shuffle=T)
p6 = DimPlot(seuratobj, reduction='protnormumap', group.by='cell_type', shuffle=T)
library(patchwork)
p1 + p2 # RNA
p3 + p4 # protein - CLR normalized
p5 + p6 # protein - dsb normalized

# weighted nearest neighbors

# RNA & CLR-normalized protein
seuratobj = FindMultiModalNeighbors(
  seuratobj, reduction.list = list("pca", "protpca"),
  dims.list = list(1:30, 1:18), # default values, not sure what is best here
  modality.weight.name = "RNA.weight"
)
seuratobj <- RunUMAP(seuratobj, nn.name = "weighted.nn", reduction.name = "wnn.umap")
p7 = DimPlot(seuratobj, reduction='wnn.umap', group.by='cell_type', shuffle=T)
p8 = DimPlot(seuratobj, reduction='wnn.umap', group.by='day', shuffle=T)
p7+p8

# RNA & dsb-normalized protein
seuratobj = FindMultiModalNeighbors(
  seuratobj, reduction.list = list("pca", "protnormpca"),
  weighted.nn.name = "weighted.nn.dsb",
  dims.list = list(1:30, 1:18), # default values, not sure what is best here
  modality.weight.name = "RNA.weight2"
)
seuratobj <- RunUMAP(seuratobj, nn.name = "weighted.nn.dsb", reduction.name = "wnn.dsb.umap")
p9 = DimPlot(seuratobj, reduction='wnn.dsb.umap', group.by='cell_type', shuffle=T)
p10 = DimPlot(seuratobj, reduction='wnn.dsb.umap', group.by='day', shuffle=T)
p9 + p10

# Other summaries

library(dplyr)
metadata %>% filter(cell_id %in% colnames(seuratobj)) %>% select(day, donor, cell_type) %>% group_by(day, donor) %>% summarise(n())
metadata %>% filter(cell_id %in% colnames(seuratobj)) %>% select(day, donor, cell_type) %>% group_by(cell_type) %>% summarise(n())
