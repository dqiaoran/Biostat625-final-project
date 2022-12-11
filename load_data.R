library(rhdf5)

## Load RNAseq data ##
h5ls("/home/asmauger/biostat625final/train_cite_inputs_raw.h5")
pathinput = "/home/asmauger/biostat625final/train_cite_inputs_raw.h5"
rna_data = h5read(pathinput, "/train_cite_inputs_raw/block0_values", index=list(1:22085, 1:70988))
rna_rows = h5read(pathinput, '/train_cite_inputs_raw/axis0')
rna_cols = h5read(pathinput, 'train_cite_inputs_raw/axis1')

## Load protein data ##
h5ls("/home/asmauger/biostat625final/train_cite_targets_raw.h5")
pathinput='/home/asmauger/biostat625final/train_cite_targets_raw.h5'
prot_data = h5read(pathinput, '/train_cite_targets_raw/block0_values')
prot_rows = h5read(pathinput, '/train_cite_targets_raw/axis0')
prot_cols = h5read(pathinput, '/train_cite_targets_raw/axis1')

all.equal(rna_cols, prot_cols)

## Load pre-normalized protein data ##

h5ls("/home/asmauger/biostat625final/train_cite_targets.h5")
pathinput='/home/asmauger/biostat625final/train_cite_targets.h5'
prot_data2 = h5read(pathinput, '/train_cite_targets/block0_values')
prot_rows2 = h5read(pathinput, '/train_cite_targets/axis0')
prot_cols2 = h5read(pathinput, '/train_cite_targets/axis1')

all.equal(prot_rows2, prot_rows)
all.equal(prot_cols2, prot_cols)

## Sparse matrix

library(Matrix)
smalldata = as(rna_data, Class='dgCMatrix')
dimnames(smalldata) = list(rna_rows, rna_cols)
smalldata_prot = as(prot_data, Class='dgCMatrix')
dimnames(smalldata_prot) = list(prot_rows, prot_cols)
smalldata_prot_norm = as(prot_data2, Class='dgCMatrix')
dimnames(smalldata_prot_norm) = list(prot_rows2, prot_cols2)

saveRDS(smalldata, file='sparse_RNA.rds')
saveRDS(smalldata_prot, file='sparse_prot.rds')
saveRDS(smalldata_prot_norm, file='sparse_prot_norm.rds')

## Use this to load sparse matrix:

sparse_RNA = readRDS('/home/asmauger/biostat625final/sparse_RNA.rds')
sparse_prot = readRDS('/home/asmauger/biostat625final/sparse_prot.rds')
sparse_prot_norm = readRDS('/home/asmauger/biostat625final/sparse_prot.rds')

all.equal(dim(sparse_RNA), c(22085, 70988))
all.equal(dim(sparse_prot), c(140, 70988))
all.equal(dim(sparse_prot_norm), c(140, 70988))

#### SPLITTING DATA ####

library(dplyr)
set.seed(347382)

## Load row names ##
pathinput = "/home/asmauger/biostat625final/train_cite_inputs_raw.h5"
rna_cols = h5read(pathinput, 'train_cite_inputs_raw/axis1')

## Load metadata
metadata = read.csv('/home/asmauger/biostat625final/metadata.csv', header=T)
# filter metadata to cell ids we have data for and add a 'row number' column
grouped_metadata = metadata %>% filter(cell_id %in% rna_cols) %>% mutate(id=1:70988) %>% group_by(donor, day)

# make sure the metadata retained the order of the data (rna_cols)
test = grouped_metadata %>% ungroup() %>% select(cell_id)
all.equal(as.matrix(test$cell_id), as.matrix(rna_cols))

# 70988 samples, 9 groups (3 donors, 3 days)
ngroup = round(70988/9)
nsplit = round(70988/18)
mysample = grouped_metadata %>% slice_sample(n=nsplit)
count(mysample) # check that the sampling worked evenly
indexes = sample %>% ungroup() %>% select(id)
length(indexes$id)

### RNA training split

pathinput = "/home/asmauger/biostat625final/train_cite_inputs_raw.h5"
rna_data_train = h5read(pathinput, "/train_cite_inputs_raw/block0_values", index=list(1:22085, indexes$id))
rna_rows_train = h5read(pathinput, '/train_cite_inputs_raw/axis0')
rna_cols_train = h5read(pathinput, '/train_cite_inputs_raw/axis1', index=list(indexes$id))

h5ls(pathinput)


h5createFile("rna_train.h5")
h5createGroup("rna_train.h5","")
h5createDataset('rna_train.h5', 'rawcounts', c(22085, 35496))
h5writeDataset(obj = rna_data_train, h5loc = "rna_train.h5", name = 'rawcounts')
h5createDataset('rna_train.h5', 'genes', 22085, chunk = 1000)
h5writeDataset(as.matrix(rna_rows_train), 'rna_train.h5', 'genes')
h5createDataset('rna_train.h5', 'cells', length(indexes$id), chunk = 1000)
h5writeDataset(as.matrix(rna_cols_train, 'rna_train.h5', 'cells')

## Seurat object

library(Seurat)
metadata = read.csv('/home/asmauger/biostat625final/metadata.csv', header=T)
seuratobj = CreateSeuratObject(counts = sparse_RNA, meta.data=metadata)
prot_assay = CreateAssayObject(counts=sparse_prot)
seuratobj[['prot']] = prot_assay

Assays(seuratobj) # should have RNA and prot

DefaultAssay(mydata2) = 'prot' # set default assay like so


