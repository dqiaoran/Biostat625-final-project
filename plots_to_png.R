## Finalizing plots to save as png

library(ggplot2)
library(patchwork)
library(forcats)

## Exploratory data analysis
cell_proportions <- readRDS('/home/asmauger/biostat625final/cell_proportions.rds')
cell_type_umap <- readRDS('/home/asmauger/biostat625final/cell_type_umap.rds')

cell_proportions_edited <- cell_proportions +
  theme(text=element_text(size=12), legend.text = element_text(size=10),legend.key.size = unit(.2, 'cm')) +
  labs(y='avg proportion of cells \n assigned to each cell type')
# brewer.pal(7, 'Set2') # figure out color names to match plots
cell_type_umap_edited <- cell_type_umap +
  theme(text=element_text(size=12), legend.position='none') +
  scale_color_discrete(type=c('#E5C494', "#FC8D62", "#66C2A5", "#E78AC3", "#A6D854", "#FFD92F", "#8DA0CB" )) +
  labs(title='')

cell_type_umap_edited + cell_proportions_edited

ggsave('/home/asmauger/biostat625final/Biostat625-final-project/exploratory.png', width=10, height=4)

## Overall umap - sciPENN embedding
overall_umap <- readRDS('/home/asmauger/biostat625final/overall_umap.rds')
overall_umap +
  scale_color_discrete(labels=c('Testing data', 'Training data')) +
  labs(title='sciPENN embedding UMAP', x='UMAP1', y='UMAP2') +
  theme(text=element_text(size=12))

ggsave('/home/asmauger/biostat625final/Biostat625-final-project/overall_umap.png', width=6, height=4)

## protein UMAPs
protein_predicted_umap <- readRDS('/home/asmauger/biostat625final/protein_predicted_umap.rds')
protein_true_umap <- readRDS('/home/asmauger/biostat625final/protein_true_umap.rds')

protein_predicted_umap_edit <- protein_predicted_umap +
  plot_annotation(title = 'Predicted relative protein expression', theme = theme(plot.title = element_text(size = 16))) &
  theme(text=element_text(size=12))
protein_true_umap_edit <- protein_true_umap +
  plot_annotation(title = 'True relative protein expression', theme = theme(plot.title = element_text(size = 16))) &
  theme(text=element_text(size=12))
wrap_elements(protein_predicted_umap_edit) + wrap_elements(protein_true_umap_edit)

ggsave('/home/asmauger/biostat625final/Biostat625-final-project/expression_umap.png', width=10, height=4)
