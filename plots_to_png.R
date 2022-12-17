## Finalizing plots to save as png

library(ggplot2)
overall_umap <- readRDS('/home/asmauger/biostat625final/overall_umap.rds')
overall_umap +
  scale_color_discrete(labels=c('Testing data', 'Training data')) +
  labs(title='sciPENN embedding UMAP', x='UMAP1', y='UMAP2') +
  theme(text=element_text(size=12))

ggsave('/home/asmauger/biostat625final/Biostat625-final-project/overall_umap.png', width=6, height=4)

