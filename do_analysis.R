library(Seurat)
library(ggplot2)

setwd(# Instert Here Path For This Script #
  )
  
# LOAD DATASETS
pbmc_list <- list()
for (sample_n in c("PAT2V1",
                   "PAT2V2",
                   "PAT3V1",
                   "PAT3V3",
                   "PAT5V1",
                   "PAT5V2")){
  pbmc <- CreateSeuratObject(counts = Read10X(data.dir = paste0("data/", sample_n, "/outs/filtered_feature_bc_matrix/")),
                             project = sample_n)
  pbmc_list[sample_n] = pbmc
  rm(pbmc)
}

pbmc_list[[1]]$condition = "baseline"
pbmc_list[[3]]$condition = "baseline"
pbmc_list[[5]]$condition = "baseline"
pbmc_list[[2]]$condition = "CBD"
pbmc_list[[4]]$condition = "CBD"
pbmc_list[[6]]$condition = "CBD"

pbmc_list[[1]]$subject = "SUB1"
pbmc_list[[2]]$subject = "SUB1"
pbmc_list[[3]]$subject = "SUB2"
pbmc_list[[4]]$subject = "SUB2"
pbmc_list[[5]]$subject = "SUB3"
pbmc_list[[6]]$subject = "SUB3"

# MERGE DATASETS
pbmc_combined = merge(pbmc_list[[1]], y = c(pbmc_list[[2]], pbmc_list[[3]], pbmc_list[[4]], pbmc_list[[5]], pbmc_list[[6]]))
pbmc_combined[["percent.mt"]] <- PercentageFeatureSet(pbmc_combined, pattern = "^MT-")
VlnPlot(pbmc_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

pbmc_combined = subset(pbmc_combined, subset =
                         nFeature_RNA > 750 & nFeature_RNA < 5000 & percent.mt < 15) 
pbmc_combined = NormalizeData(pbmc_combined)
pbmc_combined = FindVariableFeatures(pbmc_combined, selection.method = "vst", nfeatures = 2000)
pbmc_combined = ScaleData(pbmc_combined, vars.to.regress = "percent.mt")

pbmc_combined = RunPCA(pbmc_combined, features = VariableFeatures(object = pbmc_combined), npcs = 100)
ElbowPlot(object = pbmc_combined, ndims = 100) # 25

pbmc_combined = RunUMAP(pbmc_combined, reduction = "pca", dims = 1:25)
DimPlot(pbmc_combined, reduction = "umap", split.by = "orig.ident")

# FIND CLUSTERS AND CLUSTER MARKERS
pbmc_combined <- FindNeighbors(pbmc_combined, dims = 1:25)
pbmc_combined <- FindClusters(pbmc_combined, resolution = 0.1)

for(cluster in levels(pbmc_combined$seurat_clusters)){
  cluster_cells = colnames(pbmc_combined)[pbmc_combined$seurat_clusters == cluster]
  other_cells = colnames(pbmc_combined)[pbmc_combined$seurat_clusters != cluster]
  markers = FindMarkers(pbmc_combined,
                        ident.1 = cluster_cells,
                        ident.2 = other_cells,
                        logfc.threshold = 0.5,
                        test.use = "poisson",
                        min.pct = 0.25)
  markers <- markers[order(-markers$avg_log2FC),]
  write.csv(markers, paste0("data/markers_cluster_", cluster, ".csv"), quote = FALSE)
}

# PER SUBJECT, CONDITION, CLUSTER VISUALIZATION
cols_subject = c("#E41A1C", "#377EB8", "#4DAF4A")
cols_condition = c("#E41A1C", "#377EB8")
cols_cluster = c("#8DD3C7", "#A48C8C", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD")


p = DimPlot(pbmc_combined, reduction = "umap", group.by = "subject", shuffle = TRUE, pt.size = 0.01, cols = cols_subject)
ggsave(plot = p, width = 20, height = 16, units = "cm", dpi = 600, filename = "figures/umap_subject.png")

p = DimPlot(pbmc_combined, reduction = "umap", group.by = "condition", shuffle = TRUE, pt.size = 0.01, cols = cols_condition)
ggsave(plot = p, width = 20, height = 16, units = "cm", dpi = 600, filename = "figures/umap_condition.png")

p = DimPlot(pbmc_combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, shuffle = TRUE, pt.size = 0.01, cols = cols_cluster)
ggsave(plot = p, width = 20, height = 16, units = "cm", dpi = 600, filename = "figures/umap_clusters.png")

# FROM CLUSTERS TO CELL POPULATIONS
# genes of interest for cell populations
genes_of_interest = c("IL7R",
                      "LEF1",
                      "CD28",
                      "CD8A",
                      "NKG7",
                      "KLRD1",
                      "KLRG1",
                      "LYZ",
                      "CD14",
                      "S100A8",
                      "S100A9",
                      "CD79A",
                      "CD79B",
                      "JCHAIN",
                      "IGHG1",
                      "IGHA1",
                      "IGHA2")
p = VlnPlot(object = pbmc_combined, features = genes_of_interest,
            ncol = 6, pt.size = 0, cols = cols_cluster)
ggsave(plot = p, width = 40, height = 25, units = "cm", dpi = 600, filename =
         paste0("figures/cell_clusters_genes_of_interest.png"))

# Idents based on cell populations
pbmc_combined$cell_population = "none"
pbmc_combined$cell_population[pbmc_combined$seurat_clusters == "0"] = "CD4 T cells"
pbmc_combined$cell_population[pbmc_combined$seurat_clusters %in% c("1", "2")] = "CD8 T cells"
pbmc_combined$cell_population[pbmc_combined$seurat_clusters %in% c("3", "8")] = "NKT cells"
pbmc_combined$cell_population[pbmc_combined$seurat_clusters == "4"] = "NK cells"
pbmc_combined$cell_population[pbmc_combined$seurat_clusters %in% c("5", "7")] = "Myeloid cells"
pbmc_combined$cell_population[pbmc_combined$seurat_clusters == "6"] = "B cells"
pbmc_combined$cell_population[pbmc_combined$seurat_clusters == "9"] = "Plasma\ncells"

cols_cell_pop = c("#E69F00","#009E73","#D55E00","brown1",
                  "#56B4E9","#F0E442","#0072B2")

p = DimPlot(pbmc_combined, reduction = "umap", group.by = "cell_population", label = TRUE, shuffle = TRUE, pt.size = 0.01,
            cols = cols_cell_pop)
ggsave(plot = p, width = 20, height = 16, units = "cm", dpi = 600, filename = "figures/umap_cell_populations.png")

Idents(pbmc_combined) = pbmc_combined$cell_population

for(type in c("NKT cells", "CD4 T cells", "NK cells", "Myeloid cells", "CD8 T cells", "B cells")){ # Plasma cells population too small (40 cells) to be considered
  baseline_cells = colnames(pbmc_combined)[pbmc_combined$cell_population == type & pbmc_combined$condition == "baseline"]
  cbd_cells = colnames(pbmc_combined)[pbmc_combined$cell_population == type & pbmc_combined$condition == "CBD"]
  markers = FindMarkers(pbmc_combined,
                        ident.1 = cbd_cells,
                        ident.2 = baseline_cells,
                        logfc.threshold = 0.5,
                        test.use = "poisson",
                        min.pct = 0.25)
  
  markers <- markers[order(markers$avg_log2FC),]
  write.csv(markers, paste0("data/markers_condition_", type, ".csv"), quote = FALSE)
}

# SCORE BASED ON GENE MODULE
pos_reg = readRDS("data/pos_reg.rds")
pos_reg_f = unique(pos_reg$external_gene_name[pos_reg$external_gene_name %in% row.names(pbmc_combined)])

pbmc_combined <- AddModuleScore(
  object = pbmc_combined,
  features = list(pos_reg_f),
  name = 'inflammation_features'
)

selected_cells = colnames(pbmc_combined)[pbmc_combined$cell_population=="Myeloid cells"]
inflammation_scores = pbmc_combined$inflammation_features1[selected_cells]
subjects = pbmc_combined$subject[selected_cells]
conditions = pbmc_combined$condition[selected_cells]

score_figure_df = data.frame(inflammation_score = inflammation_scores, subject = subjects, condition = conditions)
score_figure_df = score_figure_df[order(score_figure_df$inflammation_score),]

p = ggplot(score_figure_df, aes(x=condition, y=inflammation_score)) +
  geom_violin(trim=FALSE) + 
  geom_boxplot(width=0.1) +
  theme_bw() + facet_wrap(~subject)
ggsave(filename = paste0("figures/score_violin_plots.png"),
       width = 20, height = 10, units = "cm", dpi = 600,
       plot = p)
