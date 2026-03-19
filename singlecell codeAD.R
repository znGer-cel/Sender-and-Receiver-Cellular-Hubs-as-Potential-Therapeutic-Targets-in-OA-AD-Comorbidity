########Data Processing, Filtering, and Integration#########
#load packages#
library(Seurat)
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(stringr)

###load_data###
#load_data+creat_seurat#
rm(list = ls())
AD1 <- Read10X(data.dir = "GSE157827/AD1/") 
AD1 <- CreateSeuratObject(counts = AD1, project = "AD1", min.cells = 3, min.features = 200)
AD2 <- Read10X(data.dir = "GSE157827/AD2/") 
AD2 <- CreateSeuratObject(counts = AD2, project = "AD2", min.cells = 3, min.features = 200)
AD4 <- Read10X(data.dir = "GSE157827/AD4/") 
AD4 <- CreateSeuratObject(counts = AD4, project = "AD4", min.cells = 3, min.features = 200)
AD5 <- Read10X(data.dir = "GSE157827/AD5/") 
AD5 <- CreateSeuratObject(counts = AD5, project = "AD5", min.cells = 3, min.features = 200)
AD9 <- Read10X(data.dir = "GSE157827/AD9/") 
AD9 <- CreateSeuratObject(counts = AD9, project = "AD9", min.cells = 3, min.features = 200)
AD10 <- Read10X(data.dir = "GSE157827/AD10/") 
AD10 <- CreateSeuratObject(counts = AD10, project = "AD10", min.cells = 3, min.features = 200)
AD20 <- Read10X(data.dir = "GSE157827/AD20/") 
AD20 <- CreateSeuratObject(counts = AD20, project = "AD20", min.cells = 3, min.features = 200)
AD21 <- Read10X(data.dir = "GSE157827/AD21/") 
AD21 <- CreateSeuratObject(counts = AD21, project = "AD21", min.cells = 3, min.features = 200)
NC3 <- Read10X(data.dir = "GSE157827/NC3/") 
NC3 <- CreateSeuratObject(counts = NC3, project = "NC3", min.cells = 3, min.features = 200)
NC7 <- Read10X(data.dir = "GSE157827/NC7/") 
NC7 <- CreateSeuratObject(counts = NC7, project = "NC7", min.cells = 3, min.features = 200)
NC12 <- Read10X(data.dir = "GSE157827/NC12/") 
NC12 <- CreateSeuratObject(counts = NC12, project = "NC12", min.cells = 3, min.features = 200)
NC14 <- Read10X(data.dir = "GSE157827/NC14/") 
NC14 <- CreateSeuratObject(counts = NC14, project = "NC14", min.cells = 3, min.features = 200)
NC15 <- Read10X(data.dir = "GSE157827/NC15/") 
NC15 <- CreateSeuratObject(counts = NC15, project = "NC15", min.cells = 3, min.features = 200)
NC17 <- Read10X(data.dir = "GSE157827/NC17/") 
NC17 <- CreateSeuratObject(counts = NC17, project = "NC17", min.cells = 3, min.features = 200)
NC18 <- Read10X(data.dir = "GSE157827/NC18/") 
NC18 <- CreateSeuratObject(counts = NC18, project = "NC18", min.cells = 3, min.features = 200)
#merge_data#
merged_seurat <- merge(AD1, y = c(AD2,AD4,AD5,AD9,AD10,AD20,AD21,
                                  NC3,NC7,NC12,NC14,NC15,NC17,NC18),
                 add.cell.ids = c("AD1","AD2","AD4","AD5","AD9","AD10","AD20","AD21",
                                 "NC3","NC7","NC12","NC14","NC15","NC17","NC18"))
#count_UMI_MT#
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
#creat_metadata_include_base&sample&group_information#
metadata <- merged_seurat@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "AD1"))] <- 'AD1'
metadata$sample[which(str_detect(metadata$cells, "AD2"))] <- 'AD2'
metadata$sample[which(str_detect(metadata$cells, "AD4"))] <- 'AD4'
metadata$sample[which(str_detect(metadata$cells, "AD5"))] <- 'AD5'
metadata$sample[which(str_detect(metadata$cells, "AD9"))] <- 'AD9'
metadata$sample[which(str_detect(metadata$cells, "AD10"))] <- 'AD10'
metadata$sample[which(str_detect(metadata$cells, "AD20"))] <- 'AD20'
metadata$sample[which(str_detect(metadata$cells, "AD21"))] <- 'AD21'
metadata$sample[which(str_detect(metadata$cells, "NC3"))] <- 'NC3'
metadata$sample[which(str_detect(metadata$cells, "NC7"))] <- 'NC7'
metadata$sample[which(str_detect(metadata$cells, "NC12"))] <- 'NC12'
metadata$sample[which(str_detect(metadata$cells, "NC14"))] <- 'NC14'
metadata$sample[which(str_detect(metadata$cells, "NC15"))] <- 'NC15'
metadata$sample[which(str_detect(metadata$cells, "NC17"))] <- 'NC17'
metadata$sample[which(str_detect(metadata$cells, "NC18"))] <- 'NC18'

metadata <- metadata %>% mutate(disease = case_when(
    sample %in% c("AD1","AD2","AD4","AD5","AD9","AD10",
                  "AD20","AD21") ~ 'AD',
    sample %in% c("NC3","NC7","NC12","NC14","NC15","NC17","NC18") ~ 'NC'))

merged_seurat@meta.data <- metadata
###QC###
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI  < 20000) & 
                            (nGene > 200) &
                            (mitoRatio < 0.20))
filtered_seurat_QC <- JoinLayers(filtered_seurat)
counts <- LayerData(object = filtered_seurat_QC, layer = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
matched_meta <- filtered_seurat@meta.data[colnames(filtered_counts), ]
filtered_seurat <- CreateSeuratObject(counts = filtered_counts,meta.data = matched_meta)
#Filter_out_samples_that_contain_fewer_than_4,000 cells#
#remove(AD6,AD8,AD13,AD19,NC11)#
table(filtered_seurat$sample)
metadata <- filtered_seurat@meta.data
#QC_Visualize#
metadata %>% ggplot(aes(x=sample, fill=sample)) + geom_bar() +
  theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

metadata %>% ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + theme_classic() +
  ylab("Cell density") +geom_vline(xintercept = 500)
#(because_shift_remove_NC16)⬇#
metadata %>% ggplot(aes(color=sample, x=nGene, fill= sample)) + geom_density(alpha = 0.2) + 
  theme_classic() +scale_x_log10() + geom_vline(xintercept = 300)

metadata %>% ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + geom_boxplot() + 
  theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +ggtitle("NCells vs NGenes")

metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +stat_smooth(method=lm) +
  scale_x_log10() + scale_y_log10() + theme_classic() +geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +facet_wrap(~sample)
#Memory_optimization##
rm(AD1,AD2,AD4,AD5,AD9,AD10,AD20,AD21,NC3,NC7,NC12,NC14,NC15,NC17,NC18,
   filtered_seurat_QC,counts,filtered_counts,nonzero,merged_seurat)
gc()
#normalize data#
scRNAlist_merge <- NormalizeData(filtered_seurat)
scRNAlist_merge <- FindVariableFeatures(scRNAlist_merge,nfeatures = 2000) 
scRNAlist_merge <- ScaleData(scRNAlist_merge,vars.to.regress = c('mitoRatio')) 
#run_PCA#
scRNAlist_merge <- RunPCA(scRNAlist_merge,npcs = 50)
print(scRNAlist_merge[["pca"]], dims = 1:5, nfeatures = 5)
#visualize_check_data#
VizDimLoadings(scRNAlist_merge, dims = 1:2, reduction = "pca")
DimPlot(scRNAlist_merge, reduction = "pca")
DimHeatmap(scRNAlist_merge, dims = 1:20, cells = 500, balanced = TRUE)
#remove_Batch_effect_by_harmony#
library(harmony)
scRNA_harmony <- RunHarmony(object = scRNAlist_merge,
                            group.by.vars = "orig.ident", 
                            reduction = "pca", 
                            dims.use = 1:20, 
                            reduction.save = "harmony") 
#joinLayers#
scRNA_harmony[['RNA']] <- JoinLayers(scRNA_harmony[['RNA']])
# chose_0.4_resolution#
scRNA_harmony <- FindNeighbors(scRNA_harmony,reduction = 'harmony',dims = 1:20)
scRNA_harmony <- FindClusters(scRNA_harmony,resolution = seq(from = 0.1,to = 1.0, by = 0.1))
scRNA_harmony <- RunUMAP(scRNA_harmony,dims = 1:20,reduction = 'harmony')
scRNA_harmony <- RunTSNE(scRNA_harmony,dims = 1:20,reduction = 'harmony')
scRNA_harmony$RNA_snn_res.0.4
Idents(scRNA_harmony) <- 'RNA_snn_res.0.4'
DimPlot(scRNA_harmony,reduction = 'umap')
#UMAP_and_TSNE#
umap_integrated_1 <- DimPlot(scRNA_harmony,reduction = 'umap',group.by = 'orig.ident',label = T) 
umap_integrated_2 <- DimPlot(scRNA_harmony,reduction = 'tsne', label = T) 
umap_integrated_3 <- DimPlot(scRNA_harmony,reduction = 'umap', label = T) 
umap_integrated_4 <- DimPlot(scRNA_harmony,reduction = 'umap',group.by = 'disease',label = T) 
umap_integrated_1
umap_integrated_2
umap_integrated_3
umap_integrated_4
# save data 
save(scRNA_harmony,file = 'scRNA_harmony.Rdata')




########Run analysis#########
rm(list = ls())
gc()
#load packages#
library(Seurat)
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(dplyr)

###load_data###
load('scRNA_harmony.Rdata')
scRNA_harmony$RNA_snn_res.0.4
Idents(scRNA_harmony) <- 'RNA_snn_res.0.4'
DimPlot(scRNA_harmony,reduction = 'umap')
###Annotation_was_based_on_marker_genes_reported_in_the_original_publication###
current_levels <- levels(scRNA_harmony)
sorted_levels <- as.character(sort(as.numeric(current_levels)))
Idents(scRNA_harmony) <- factor(Idents(scRNA_harmony), levels = sorted_levels)
levels(scRNA_harmony)
new.cluster.ids <- c("oligodendrocyte_1",   # "MBP"0
                     "astrocytes_1",        # "AQP4"1
                     "excitatory_neuron_1", # "CAMK2A"2
                     "excitatory_neuron_2", # "CAMK2A"3
                     "inhibitory_neuron_1", # "GAD1"
                     "inhibitory_neuron_2", # "GAD1"8
                     "excitatory_neuron_3", # "CAMK2A"6
                     "microglia_1",         # "C3"    7
                     "inhibitory_neuron_3", # "GAD1"8
                     "excitatory_neuron_4", # "CAMK2A"9
                     "excitatory_neuron_5", # "CAMK2A"10
                     "inhibitory_neuron_4", # "GAD1"11
                     "excitatory_neuron_6", # "CAMK2A"12
                     "inhibitory_neuron_5", # "GAD1"
                     "excitatory_neuron_7", # "CAMK2A"
                     "excitatory_neuron_8", # "CAMK2A"
                     "excitatory_neuron_9", # "CAMK2A"
                     "endothelial_cell_1", #  "CLDN5"
                     "astrocytes_2",        # "AQP4"
                     "excitatory_neuron_10", # "CAMK2A"
                     "excitatory_neuron_11", # "CAMK2A"
                     "oligodendrocyte_2",   # "MBP"
                     "excitatory_neuron_12", # "CAMK2A"
                     "excitatory_neuron_13", # "CAMK2A"
                     "inhibitory_neuron_6", # "GAD1"
                     "oligodendrocyte_3",   # "MBP"
                     "astrocytes_3")        # "AQP4"
names(new.cluster.ids) <- levels(scRNA_harmony)
scRNA_harmony <- RenameIdents(scRNA_harmony, new.cluster.ids)
scRNA_harmony<- AddMetaData(object = scRNA_harmony, 
                            metadata = scRNA_harmony@active.ident,
                            col.name = "celltype")
#visualize_Annotation_marker_genes#
cell_mark <- c("ADGRV1","GPC5","RYR3", "AQP4", # astrocytes 
               "ABCB1","EBF1", "CLDN5", # endothelial_cell
               "CBLN2","LDB2", "CAMK2A", # excitatory_neuron
               "LHFPL3","PCDH15","GAD1",   # inhibitory_neuron
               "LRMDA","DOCK8","C3",       # microglia
               "PLP1","ST18","MBP")    # oligodendrocyte
DotPlot(scRNA_harmony, 
        features = cell_mark,
        assay = "RNA",
        cluster.idents = TRUE,  
        scale.by = "size",      
        scale = TRUE,           
        col.min = -2,            
        col.max = 2) +
  coord_flip() +theme_bw() +  labs(x = "Genes", y = "Cell Types") +  
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1, size=8),
        axis.text.y = element_text(size=8)) +
  scale_color_gradient2(low = "#2166AC",mid = "white",high = "#B2182B")

###AUC###
seurat_integrated<-scRNA_harmony
DefaultAssay(seurat_integrated) <- "RNA"
library(readxl)
library("Matrix")
library("AUCell")
terms<- read_excel("list.xlsx")
print(terms)
geneset<- terms[['comorbidity_associated_upregulated_genes']]
geneset <- geneset[!is.na(geneset)]
valid_genes <- intersect(geneset, rownames(seurat_integrated))
if (length(valid_genes) < 5) stop("Too few valid genes).")
gene_sets <- list("comorbidity_associated_upregulated_genes" = valid_genes)
seurat_integrated <- JoinLayers(seurat_integrated)
expr_matrix <- LayerData(seurat_integrated, layer = "data")
expr_matrix <- as(expr_matrix, "dgCMatrix") 
class(expr_matrix)
cells_rankings <- AUCell_buildRankings(expr_matrix,plotStats = FALSE,  splitByBlocks = TRUE)
cells_AUC <- AUCell_calcAUC(gene_sets,cells_rankings,
                            aucMaxRank = ceiling(0.05 * nrow(cells_rankings)),nCores = 4)
seurat_integrated$upregulated <- as.numeric(getAUC(cells_AUC)["comorbidity_associated_upregulated_genes", ])
##visualize_AUC_umap##
umap_df <- FetchData(seurat_integrated, vars = c("umap_1", "umap_2", "upregulated"))
ggplot(umap_df, aes(x = umap_1, y = umap_2, color = upregulated)) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_color_viridis_c(option = "viridis",name = "AUC Score" ) +
  labs(x = "UMAP1", y = "UMAP2", 
       title = "" ) +
  theme_classic()
##visualize_AUC_Bubble_Plot## 
auc_df <- seurat_integrated@meta.data %>%group_by(celltype) %>%
  summarise(auc_mean = mean(upregulated),n_cells = n())
ggplot(auc_df, aes(x = celltype, y = 1,size = n_cells,color = auc_mean)) +
  geom_point(alpha = 0.9) +
  scale_size(range = c(3, 15)) +
  scale_color_gradientn(colors = c("#2166AC", "white", "#B2182B")) +
  theme_bw() +
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),axis.text.x = element_text(angle =90, hjust = 1, size = 10)) +
  labs(x = "Cell Type",color = "AUC Score",size = "Number of Cells",
       title = "AUC Activity Across Cell Types")
##visualize_AUC_BOX_plot##
#chose_cell_type#
target_cells <- c('astrocytes_1',"microglia_1","oligodendrocyte_1","oligodendrocyte_3",
                  "excitatory_neuron_5","endothelial_cell_1") 
seurat_subset <- subset(seurat_integrated, subset = celltype %in% target_cells) 
plot_data <- FetchData(seurat_subset, vars = c("upregulated","disease","celltype"))
#disease_group_in_order#
plot_data$disease <- factor(plot_data$disease, levels = c("AD", "NC"))
#cell_type_in_order#
plot_data$cell_type <- factor(plot_data$celltype,
                              levels = c('astrocytes_1',"microglia_1","oligodendrocyte_1", 
                         "oligodendrocyte_3","excitatory_neuron_5","endothelial_cell_1"))
#visualize#
ggplot(plot_data, aes(x = cell_type,y = upregulated,fill = disease)) +     
  geom_boxplot(position = position_dodge(0.8),width = 0.7,outlier.size = 0.5)+
  scale_fill_manual(values = c("NC" = "#1f77b4", "AD" = "#ff7f0e"),name = " ")+
  labs(x = "Cell Type",y = "AUC Score",title = " ",subtitle = " ") +
  theme_classic() +
  theme(legend.position = "top",legend.title = element_text(size = 16, face = "bold"),  
    legend.text  = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
    axis.text.y = element_text(size = 16),  axis.title.x = element_text(size = 16, face = "bold"),        
    axis.title.y = element_text(size = 16, face = "bold"),plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5))

###visualize_proportions_Bar_plot###
library(forcats)   
library(scales)    
metadata <- seurat_integrated@meta.data
#chose_cell_type#
selected_celltypes <- c('astrocytes_1',"microglia_1","oligodendrocyte_1", "oligodendrocyte_3",
                        "excitatory_neuron_5","endothelial_cell_1")
#Calculate_cell-type_proportions#
total_per_disease <- metadata %>%group_by(disease) %>%summarise(total_cells = n(), .groups = "drop")
celltype_percent <- metadata %>%filter(celltype %in% selected_celltypes) %>%
  group_by(disease, celltype) %>%summarise(n = n(), .groups = "drop") %>%
  left_join(total_per_disease, by = "disease") %>%mutate(percent = 100 * n / total_cells)
#cell_type_in_order#
celltype_percent <- celltype_percent %>%mutate(celltype = factor(celltype, levels = selected_celltypes))
#disease_group_in_order#
celltype_percent <- celltype_percent %>%mutate(disease = factor(disease, levels = c("NC", "AD")))
#visualize#
p <- ggplot(celltype_percent, aes(x = celltype, y = percent, fill = disease)) +
  geom_bar(stat = "identity",position = position_dodge(width = 0.8),color = "black",width = 0.7)+
  geom_text(aes(label = sprintf("%.1f", percent)),position = position_dodge(width = 0.8),
    vjust = -0.4, size = 4.5, color = "black")+
  scale_fill_manual(name = " ",values = c("NC" = "#1f77b4", "AD" ="#ff7f0e")) +
  labs(x = "",y = "Proportion of cells (%)",title = " "
  ) +theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "top",legend.text = element_text(size = 14),    
        legend.title = element_text(size = 15),legend.key.size = unit(0.9, "cm"),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.2))+
  scale_y_continuous(limits = c(0, max(celltype_percent$percent, na.rm = TRUE) * 1.15),
                     expand = expansion(mult = c(0, 0.05)))
p

###donor-aware differential abundance testing###   
#chose_cell_type#
target_celltypes <- c("oligodendrocyte_1",  "astrocytes_1","excitatory_neuron_1","excitatory_neuron_2",
                      "inhibitory_neuron_1", "inhibitory_neuron_2","excitatory_neuron_3",
                      "microglia_1","inhibitory_neuron_3","excitatory_neuron_4", "excitatory_neuron_5", 
                      "inhibitory_neuron_4", "excitatory_neuron_6", "inhibitory_neuron_5", 
                      "excitatory_neuron_7", "excitatory_neuron_8", "excitatory_neuron_9","endothelial_cell_1", 
                      "astrocytes_2", "excitatory_neuron_10", "excitatory_neuron_11", "oligodendrocyte_2",   
                      "excitatory_neuron_12", "excitatory_neuron_13","inhibitory_neuron_6", 
                      "oligodendrocyte_3","astrocytes_3")
#check_data#
md <- scRNA_harmony@meta.data
donor_group <- tapply(md$disease, md$sample, function(x) unique(x)[1])
#run#
one_ct <- function(ct) { n_ct<- tapply(md$celltype == ct, md$sample, sum)
                         n_total  <- tapply(rep(1, nrow(md)), md$sample, sum)
                         prop_ct  <- n_ct / n_total
                         g <- donor_group[names(prop_ct)]
                         grp_levels <- unique(md$disease)
                         if (length(grp_levels) < 2) return(NULL)
                         g1 <- prop_ct[g == grp_levels[1]]
                         g2 <- prop_ct[g == grp_levels[2]]
                         if (sum(!is.na(g1)) < 2 || sum(!is.na(g2)) < 2) {p <- NA_real_} 
                         else {p <- suppressWarnings(wilcox.test(g1, g2, exact=FALSE)$p.value)}
                   data.frame(celltype = ct,n_donors_group1 = sum(!is.na(g1)),
                              n_donors_group2 = sum(!is.na(g2)),median_group1=median(g1, na.rm = TRUE),
                              median_group2= median(g2, na.rm = TRUE),p.value = p,
                              stringsAsFactors = FALSE)}
res_list <- lapply(target_celltypes, one_ct)
res_DA <- do.call(rbind, res_list)
res_DA$FDR <- p.adjust(res_DA$p.value, method = "BH")
write.csv(res_DA, "DA_donor_wilcoxon_SELECTED_celltypes.csv", row.names = FALSE)
print(res_DA)

###CellChat ####
library(CellChat)
library(ggalluvial)
library(patchwork)
library(RColorBrewer)
options(stringsAsFactors = FALSE)
#Prepare_the_data_for_CellChat#
seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])
data.input <- LayerData(seurat_integrated, assay = "RNA", layer = "data")
metadata <- data.frame(cell_type = Idents(seurat_integrated),
                       disease = seurat_integrated$disease,    
                       row.names = colnames(seurat_integrated))
print(unique(metadata$disease))
#Cell–cell_communication_is_shown_for_AD_only#
selected_diseases <- c('AD')  
disease_cells <- rownames(metadata)[metadata$disease %in% selected_diseases]
#chose_cell_type#
selected_celltypes <- c('astrocytes_1',"microglia_1","oligodendrocyte_1", "oligodendrocyte_3",
                        "excitatory_neuron_5","endothelial_cell_1")

celltype_cells <- rownames(metadata)[metadata$cell_type %in% selected_celltypes]
selected_cells <- intersect(disease_cells, celltype_cells)
data.input.subset <- data.input[, selected_cells]
metadata.subset <- metadata[selected_cells, , drop = FALSE] 
cellchat <- createCellChat(object = data.input.subset)
cell_type_df <- data.frame(cell_type = metadata.subset$cell_type,
                           row.names = rownames(metadata.subset))
full_meta_df <- data.frame(cell_type = metadata.subset$cell_type,
                           disease = metadata.subset$disease,
                           row.names = rownames(metadata.subset))
cellchat <- addMeta(cellchat, meta = full_meta_df)
cellchat <- setIdent(cellchat, ident.use = "cell_type")
cellchat@idents <- factor(as.character(cellchat@meta$cell_type), 
                          levels = selected_celltypes)
print(levels(cellchat@idents))
groupSize <- as.numeric(table(cellchat@idents))
print(groupSize)
print(table(cellchat@meta$disease))
#Set_database#
CellChatDB <- CellChatDB.human 
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
#run_analysis#
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)  
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
#quantitative_centrality_metrics#
library(igraph)
W <- cellchat@net$weight
if (is.null(W) || all(W == 0)) stop("CellChat no network")
g <- graph_from_adjacency_matrix(W, mode = "directed", weighted = TRUE, diag = FALSE)
outgoing <- strength(g, mode = "out", weights = E(g)$weight)
incoming <- strength(g, mode = "in", weights = E(g)$weight)
bet <- betweenness(g, directed = TRUE, weights = 1/E(g)$weight)
clo <- closeness(g, mode = "all", weights = 1/E(g)$weight)
centrality <- data.frame(celltype   = rownames(W),
                         outgoing   = outgoing,
                         incoming   = incoming,
                         betweenness = bet,
                         closeness   = clo)
write.csv(centrality, "CellChat_global_centrality.csv", row.names = FALSE)
centrality

#visualize#
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
print("pathways")
pathways <- cellchat@netP$pathways
print(pathways)
#visualize#
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
colors <- RColorBrewer::brewer.pal(length(selected_celltypes), "Set1")
names(colors) <- selected_celltypes
#visualize_communication_weight _matrix_plot#
par(mar = c(0, 0, 2.5, 0))  #
netVisual_circle(cellchat@net$weight,
                 color.use=colors,
                 vertex.weight=groupSize,
                 weight.scale=TRUE,
                 edge.weight.max=max(cellchat@net$weight),
                 title.name="Communication Weight Matrix",
                 vertex.label.cex = 1.4)
#visualize_communication_weight _heatmap_plot#
netVisual_heatmap(cellchat,
                  signaling = NULL,
                  measure = "weight",
                  color.heatmap = "Reds",
                  title.name = "Communication Weight Heatmap",
                  font.size = 14,
                  font.size.title = 16)
#Prioritization_of_sender_and_receiver_signals_of_interest#
sender_cell <- "astrocytes_1"
receiver_cell <- "oligodendrocyte_3"
pathway_comm <- subsetCommunication(cellchat, slot.name = "netP",
                                    sources.use = sender_cell,
                                    targets.use = receiver_cell)
pathway_comm_sorted <- pathway_comm[order(-pathway_comm$prob), ]
print("Top 10 signaling pathways:")
print(head(pathway_comm_sorted, 10))
#Visualize_specific_pathways#
colors <- brewer.pal(length(selected_celltypes), "Dark2")
names(colors) <- selected_celltypes
pathways.show<-'NRXN'
netVisual_aggregate(cellchat, 
                    signaling = pathways.show,
                    layout = "chord",
                    vertex.size = groupSize,
                    color.use = colors,
                    title.name = paste(pathway, "Signaling"),
                    arrow.size = 0.02,
                    edge.width.max = 15)

###overlap_between_SUGS_and_ligand–receptor###
caup_genes <- unique(valid_genes)
caup_genes <- caup_genes[!is.na(caup_genes)]
UP <- function(x) unique(toupper(na.omit(x)))
split_genes <- function(x) {if (is.null(x)) return(character(0))
  x <- na.omit(x)
  unlist(strsplit(x, "\\s*[+|_,/; ]\\s*"))}
#check_overlap_in_Global#
comm_all <- subsetCommunication(cellchat, slot.name = "netP")
lr_all <- unique(c(split_genes(comm_all$ligand), split_genes(comm_all$receptor)))
overlap_global <- intersect(UP(caup_genes), UP(lr_all))
cat("Global overlap size =", length(overlap_global), "\n")
print(sort(overlap_global))
#check_overlap_in_one_sender#
sender_only <- "astrocytes_1"
comm_sender <- subsetCommunication(cellchat, slot.name = "netP", sources.use = sender_only)
lr_sender <- unique(c(split_genes(comm_sender$ligand), split_genes(comm_sender$receptor)))
overlap_sender <- intersect(UP(caup_genes), UP(lr_sender))
cat(sprintf("Overlap (sender=%s) size = %d\n", sender_only, length(overlap_sender)))
print(sort(overlap_sender))
#check_overlap_in_one_receiver#
receiver_only <- "oligodendrocyte_1"
comm_receiver <- subsetCommunication(cellchat, slot.name = "netP", targets.use = receiver_only)
lr_receiver <- unique(c(split_genes(comm_receiver$ligand), split_genes(comm_receiver$receptor)))
overlap_receiver <- intersect(UP(caup_genes), UP(lr_receiver))
cat(sprintf("Overlap (receiver=%s) size = %d\n", receiver_only, length(overlap_receiver)))
print(sort(overlap_receiver))
#check_overlap_in_one_receiver#
receiver_only <- "oligodendrocyte_3"
comm_receiver <- subsetCommunication(cellchat, slot.name = "netP", targets.use = receiver_only)
lr_receiver <- unique(c(split_genes(comm_receiver$ligand), split_genes(comm_receiver$receptor)))
overlap_receiver <- intersect(UP(caup_genes), UP(lr_receiver))
cat(sprintf("Overlap (receiver=%s) size = %d\n", receiver_only, length(overlap_receiver)))
print(sort(overlap_receiver))

###Find_top10_Markers_for_every_cell-subtype###
pbmc.markers <- FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
topmaker<-pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(topmaker, file = "topmaker.csv", row.names = FALSE)

###pseudobulk###
target_celltypes <-c("oligodendrocyte_1",  "astrocytes_1","excitatory_neuron_1","excitatory_neuron_2",
                     "inhibitory_neuron_1", "inhibitory_neuron_2","excitatory_neuron_3",
                     "microglia_1","inhibitory_neuron_3","excitatory_neuron_4", "excitatory_neuron_5", 
                     "inhibitory_neuron_4", "excitatory_neuron_6", "inhibitory_neuron_5", 
                     "excitatory_neuron_7", "excitatory_neuron_8", "excitatory_neuron_9","endothelial_cell_1", 
                     "astrocytes_2", "excitatory_neuron_10", "excitatory_neuron_11", "oligodendrocyte_2",   
                     "excitatory_neuron_12", "excitatory_neuron_13","inhibitory_neuron_6", 
                     "oligodendrocyte_3","astrocytes_3")

suppressWarnings(suppressMessages({
  edgeR_ok <- requireNamespace("edgeR", quietly = TRUE) 
  library(Matrix)}))
.assay <- DefaultAssay(scRNA_harmony)
get_counts <- function(obj, assay) {
  lyr_names <- tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) character(0))
  if ("counts" %in% lyr_names) {
    return(SeuratObject::GetAssayData(obj, assay = assay, layer = "counts"))}
  if ("data"   %in% lyr_names) {
    m <- SeuratObject::GetAssayData(obj, assay = assay, layer = "data")
    m <- as(m, "dgCMatrix")
    return(m)}
  m <- tryCatch(SeuratObject::GetAssayData(obj, assay = assay), error = function(e) NULL)
  if (!is.null(m)) return(as(m, "dgCMatrix"))
  stop("Seurat v5 no counts/data")}
cm <- get_counts(scRNA_harmony, .assay)  
md <- scRNA_harmony@meta.data
out_dir <- "PB_validation_SELECTED"
dir.create(out_dir, showWarnings = FALSE)
run_one_ct <- function(ct) {
  cells_ct <- rownames(md)[md$celltype == ct]
  if (length(cells_ct) == 0) { message("jump over：", ct, "nocells"); return(invisible(NULL)) }
  donors <- unique(md$sample[md$celltype == ct])
  grp_by_donor <- tapply(md$disease, md$sample, function(x) unique(x)[1])[donors]
  grp_levels   <- unique(md$disease)
  if (length(grp_levels) < 2) { message("jump over：", ct, "groupproblem"); return(invisible(NULL)) }
  donors_g1 <- donors[grp_by_donor == grp_levels[1]]
  donors_g2 <- donors[grp_by_donor == grp_levels[2]]
  if (length(donors_g1) < 2 || length(donors_g2) < 2) {
    message("jump over：", ct, "sample<2 per group"); return(invisible(NULL))}
  pb_cols <- lapply(donors, function(dn) {
    cs <- rownames(md)[md$celltype == ct & md$sample == dn]
    if (length(cs) == 0) return(Matrix(0, nrow = nrow(cm), ncol = 1, sparse = TRUE))
    Matrix::rowSums(cm[, cs, drop = FALSE])
  })
  pb <- do.call(cbind, pb_cols)
  colnames(pb) <- donors
  if (edgeR_ok) {y <- edgeR::DGEList(counts = pb,
                      samples = data.frame(sample = donors,
                                           group  = grp_by_donor[donors],
                                           row.names = donors))
    keep <- rowSums(y$counts) > 0
    if (sum(keep) < 10) { message("jump over", ct, "genes < 10"); return(invisible(NULL)) }
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)
    design <- model.matrix(~ y$samples$group)
    y <- edgeR::estimateDisp(y, design)
    fit <- edgeR::glmQLFit(y, design)
    qlf <- edgeR::glmQLFTest(fit, coef = 2)
    tt  <- edgeR::topTags(qlf, n = Inf)$table
    tt$gene <- rownames(tt)
    out <- file.path(out_dir, paste0("PB_edgeR_DEGs_", make.names(ct), ".csv"))
    write.csv(tt, out, row.names = FALSE)
    message("[edgeR]：", out)} 
  invisible(TRUE)}
invisible(lapply(target_celltypes, run_one_ct))

