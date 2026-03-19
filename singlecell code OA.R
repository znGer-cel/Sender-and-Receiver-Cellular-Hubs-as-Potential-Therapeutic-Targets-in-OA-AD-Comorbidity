########Data Processing, Filtering, and Integration#########
#load packages#
library(Seurat)
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(stringr)
library(readxl)
###load_data###
#load_data+creat_seurat#
rm(list = ls())
Normal2 <- Read10X(data.dir = "GSE255460/C2/") 
Normal2 <- CreateSeuratObject(counts = Normal2, project = "Normal2", min.cells = 3, min.features = 200)
Normal3 <- Read10X(data.dir = "GSE255460/C3/") 
Normal3 <- CreateSeuratObject(counts = Normal3, project = "Normal3", min.cells = 3, min.features = 200)
OA1 <- Read10X(data.dir = "GSE255460/OA1/") 
OA1 <- CreateSeuratObject(counts = OA1, project = "OA1", min.cells = 3, min.features = 200)
OA1_1 <- Read10X(data.dir = "GSE255460/OA1_1/") 
OA1_1 <- CreateSeuratObject(counts = OA1_1, project = "OA1_1", min.cells = 3, min.features = 200)
OA2 <- Read10X(data.dir = "GSE255460/OA2/") 
OA2 <- CreateSeuratObject(counts = OA2, project = "OA2", min.cells = 3, min.features = 200)
OA2_1 <- Read10X(data.dir = "GSE255460/OA2_1/") 
OA2_1 <- CreateSeuratObject(counts = OA2_1, project = "OA2_1", min.cells = 3, min.features = 200)
OA3 <- Read10X(data.dir = "GSE255460/OA3/") 
OA3 <- CreateSeuratObject(counts = OA3, project = "OA3", min.cells = 3, min.features = 200)
OA3_1 <- Read10X(data.dir = "GSE255460/OA3_1/") 
OA3_1 <- CreateSeuratObject(counts = OA3_1, project = "OA3_1", min.cells = 3, min.features = 200)
OA4 <- Read10X(data.dir = "GSE255460/OA4/") 
OA4 <- CreateSeuratObject(counts = OA4, project = "OA4", min.cells = 3, min.features = 200)
OA4_1 <- Read10X(data.dir = "GSE255460/OA4_1/") 
OA4_1 <- CreateSeuratObject(counts = OA4_1, project = "OA4_1", min.cells = 3, min.features = 200)
OA5 <- Read10X(data.dir = "GSE255460/OA5/") 
OA5 <- CreateSeuratObject(counts = OA5, project = "OA5", min.cells = 3, min.features = 200)
OA5_1 <- Read10X(data.dir = "GSE255460/OA5_1/") 
OA5_1 <- CreateSeuratObject(counts = OA5_1, project = "OA5_1", min.cells = 3, min.features = 200)
OA6 <- Read10X(data.dir = "GSE255460/OA6/") 
OA6 <- CreateSeuratObject(counts = OA6, project = "OA6", min.cells = 3, min.features = 200)
OA6_1 <- Read10X(data.dir = "GSE255460/OA6_1/") 
OA6_1 <- CreateSeuratObject(counts = OA6_1, project = "OA6_1", min.cells = 3, min.features = 200)
OA7 <- Read10X(data.dir = "GSE255460/OA7/") 
OA7 <- CreateSeuratObject(counts = OA7, project = "OA7", min.cells = 3, min.features = 200)
OA8 <- Read10X(data.dir = "GSE255460/OA8/") 
OA8 <- CreateSeuratObject(counts = OA8, project = "OA8", min.cells = 3, min.features = 200)
OA8_1 <- Read10X(data.dir = "GSE255460/OA8_1/") 
OA8_1 <- CreateSeuratObject(counts = OA8_1, project = "OA8_1", min.cells = 3, min.features = 200)
#merge_data#
merged_seurat <- merge(Normal2, y = c(Normal3,OA1,OA1_1,OA2,OA2_1,OA3,OA3_1,
                                      OA4,OA4_1,OA5,OA5_1,OA6,OA6_1,OA7,OA8,OA8_1),
                       add.cell.ids = c("Normal2","Normal3","OA1","OA1_1",
                                        "OA2","OA2_1","OA3","OA3_1","OA4","OA4_1","OA5","OA5_1",
                                        "OA6","OA6_1","OA7","OA8","OA8_1"))
#count_UMI_MT#
merged_seurat$Novelty_score <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
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
metadata$sample[which(str_detect(metadata$cells, "Normal2"))] <- 'Normal2'
metadata$sample[which(str_detect(metadata$cells, "Normal3"))] <- 'Normal3'
metadata$sample[which(str_detect(metadata$cells, "OA1"))] <- 'OA1'
metadata$sample[which(str_detect(metadata$cells, "OA1_1"))] <- 'OA1_1'
metadata$sample[which(str_detect(metadata$cells, "OA2"))] <- 'OA2'
metadata$sample[which(str_detect(metadata$cells, "OA2_1"))] <- 'OA2_1'
metadata$sample[which(str_detect(metadata$cells, "OA3"))] <- 'OA3'
metadata$sample[which(str_detect(metadata$cells, "OA3_1"))] <- 'OA3_1'
metadata$sample[which(str_detect(metadata$cells, "OA4"))] <- 'OA4'
metadata$sample[which(str_detect(metadata$cells, "OA4_1"))] <- 'OA4_1'
metadata$sample[which(str_detect(metadata$cells, "OA5"))] <- 'OA5'
metadata$sample[which(str_detect(metadata$cells, "OA5_1"))] <- 'OA5_1'
metadata$sample[which(str_detect(metadata$cells, "OA6"))] <- 'OA6'
metadata$sample[which(str_detect(metadata$cells, "OA6_1"))] <- 'OA6_1'
metadata$sample[which(str_detect(metadata$cells, "OA7"))] <- 'OA7'
metadata$sample[which(str_detect(metadata$cells, "OA8"))] <- 'OA8'
metadata$sample[which(str_detect(metadata$cells, "OA8_1"))] <- 'OA8_1'

metadata <- metadata %>%
  mutate(disease = case_when(
    sample %in% c("Normal2","Normal3") ~ 'Normal',
    sample %in% c("OA1","OA1_1","OA2","OA2_1","OA3","OA3_1","OA4","OA4_1",
                  "OA5","OA5_1","OA6","OA6_1","OA7","OA8","OA8_1") ~ 'OA'
  ))
merged_seurat@meta.data <- metadata
###QC###
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 200) &
                            (nGene <= 4000) &
                            (Novelty_score > 0.80) &
                            (mitoRatio < 0.10))
filtered_seurat_QC <- JoinLayers(filtered_seurat)
counts <- LayerData(object = filtered_seurat_QC, layer = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
matched_meta <- filtered_seurat@meta.data[colnames(filtered_counts), ]
filtered_seurat <- CreateSeuratObject(counts = filtered_counts,meta.data = matched_meta)
#Filter_out_samples_that_contain_fewer_than_4,000 cells#
#remove(Normal_1,OA7)#
table(filtered_seurat$sample)
metadata <- filtered_seurat@meta.data
#QC_Visualize#
metadata %>% ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(plot.title = element_text(hjust=0.5, face="bold"))+
  ggtitle("NCells")
metadata %>% ggplot(aes(color=sample, x=nUMI, fill= sample))+ 
  geom_density(alpha = 0.2) + scale_x_log10()+ 
  theme_classic() +ylab("Cell density")+geom_vline(xintercept = 500)
metadata %>% ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + theme_classic() +
  scale_x_log10() + geom_vline(xintercept = 300)
metadata %>% ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +ggtitle("NCells vs NGenes")
metadata %>% ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +scale_x_log10() + scale_y_log10() + theme_classic() +
  geom_vline(xintercept = 500) +geom_hline(yintercept = 250) +facet_wrap(~sample)
##Memory_optimization##
rm(Normal2,Normal3,OA1,OA1_1,OA2,OA2_1,OA3,OA3_1,OA4,OA4_1,OA5,OA5_1,OA6,OA6_1,
   OA7,OA8,OA8_1,filtered_seurat_QC,counts,filtered_counts,nonzero,merged_seurat)
gc()
#normalize data#
scRNAlist_merge <- NormalizeData(filtered_seurat)
scRNAlist_merge <- FindVariableFeatures(scRNAlist_merge,nfeatures = 2000) 
scRNAlist_merge <- ScaleData(scRNAlist_merge,vars.to.regress = c('mitoRatio')) 
#run_PCA#
scRNAlist_merge <- RunPCA(scRNAlist_merge, npcs = 50)
print(scRNAlist_merge[["pca"]], dims = 1:40, nfeatures = 5)
#visualize_check_data#
VizDimLoadings(scRNAlist_merge, dims = 1:5, reduction = "pca")
DimPlot(scRNAlist_merge, reduction = "pca")
DimHeatmap(scRNAlist_merge, dims = 1:20, cells = 500, balanced = TRUE)
#remove_Batch_effect_by_harmony#
library(harmony)
scRNA_harmony <- RunHarmony(object = scRNAlist_merge,
                            group.by.vars ="sample", 
                            reduction = "pca", 
                            dims.use = 1:30, 
                            reduction.save = "harmony") 
#joinLayers#
scRNA_harmony[['RNA']] <- JoinLayers(scRNA_harmony[['RNA']])
#chose_0.4_resolution#
scRNA_harmony <- FindNeighbors(scRNA_harmony,reduction = 'harmony',dims = 1:30)
scRNA_harmony <- FindClusters(scRNA_harmony,resolution = seq(from = 0.1,to = 1, by = 0.1))
scRNA_harmony <- RunUMAP(scRNA_harmony,dims = 1:30,reduction = 'harmony')
#scRNA_harmony <- RunTSNE(scRNA_harmony,dims = 1:30,reduction = 'harmony')
scRNA_harmony$RNA_snn_res.0.7
Idents(scRNA_harmony) <- 'RNA_snn_res.0.7'
DimPlot(scRNA_harmony,reduction = 'umap')
#UMAP_and_TSNE#
umap_integrated_1 <- DimPlot(scRNA_harmony,reduction = 'umap',group.by = 'sample',label = T) 
#umap_integrated_2 <- DimPlot(scRNA_harmony,reduction = 'tsne', label = T) 
umap_integrated_3 <- DimPlot(scRNA_harmony,reduction = 'umap', label = T) 
umap_integrated_4 <- DimPlot(scRNA_harmony,reduction = 'umap',group.by = 'disease',label = T) 
umap_integrated_1
#umap_integrated_2
umap_integrated_3
umap_integrated_4
# save data #
save(scRNA_harmony,file = 'scRNA_harmony2.Rdata')


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
load('scRNA_harmony2.Rdata')
scRNA_harmony$RNA_snn_res.0.7
Idents(scRNA_harmony) <- 'RNA_snn_res.0.7'
DimPlot(scRNA_harmony,reduction = 'umap')
umap_integrated_3 <- DimPlot(scRNA_harmony,reduction = 'umap', label = T,label.size = 6)+NoLegend()
umap_integrated_3 
###Annotation_was_based_on_marker_genes_reported_in_the_original_publication###
current_levels <- levels(scRNA_harmony)
sorted_levels <- as.character(sort(as.numeric(current_levels)))
Idents(scRNA_harmony) <- factor(Idents(scRNA_harmony), levels = sorted_levels)
levels(scRNA_harmony)
new.cluster.ids <- c("RepC_1",         # "CILP2","CILP","OGN"     0
                     "preHTC_1",       # "PRG4","ABI3BP","CRTAC1' 1
                     "RegC_1",         # "CHI3L1","CHI3L2"        2
                     "ProC_1",         # "C11orf96","BMP2","HMGA1"3
                     "FC_1",           # "MMP2","COL1A2","COL1A1" 4
                     "EC_1",           # "FRZB","CYTL1"           5
                     "HomC_1",         # "HSPA1A","HSPA1B","HSPA6","DDIT3","JUN"6
                     "HomC_2",         # "HSPA1A","HSPA1B","HSPA6","DDIT3","JUN"7
                     "preHTC_2",       # "PRG4","ABI3BP","CRTAC1'     8
                     "preFC_1",        # "PLCG2","COL27A1","WWP2"     9
                     "ProC_2",         # "C11orf96","BMP2","HMGA1"3  10
                     "HTC_1",          # "COL10A1","SPP1"            11
                     "InfC_1",         # "CD74","GPR183"             12
                     "FC_2",           # "MMP2","COL1A2","COL1A1"    13
                     "preInfC_1" )     # "IFI16","IFI27              14
names(new.cluster.ids) <- levels(scRNA_harmony)
scRNA_harmony <- RenameIdents(scRNA_harmony, new.cluster.ids)
scRNA_harmony<- AddMetaData(object = scRNA_harmony, 
                            metadata = scRNA_harmony@active.ident,
                            col.name = "celltype")
#visualize_Annotation_marker_genes#
cell_mark <- c("HSPA1A","HSPA1B","HSPA6","DDIT3","JUN",    #HomC (homeostasis chondrocytes)
               "CHRDL2","FRZB","CYTL1",   # EC (effector chondrocytes)
               "COL10A1","IBSP","SPP1",   # HTC
               "CILP2","CILP","OGN",      # RepC (reparative chondrocytes)
               "C11orf96","BMP2","HMGA1", # ProC (proliferation chondrocytes)
               "PRG4","ABI3BP","CRTAC1",  # preHTC 
               "MMP2","COL1A2","COL1A1",  # FC (fibrocartilage chondrocytes)
               "PLCG2","COL27A1","WWP2",  # preFC 
               "CHI3L1","CHI3L2",         # RegC (regulator chondrocytes)
               "CXCL8","CD74","GPR183",   # InfC 
               "IFI16","IFI27")          # preInfC
DotPlot(scRNA_harmony, 
        features = cell_mark,
        assay = "RNA",
        cluster.idents = TRUE,  
        scale.by = "size",      
        scale = TRUE,           
        col.min = -2,            
        col.max = 2) +coord_flip()+theme_bw() +
  labs(x = "Genes", y = "Cell Types") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1, size=8),
        axis.text.y = element_text(size=8)) +
  scale_color_gradient2(low = "#2166AC",mid = "white",high = "#B2182B")


###AUC###
library(readxl)
seurat_integrated<-scRNA_harmony
terms<- read_excel("list.xlsx")
print(terms)
geneset<- terms[['comorbidity_associated_upregulated_genes']]
geneset <- geneset[!is.na(geneset)]
library("Matrix")
library("AUCell")
valid_genes <- intersect(geneset, rownames(seurat_integrated))
if (length(valid_genes) < 5) stop("Too few valid genes).")
gene_sets <- list("comorbidity_associated_upregulated_genes" = valid_genes)
seurat_integrated <- JoinLayers(seurat_integrated)
expr_matrix <- LayerData(seurat_integrated, layer = "data")
expr_matrix <- as(expr_matrix, "dgCMatrix")              
cells_rankings <- AUCell_buildRankings(expr_matrix,plotStats = FALSE,  splitByBlocks = TRUE)
cells_AUC <- AUCell_calcAUC(gene_sets,cells_rankings,
                            aucMaxRank = ceiling(0.05 * nrow(cells_rankings)),nCores = 4)
seurat_integrated$upregulated <- as.numeric(getAUC(cells_AUC)["comorbidity_associated_upregulated_genes", ])
##visualize_AUC_umap##
umap_df <- FetchData(seurat_integrated, vars = c("umap_1", "umap_2", "upregulated"))
ggplot(umap_df, aes(x = umap_1, y = umap_2, color = upregulated)) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_color_viridis_c(option = "viridis",name = "AUC Score" ) +
  labs(x = "UMAP1", y = "UMAP2",title = " " )+theme_classic()

##visualize_AUC_Bubble_Plot## 
auc_df <- seurat_integrated@meta.data %>%
  group_by(celltype) %>%
  summarise(auc_mean = mean(upregulated),n_cells = n())
ggplot(auc_df, aes(x = celltype, y = 1,size = n_cells,color = auc_mean)) +
  geom_point(alpha = 0.9)+scale_size(range = c(3, 15)) +
  scale_color_gradientn(colors = c("#2166AC", "white", "#B2182B"))+theme_bw() +
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  labs(x ="Cell Type",color = "AUC Score",size = "Number of Cells",title = "AUC Activity Across Cell Types")

##visualize_AUC_BOX_plot##
#chose_cell_type#
target_cells <- c("preHTC_1", "InfC_1","FC_1") 
seurat_subset <- subset(seurat_integrated, subset = celltype %in% target_cells)  
plot_data <- FetchData(seurat_subset, vars = c("upregulated", "disease", "celltype"))
#disease_group_in_order#
plot_data$disease <- factor(plot_data$disease, levels = c("OA", "Normal"))
#cell_type_in_order#
plot_data$cell_type <- factor(plot_data$celltype,levels = c("preHTC_1", "InfC_1", "RepC_1","FC_1"))
#visualize#
ggplot(plot_data,aes(x = cell_type, y = upregulated,fill = disease)) +     
  geom_boxplot(position = position_dodge(0.8), width = 0.7,outlier.size = 0.5) +
  scale_fill_brewer(palette = "Set1",name = "") +
  labs(x = "Cell Type",y = "AUC Score",title = " ",subtitle = " ") +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_text(size = 16, face = "bold"),  
        legend.text  = element_text(size = 16),                 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
        axis.text.y = element_text(size = 16),                         
        axis.title.x = element_text(size = 16, face = "bold"),         
        axis.title.y = element_text(size = 16, face = "bold"),         
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5))

###visualize_proportions_Bar_plot###
library(forcats)   
library(scales)    
metadata <- seurat_integrated@meta.data
#chose_cell_type#
selected_celltypes <- c("RepC_1","preHTC_1","RegC_1","ProC_1","FC_1","EC_1",           
                        "HomC_1","HomC_2","preHTC_2", "preFC_1","ProC_2",         
                        "HTC_1","InfC_1","FC_2","preInfC_1")
#Calculate_cell-type_proportions#
total_per_disease <- metadata %>%group_by(disease) %>%summarise(total_cells = n(), .groups = "drop")
celltype_percent <- metadata %>%filter(celltype %in% selected_celltypes) %>%
  group_by(disease, celltype) %>%summarise(n = n(), .groups = "drop") %>%
  left_join(total_per_disease, by = "disease") %>%mutate(percent = 100 * n / total_cells)
#cell_type_in_order#
celltype_percent <- celltype_percent %>%mutate(celltype = factor(celltype, levels = selected_celltypes))
#disease_group_in_order#
celltype_percent <- celltype_percent %>%mutate(disease = factor(disease, levels = c("Normal", "OA")))
#visualize#
ggplot(celltype_percent, aes(x = celltype, y = percent, fill = disease)) +
  geom_bar(stat = "identity",position = position_dodge(width = 0.8),color = "black",width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", percent)),position = position_dodge(width = 0.8),
            vjust = -0.4, size = 5, color = "black") +
  scale_fill_manual(name = "",values = c("Normal" = "#1f77b4", "OA" ="#d62728"))+
  labs(x = "",y = "Proportion of cells (%)",title = " ") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "top",legend.text = element_text(size = 14),    
        legend.title = element_text(size = 15),    
        legend.key.size = unit(0.9, "cm"),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.2))+
  scale_y_continuous(limits = c(0, max(celltype_percent$percent, na.rm = TRUE) * 1.15),
                     expand = expansion(mult = c(0, 0.05)))

###donor-aware differential abundance testing###   
#chose_cell_type#
target_celltypes <- c("RepC_1","preHTC_1","RegC_1","ProC_1","FC_1","EC_1",           
                      "HomC_1","HomC_2","preHTC_2", "preFC_1","ProC_2",         
                      "HTC_1","InfC_1","FC_2","preInfC_1")
#check_data#
md <- scRNA_harmony@meta.data
donor_group <- tapply(md$disease, md$sample, function(x) unique(x)[1])
#run#
one_ct <- function(ct) {n_ct <- tapply(md$celltype == ct, md$sample, sum)  
n_total<- tapply(rep(1, nrow(md)), md$sample, sum)   
prop_ct<- n_ct / n_total                             
g <- donor_group[names(prop_ct)]
grp_levels <- unique(md$disease)
if (length(grp_levels) < 2) return(NULL)
g1 <- prop_ct[g == grp_levels[1]]
g2 <- prop_ct[g == grp_levels[2]]
if (sum(!is.na(g1)) < 2 || sum(!is.na(g2)) < 2) {
  p <- NA_real_} else {p <- suppressWarnings(wilcox.test(g1, g2, exact = FALSE)$p.value)}
data.frame(celltype = ct,n_donors_group1 = sum(!is.na(g1)),n_donors_group2 = sum(!is.na(g2)),
           median_group1   = median(g1, na.rm = TRUE),median_group2   = median(g2, na.rm = TRUE),
           p.value = p,stringsAsFactors = FALSE)}
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
#Cell–cell_communication_is_shown_for_OA_only#
selected_diseases <- c('OA')  
disease_cells <- rownames(metadata)[metadata$disease %in% selected_diseases]
#chose_cell_type#
selected_celltypes <- c("RepC_1","preHTC_1","preHTC_2","ProC_1","FC_1","FC_2","InfC_1","preInfC_1")
celltype_cells <- rownames(metadata)[metadata$cell_type %in% selected_celltypes]
selected_cells <- intersect(disease_cells, celltype_cells)
data.input.subset <- data.input[, selected_cells]
metadata.subset <- metadata[selected_cells, , drop = FALSE] 
cellchat <- createCellChat(object = data.input.subset)
cell_type_df <- data.frame(cell_type = metadata.subset$cell_type,row.names = rownames(metadata.subset))
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
centrality <- data.frame(celltype=rownames(W),outgoing=outgoing,
                         incoming=incoming,betweenness=bet,closeness=clo)
write.csv(centrality, "CellChat_global_centrality.csv", row.names = FALSE)
centrality
#visualize#
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
print("pathways")
pathways <- cellchat@netP$pathways
print(pathways)
#visualize#
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
sender_cell <- "FC_1"
receiver_cell <- "InfC_1"
pathway_comm <- subsetCommunication(cellchat,slot.name = "netP",
                                    sources.use = sender_cell,targets.use = receiver_cell)
pathway_comm_sorted <- pathway_comm[order(-pathway_comm$prob), ]
print("Top 10 signaling pathways:")
print(head(pathway_comm_sorted, 10))
#Visualize_specific_pathways#
colors <- brewer.pal(length(selected_celltypes), "Dark2")
names(colors) <- selected_celltypes
pathways.show<-'COLLAGEN'
netVisual_aggregate(cellchat, 
                    signaling = pathways.show,
                    layout = "chord",
                    vertex.size = groupSize,
                    color.use = NULL,
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
lr_all <- unique(c(split_genes(comm_all$ligand), split_genes(comm_all$receptor)))
#check_overlap_in_Global#
comm_all <- subsetCommunication(cellchat, slot.name = "netP")
lr_all <- unique(c(split_genes(comm_all$ligand), split_genes(comm_all$receptor)))
overlap_global <- intersect(UP(caup_genes), UP(lr_all))
cat("Global overlap size =", length(overlap_global), "\n")
print(sort(overlap_global))
#check_overlap_in_one_sender#
sender_only <- "FC_1"
comm_sender <- subsetCommunication(cellchat, slot.name = "netP", sources.use = sender_only)
lr_sender <- unique(c(split_genes(comm_sender$ligand), split_genes(comm_sender$receptor)))
overlap_sender <- intersect(UP(caup_genes), UP(lr_sender))
cat(sprintf("Overlap (sender=%s) size = %d\n", sender_only, length(overlap_sender)))
print(sort(overlap_sender))

###Find_top10_Markers_for_every_cell-subtype###
pbmc.markers <- FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
topmaker<-pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(topmaker, file = "topmaker.csv", row.names = FALSE)


###pseudobulk###
target_celltypes <-c("RepC_1","preHTC_1","RegC_1","ProC_1","FC_1","EC_1",           
                     "HomC_1","HomC_2","preHTC_2", "preFC_1","ProC_2",         
                     "HTC_1","InfC_1","FC_2","preInfC_1")

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
