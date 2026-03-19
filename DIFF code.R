#####****Different_analysis****#####
####for_AD####
##load packages##
library(data.table)  
library(dplyr)      
library(ggplot2)     
library(DESeq2)      
library(edgeR)       
library(limma)
library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggforce)
library(ggrepel)
#load_data#
ADdata <- read_excel("AD_Microarray_data.xlsx",col_names = TRUE)
#pre-differential_analysis_conduct#
ADdata_mean<-aggregate(.~SYMBOL,mean,data=ADdata)
rownames(ADdata_mean) <- ADdata_mean$SYMBOL
ADdata_mean <- ADdata_mean[ , -c(1)]
#grouping#
expr_AD <- ADdata_mean[, 1:28]
expr_Norm <- ADdata_mean[, 29:50]
ADdata_group <- cbind(expr_AD, expr_Norm)  
group4 <- c(rep('AD',ncol(expr_AD)), rep('Norm',ncol(expr_Norm)))
group4 <- factor(group4, levels = c("AD",'Norm'))
#Data filtering#
keep_abundance <- rowSums(ADdata_group) >= 5
ADdata_group_filtered <- ADdata_group[keep_abundance, ]
variances <- apply(ADdata_group_filtered, 1, var) 
var_threshold <- quantile(variances, probs = 0.15) 
keep_variance <- variances >= var_threshold 
ADdata_group_filtered <- ADdata_group_filtered[keep_variance, ]
#*#BOXPOT_check_data#*#
any(ADdata_group_filtered < 0)
group_colors <- c("#FF9999", "#99CCFF")
colors <- rep(group_colors, times = c(28, 22)) 
par(cex = 0.7)
if(ncol(ADdata_group_filtered) > 40) par(cex = 0.5)
boxplot(as.matrix(ADdata_group_filtered), las = 2, col = colors)
#*#PCA#*#
after <- t(ADdata_group_filtered)
pr.out<-prcomp(after)
pcadat<-pr.out
plotdat <- as.data.frame(pcadat$x[,1:2])
rownames(plotdat)
plotdat$group <- c('AD','AD','AD','AD','AD','AD','AD','AD','AD','AD','AD','AD',
                   'AD','AD','AD','AD','AD','AD','AD','AD','AD','AD','AD','AD',
                   'AD','AD','AD','AD',
                   'Norm','Norm','Norm','Norm','Norm','Norm','Norm','Norm',
                   'Norm','Norm','Norm','Norm','Norm','Norm','Norm','Norm',
                   'Norm','Norm','Norm','Norm','Norm','Norm')
group_colors <- c("AD" = "#FF9999", "Norm" = "#99CCFF")
ggplot(plotdat, aes(x = PC1, y = PC2, colour = group, group = group)) + 
  geom_point(size = 6) +
  scale_color_manual(values = group_colors) +
  theme_bw() + 
  theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) +
  stat_ellipse(type = "t", linetype = 2) +
  geom_hline(aes(yintercept = 0), colour = "gray88", linetype = "dashed") + 
  geom_vline(aes(xintercept = 0), colour = "gray88", linetype = "dashed") 
#Different_analysis_by_limma#
library(limma)
design <- model.matrix(~0 + factor(group4))
colnames(design) = levels(factor(group4))
rownames(design) = colnames(ADdata_group_filtered)
fit <- lmFit(ADdata_group_filtered, design)
cont.matrix <- makeContrasts(ADvsNorm = AD- Norm,levels = design)
fit1 <- contrasts.fit(fit, cont.matrix)
fit1 <- eBayes(fit1)
#get result#
result <- topTable(fit1, coef = 'ADvsNorm',adjust="fdr", number=Inf)
DEG_limma_voom <- na.omit(result)
write.table(DEG_limma_voom,"AD_vs_N_result(RNA_bulk).csv",row.names=TRUE,col.names=TRUE,sep=",")

####for_OA####
#load_data#
OAdata <- read_excel("OA_RNA_bulk_data.xlsx",col_names = TRUE)
#pre-differential_analysis_conduct#
#Data filtering#
OAdata_only <- OAdata%>%select(-GeneSymbol) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.))))
keep <- rowSums(OAdata_only >= 10, na.rm = TRUE) >= 2
filtered_df <- OAdata[keep, , drop = FALSE]
OAdata_filtered <- as.data.frame(filtered_df) 
rownames(OAdata_filtered ) <- OAdata_filtered $GeneSymbol
OAdata_filtered <- OAdata_filtered [ , -c(1)]
#grouping#
expr_OA <- OAdata_filtered [, 1:10]
expr_Norm <- OAdata_filtered [, 11:18]
ADdata_group_filtered <- cbind(expr_OA, expr_Norm)  
group3 <- c(rep('OA', ncol(expr_OA)), rep('Norm', ncol(expr_Norm)))
group3 <- factor(group3, levels = c("OA",'Norm'))
#normalize_by_DESeq2#
colData <- data.frame(row.names = colnames(ADdata_group_filtered),group = group3)
ADdata_donducted<-apply(ADdata_group_filtered,2, as.integer)
rownames(ADdata_donducted)<- rownames(ADdata_group_filtered)
dds <- DESeqDataSetFromMatrix(countData = ADdata_donducted, 
                              colData = colData,        
                              design = ~ group)    
#got_vstdata_Visualization#
dds <- DESeq(dds)
vst <- vst(dds) 
vst_data<-assay(vst)
#*#BOXPOT_check_data#*#
any(vst_data<0)
group_colors <- c("#FF9999", "#99CCFF")
colors <- rep(group_colors, times = c(10, 8))
par(cex = 0.7)
if(ncol(vst_data) > 40) par(cex = 0.5)
boxplot(as.matrix(vst_data), las = 2, col = colors)
#*#PCA#*#
after1 <- t(vst_data)
pr.out1<-prcomp(after1)
pcadat1<-pr.out1
plotdat1 <- as.data.frame(pcadat1$x[,1:3])
plotdat1$group <- c('OA','OA','OA','OA','OA','OA','OA','OA','OA','OA',
                   'Normal','Normal','Normal','Normal','Normal','Normal','Normal','Normal')
group_colors <- c("OA" = "#FF9999", "Normal" = "#99CCFF")
ggplot(plotdat1, aes(x = PC1, y = PC2, colour = group, group = group)) + 
  geom_point(size = 6) +
  scale_color_manual(values = group_colors) +
  theme_bw() + 
  theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_blank()) +
  stat_ellipse(type = "t", linetype = 2) +
  geom_hline(aes(yintercept = 0), colour = "gray88", linetype = "dashed") + 
  geom_vline(aes(xintercept = 0), colour = "gray88", linetype = "dashed") +
  geom_text_repel(aes(label = row.names(plotdat1)), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50')
#Different_analysis_by_DESeq2#
resultsNames(dds)
res1 <- results(dds, contrast = c("group", levels(group3)[1], levels(group3)[2]))
resOrdered1 <- res1[order(res1$padj), ]
DEG1 <- as.data.frame(resOrdered1)
write.table(DEG1,"OA_vs_N_result(RNA_bulk).csv",row.names=TRUE,col.names=TRUE,sep=",")

#####****hypergeometric_test****#####
OA_DEGs <- data.frame(genes = rownames(DEG1), DEG1, row.names = NULL)
AD_DEGs <- data.frame(genes = rownames(DEG_limma_voom), DEG_limma_voom, row.names = NULL)
AD_DEGs <- AD_DEGs %>%
  mutate(tendency = case_when(
    adj.P.Val < 0.05 & logFC > 1  ~ "UP",
    adj.P.Val < 0.05 & logFC < -1 ~ "DOWN",
    TRUE ~ NA_character_
  ))

OA_DEGs <- OA_DEGs %>%
  mutate(tendency = case_when(
    padj < 0.05 & log2FoldChange > 1  ~ "UP",
    padj < 0.05 & log2FoldChange < -1 ~ "DOWN",
    TRUE ~ NA_character_
  ))
library(DescTools)
U <- intersect(AD_DEGs$genes, OA_DEGs$genes)
AD_UP   <- AD_DEGs %>% filter(genes %in% U, tendency == "UP")   %>% pull(genes) %>% unique()
AD_DOWN <- AD_DEGs %>% filter(genes %in% U, tendency == "DOWN") %>% pull(genes) %>% unique()
OA_UP   <- OA_DEGs %>% filter(genes %in% U, tendency == "UP")   %>% pull(genes) %>% unique()
OA_DOWN <- OA_DEGs %>% filter(genes %in% U, tendency == "DOWN") %>% pull(genes) %>% unique()
#for_up#
N  <- length(U)
n1 <- length(AD_UP)
n2 <- length(OA_UP)
k  <- length(intersect(AD_UP, OA_UP))
p_hyper_up <- phyper(k - 1, n1, N - n1, n2, lower.tail = FALSE)
p_hyper_up
#for_down#
n1_down <- length(AD_DOWN)
n2_down <- length(OA_DOWN)
k_down  <- length(intersect(AD_DOWN, OA_DOWN))
p_hyper_down <- phyper(k_down - 1, n1_down, N - n1_down, n2_down, lower.tail = FALSE)
p_hyper_down
#####****effect_size_concordance_analysis(spearman)****#####
dfU <- AD_DEGs %>%
  filter(genes %in% U) %>%
  select(genes, lfc1 = logFC) %>%
  inner_join(
    OA_DEGs %>% select(genes, lfc2 = log2FoldChange),
    by = "genes"
  )
spearman <- cor.test(dfU$lfc1, dfU$lfc2, method = "spearman")
spearman
library(ggplot2)
ggplot(dfU, aes(x = lfc1, y = lfc2)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal() +
  labs(
    x = "AD logFC",
    y = "OA log2FoldChange",
    title = "Effect size concordance (Spearman correlation)"
  )
same <- sign(dfU$lfc1) == sign(dfU$lfc2)
mean(same, na.rm = TRUE)