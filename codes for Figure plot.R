library(readxl)
library(dplyr)
library(tidyverse)

proteins_raw <- openxlsx::read.xlsx('./proteomics/protein.xlsx') %>% rowid_to_column()

metadata <- openxlsx::read.xlsx('./metadata.xlsx')

colnames(proteins_raw)[8:15] <- metadata$prot_sampleID %>% na.exclude()

proteins_nor <- proteins_raw[,8:15]

# 数据过滤
kept.inx <- apply(proteins_nor, MARGIN = 1, function(x) {
  sum(is.na(x)) < 2  
}) 

proteins_filed <- proteins_nor[kept.inx, ]

# openxlsx::write.xlsx(proteins_filed,'./proteomics/proteins_filtered.xlsx')


# 缺失值填充

# 定义处理组列名（使用标准化后的列名）
group_cols <- list(
  group1 = colnames(proteins_filed)[1:4],
  group2 = colnames(proteins_filed)[5:8]
)

# 优化后的缺失值填充函数
fill_group_na <- function(df, cols) {
  # 提取处理组数据副本
  group_data <- df[, cols]
  
  # 计算行均值（使用原始数据）
  row_means <- rowMeans(group_data, na.rm = TRUE) %>% replace(., is.nan(.), 0.1)
  
  # 创建缺失值掩模
  na_mask <- is.na(group_data)
  
  # 批量替换缺失值
  group_data[na_mask] <- matrix(rep(row_means, ncol(group_data)), 
                                ncol = ncol(group_data))[na_mask]
  
  # 返回更新后的数据
  df[, cols] <- group_data
  return(df)
}

# 应用处理
for (group in names(group_cols)) {
  proteins_filed <- fill_group_na(proteins_filed, group_cols[[group]])
}

# 验证结果
cat("剩余NA数量:", sum(is.na(proteins_filed)), "\n")

# 数据归一化

# 计算每个样本的总丰度（列方向）
sample_sums <- colSums(proteins_filed, na.rm = TRUE)

# 找到总丰度最低的样本索引
min_sum_index <- which.min(sample_sums)
min_sum_value <- sample_sums[min_sum_index]

normalize_to_min_abundance <- function(data_matrix) {
  # 计算各样本总丰度
  sums <- colSums(data_matrix, na.rm = TRUE)
  
  # 定位最小总丰度样本
  ref_sum <- min(sums)
  
  # 计算缩放因子：min_sum / current_sum
  scaling_factors <- ref_sum / sums
  
  # 应用归一化（行方向按样本缩放）
  normalized_data <- sweep(data_matrix, 2, scaling_factors, FUN = "*")
  
  return(normalized_data)
}

# 应用归一化函数
proteins_normalized <- normalize_to_min_abundance(proteins_filed)

table(is.na(proteins_normalized))

# 验证结果（所有样本总丰度应等于最小值）
colSums(proteins_normalized, na.rm = TRUE)  # 输出应全为 min_sum_value


proteins_normalized <- proteins_normalized %>% rowid_to_column()

proteins_normalized <- right_join(proteins_raw[,1:7],proteins_normalized,by='rowid') 

rownames(proteins_normalized) <- paste0('protein',proteins_normalized$rowid)

proteins_normalized$rowid <- paste0('protein',proteins_normalized$rowid)

# PcoA 以及adonis分析
library(vegan)
library(tidyverse)
library(colorspace)

# 计算Bray-Curtis距离
bray_dist <- vegdist(t(proteins_normalized[,8:15]), method = "bray")

# 执行PCoA分析（使用cmdscale函数）
pcoa_result <- cmdscale(bray_dist, k = 3, eig = TRUE)

# 提取主坐标轴及解释率
pcoa_points <- as.data.frame(pcoa_result$points)
colnames(pcoa_points) <- c("PC1", "PC2", "PC3")
eig_ratio <- round(pcoa_result$eig / sum(pcoa_result$eig) * 100, 1)  # 各轴解释率

# 合并分组信息
pcoa_points$group <- c(rep('HS',4),rep("PF",4))

# 计算每个分组的凸包顶点（基于PC1和PC2）
hull_data <- pcoa_points %>%
  group_by(group) %>%          # 按分组变量（如group）分组
  slice(chull(PC1, PC2))       # 计算凸包顶点索引

# 绘制PCoA图（凸包代替椭圆）
ggplot(pcoa_points, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size=4) +
  geom_polygon(                # 添加凸包多边形
    data = hull_data, 
    aes(fill = group, color = group), 
    alpha = 0.1,               # 多边形透明度
    linetype = "dashed",       # 边界线类型（可选）
    linewidth = 0.5            # 边界线粗细（可选）
  ) +
  labs(x = paste0("PCoA1 (", eig_ratio[1], "%)"),
       y = paste0("PCoA2 (", eig_ratio[2], "%)")) +
  theme_bw() +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.8,0.2),
        panel.grid = element_blank(),
        axis.text = element_blank())+
  scale_color_manual(values = c('PF'='#66C2A5',
                                'HS'='#FC8D62'))+
  scale_fill_manual(values = c('PF'='#66C2A5',
                               'HS'='#FC8D62'))

ggsave('./Figures/protemics_pcoa.pdf', width = 4, height = 4)

# 使用adonis2函数进行分组差异检验
adonis_result <- adonis2(
  formula = bray_dist ~ pcoa_points$group, 
  permutations = 999,  # 置换次数
  method = "bray"
)

# 提取关键结果
r_squared <- round(adonis_result$R2[1], 2)
p_value <- adonis_result$`Pr(>F)`[1]

# 输出结果
cat("Adonis R²:", r_squared, "\np-value:", p_value)




# deseq2分析差异蛋白

library(DESeq2)

countData <- proteins_normalized[,8:15] %>% round()

condition <- factor(c(rep("HS",4),rep("PF",4)),levels = c("PF","HS"))
# 定义condition
condition
colData <- data.frame(row.names=colnames(countData), condition) # 样品信息矩阵
colData
#构建dds对象,开始DESeq流程
dds <- DESeqDataSetFromMatrix(countData, colData, design= ~ condition)
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
dds<- DESeq(dds)# 对dds进行Normalize
resultsNames(dds) # 查看结果的名称

res <- results(dds,name = "condition_HS_vs_PF") # 使用results()函数获取结果，并赋值给res
res <- na.exclude(as.data.frame(res))
res <- res[order(res$padj),]
resdata <- merge (as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)

table(resdata$padj <0.05 & resdata$log2FoldChange>0)
table(resdata$padj <0.05 & resdata$log2FoldChange<0)

resdata <- right_join(proteins_normalized[,1:7], resdata, by=c('rowid'='Row.names'))

openxlsx::write.xlsx(resdata, file="./proteomics/deseq2_results.xlsx",rowNames = F)

# 差异蛋白火山图

plot_data <- resdata

plot_data$color <- ifelse(plot_data$padj<0.05 & plot_data$log2FoldChange < 0, 'Down',
                          ifelse(plot_data$padj<0.05 & plot_data$log2FoldChange > 0, 'Up', 'N.S.'
                                 ))

plot_data$symbol <- ifelse(plot_data$padj<0.05 & plot_data$log2FoldChange < -2.5, 'Down_out',
                           ifelse(plot_data$padj<0.05 & plot_data$log2FoldChange > 2.5, 'Up_out', 'N.S.'
                           ))

plot_data$log2FoldChange[plot_data$log2FoldChange>2.5] <- 2.5
plot_data$log2FoldChange[plot_data$log2FoldChange< -2.5] <- -2.5



ggplot(plot_data, aes(log2FoldChange, -log10(padj), color=color))+
  geom_point(aes(shape = symbol),size=2, fill="white")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none')+
  scale_color_manual(values = c('Down'='#66C2A5',
                                'Up'='#FC8D62',
                                'N.S.'='#e8e8e8'))+
  scale_shape_manual(
    values = c("Up_out" = 17, "N.S." = 1, "Down_out" = 6),  # pch=17（实心右三角）, pch=1（空心圆）, pch=6（倒三角）
    labels = c("Up_out", "N.S.", "Down_out")
  )

ggsave('./Figures/diff_proteins_volcanol.pdf',height = 4, width = 4)


# 蛋白质功能富集

metadata$SCC <- log10(metadata$`Somatic_Cells(thousands/mL)`/1000)+3

# 体细胞数结果
library(ggsignif)

ggplot(metadata, aes(group, SCC, color=group,fill = group))+
  geom_boxplot(outlier.shape=NA, size=2) +
  geom_jitter(shape=21, size=4, fill="white", width=0.1, height=0)+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background =element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  scale_color_manual(values = c('PF'=darken('#66C2A5',0.5),
                                'HS'=darken('#FC8D62',0.5)))+
  scale_fill_manual(values = c('PF'='#66C2A5',
                               'HS'='#FC8D62'))+
  geom_signif(comparisons = list(c('PF','HS')), test = 't.test',color='black',tip_length = 0)


ggsave('./Figures/Somatic_Cells.pdf',height = 4, width = 4)

aov1 <- aov(log10(`Somatic_Cells(thousands/mL)`)~day+group+day*group+(1|cowID2), metadata)

summary(aov1)


# 加载必要的库
library(tidyverse)
library(vegan)
library(ggsignif)
library(Mfuzz)


# metaphlan_hostfree
metaphlan <- read.delim2("./tables/metaphlan_table.txt",skip = 1)

metaphlan <- metaphlan %>% separate(col = 'clade_name',remove = F,
                                    into = c('kindom','phylum','class','order','family','genus','species','strains'),
                                    sep = "\\|")

#species

meta_species <- metaphlan[!is.na(metaphlan$species) & is.na(metaphlan$strains),]

meta_species <- meta_species[,-9] %>%
  mutate_at(vars(9:ncol(.)), as.numeric)

############alpha diversity analysis##################

shannon <- diversity(t(meta_species[,9:48]), index = 'shannon')
simpson <- diversity(t(meta_species[,9:48]), "simpson")
invsimpson <- diversity(t(meta_species[,9:48]), "inv")

alpha_diversity <- as.data.frame(cbind(shannon, simpson, invsimpson)) 
alpha_diversity <- cbind(alpha_diversity,metadata)

# shannon
ggplot(data = alpha_diversity, aes(x=group,y=simpson, fill=group,colour = group)) +
  lvplot::geom_lv(k=5, outlier.shape = NA)+
  geom_boxplot(outlier.shape=NA, coef=0, fill="#00000000", size=2,) +
  geom_jitter(shape=21, size=4, fill="white", width=0.1, height=0)+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background =element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.text = element_blank(),
        # axis.title.x = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  scale_color_manual(values = c('PF'=darken('#66C2A5',0.5),
                                'HS'=darken('#FC8D62',0.5)))+
  scale_fill_manual(values = c('PF'=darken('#66C2A5',0),
                               'HS'=darken('#FC8D62',0)))+
  facet_grid(~day)+
  geom_signif(comparisons = list(c('PF','HS')), color='black',tip_length = 0)


ggsave("./Figures/alpha_simpson_连续.pdf", width = 10, height = 4)


# 整体比较

ggplot(data = alpha_diversity, aes(x=group,y=shannon, fill=group,colour = group)) +
  lvplot::geom_lv(k=5, outlier.shape = NA)+
  geom_boxplot(outlier.shape=NA, coef=0, fill="#00000000", size=2,) +
  geom_jitter(shape=21, size=4, fill="white", width=0.1, height=0)+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background =element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.text = element_blank(),
        # axis.title.x = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  scale_color_manual(values = c('PF'=darken('#66C2A5',0.5),
                                'HS'=darken('#FC8D62',0.5)))+
  scale_fill_manual(values = c('PF'=darken('#66C2A5',0),
                               'HS'=darken('#FC8D62',0)))+
  geom_signif(comparisons = list(c('PF','HS')), color='black',tip_length = 0)


ggsave("./Figures/alpha_shannon_整体.pdf", width = 4, height = 4)


aov1 <- aov(simpson~day+group+day*group+(1|cowID2), alpha_diversity)

summary(aov1)

aov1 <- aov(shannon~day+group+day*group+(1|cowID2), alpha_diversity)

summary(aov1)


# species过滤
# low count filter
rmn_feat <- nrow(meta_species)
minLen <- 0.2 * (ncol(meta_species)-9)
kept.inx <- apply(meta_species[,c(9:48)], MARGIN = 1, function(x) {
  sum(x > 0.01) >= minLen 
}) 

meta_species_remain <- meta_species[kept.inx,]

openxlsx::write.xlsx(meta_species_remain,'./Figures/meta_species_remain.xlsx')

# 使用PHATE降维分析

# 安装并加载reticulate包
library(reticulate)
library(phateR)

Sys.setenv(RETICULATE_CONDA = "C:/Users/sheng/anaconda3/Scripts/conda.exe")

use_condaenv("r-phate", required = TRUE)
phate <- import("phate")

bc_dist <- vegan::vegdist(t(meta_species_remain[,-c(1:9)]), method = "bray")


# 安装并加载必要包
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("mixOmics")
library(mixOmics)
library(ggplot2)

# 1. 数据准备
# 假设meta_species_remain是物种丰度矩阵，前9列为分类信息
X <- t(meta_species_remain[, -c(1:8)])  # 转置为样本×特征矩阵
Y <- factor(rep(c("HS", "PF"), each = nrow(X)/2))  # 替换为实际分组变量

# 2. 数据预处理（对数转换+标准化）
X_log <- log1p(X)  # 处理零值
X_scaled <- scale(X_log, center = TRUE, scale = TRUE)

# 3. sPLS-DA参数调优（交叉验证）
set.seed(123)
list.keepX <- c(seq(10, 100, 10))  # 测试不同特征选择数量
tune.splsda <- tune.splsda(
  X = X_scaled, 
  Y = Y,
  ncomp = 5,  # 初始尝试5个成分
  validation = 'Mfold',
  folds = 5,
  dist = 'max.dist',
  measure = "BER",  # 平衡错误率
  test.keepX = list.keepX,
  nrepeat = 10,
  progressBar = TRUE
)

# 4. 提取最优参数
ncomp <- tune.splsda$choice.ncomp$ncomp  # 最佳成分数
select.keepX <- tune.splsda$choice.keepX[1:5]  # 各成分最优特征数

# 5. 构建最终sPLS-DA模型
final.splsda <- splsda(
  X = X_scaled,
  Y = Y,
  ncomp = 5,
  keepX = select.keepX
)

# 6. 模型评估
perf.splsda <- perf(final.splsda, 
                    validation = 'Mfold',
                    folds = 5, 
                    nrepeat = 10,
                    auc = TRUE)
plot(perf.splsda, col = color.mixo(1:3), sd = TRUE)  # 显示分类错误率

# 7. 可视化结果
# 样本分布图
plotIndiv(final.splsda, 
          comp = c(1,2),
          ind.names = FALSE,
          ellipse = TRUE,
          legend = TRUE,
          col = c('PF'='#66C2A5','HS'='#FC8D62'),
          centroid = F,
          title = 'none',
          style = 'ggplot2')


ggsave("./Figures/splsda.pdf", width = 7, height = 5)


################# 组间距离和组内距离之间的差异性分析 #################
# 检查样本数量是否足够
if (attr(bc_dist, "Size") < 2) {
  stop("距离矩阵需要至少 2 个样本以生成组合")
}

# 提取下三角部分和样本对名称
lower_tri <- bc_dist[lower.tri(bc_dist, diag = FALSE)]
sample_pairs <- utils::combn(attr(bc_dist, "Labels"), 2, simplify = FALSE)  # 使用距离矩阵的标签属性
pair_names <- sapply(sample_pairs, paste, collapse = "-")

# 创建结果数据框
distance_list <- data.frame(Sample_Pair = pair_names, Distance = lower_tri)
head(distance_list, 10)  # 预览前 10 个样本对

distance_list <- separate(distance_list, col = 'Sample_Pair', into = c('SampleID1','SampleID2'), sep = '-')

distance_list <- left_join(distance_list, metadata[,1:3], by=c('SampleID1'='SampleID'))
distance_list <- left_join(distance_list, metadata[,1:3], by=c('SampleID2'='SampleID'))


distance_list_zujian <- distance_list %>% subset(day.x == day.y & day_group.x != day_group.y)
distance_list_zujian$group <- 'zujian'


distance_list_zunei <- distance_list %>% subset(day.x == day.y & day_group.x == day_group.y)
distance_list_zunei$group <- 'zunei'

distance_list <- rbind(distance_list_zujian, distance_list_zunei)


ggplot(data = distance_list, aes(x=group,y=Distance , fill=group,colour = group)) +
  lvplot::geom_lv(k=5, outlier.shape = NA)+
  geom_boxplot(outlier.shape=NA, coef=0, fill="#00000000", size=2,) +
  geom_jitter(shape=21, size=4, fill="white", width=0.1, height=0)+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background =element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.text = element_blank(),
        # axis.title.x = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  scale_color_manual(values = c('zunei'=darken('#66C2A5',0.5),
                                'zujian'=darken('#FC8D62',0.5)))+
  scale_fill_manual(values = c('zunei'=darken('#66C2A5',0),
                               'zujian'=darken('#FC8D62',0)))+
  geom_signif(comparisons = list(c('zunei','zujian')), test = 't.test', color='black',tip_length = 0)+
  facet_grid(~day.x)


ggsave("./Figures/alpha_shannon_整体.pdf", width = 4, height = 4)



# 使用预先计算的距离矩阵，设置knn.dist.method为"precomputed"
phate_result <- phate(
  bc_dist,
  knn.dist.method = "precomputed",  # 指定使用预计算距离
  knn = 5,                          # 近邻数，根据数据调整
  ndim = 3,                         # 输出维度，通常2或3
  t = "auto",                       # 自动选择扩散时间
  n.landmark = 2000,                # 使用地标点加速计算，可选
  gamma = 1,                        # 熵正则化参数
  verbose = TRUE                    # 显示进度信息
)

# 获取PHATE坐标
phate_coords <- phate_result$embedding %>% as.data.frame() %>% rownames_to_column(var = 'SampleID')

metadata <- unite(metadata, col = day_group, c(2,3), remove = F)

phate_coords <- left_join(phate_coords, metadata,by=c('SampleID'='SampleID'))

# 可视化
ggplot(phate_coords, aes(PHATE1, PHATE2, color = group)) +
  theme_bw()+
  theme(axis.title = element_text(size = 10),
        axis.title.y = element_text(angle = 90),
        # axis.text = element_blank(),
        legend.position = 'none',
        legend.title = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_color_manual(values = c('PF'=darken('#66C2A5',0),
                                'HS'=darken('#FC8D62',0)))+
  scale_fill_manual(values = c('PF'=darken('#66C2A5',-0),
                                'HS'=darken('#FC8D62',-0)))+
  stat_ellipse(
    geom = "polygon", 
    level = 0.9,      # 置信水平（可调为0.9或其他）
    alpha = 0.3,       # 填充透明度
    aes(fill = group), # 按分组填充颜色
    show.legend = FALSE)+
  geom_point(size=4,shape=21,aes(fill = group), color='black') +
  facet_grid(~day)

ggsave('./Figures/PHATE降维_动态.pdf',width = 12, height = 4)

bc_dist_df <- as.data.frame(bc_dist)

for (i in seq(1,9,2)) {
  subsamples <- metadata %>% filter(day == i) %>% pull(SampleID)
  
  submetadata <- metadata[metadata$day==i,]
  
  sub_bc_dist <- vegan::vegdist(t(meta_species[,subsamples]), method = "bray")
  
  # 基于BC距离计算HS和PF的差异性
  adonis_result <-  adonis2(sub_bc_dist ~ group, data = submetadata)
  
  print(paste0("###### Day", i," ###############################"))
  print(adonis_result)
  print("#####################################")
}

adonis2(bc_dist ~ day+group+day*group+(1|cowID2), data = metadata)

# 计算PHATE2维度的差异性

ggplot(phate_coords, aes(day, PHATE1, color = group)) +
  geom_smooth(aes(fill = group))+
  geom_point(shape=21, size=4, fill="white") +
  theme_bw()+
  theme(
    # axis.text = element_blank(),
    legend.position = c(0.1,0.8),
    legend.title = element_blank(),
    axis.ticks.length = unit(2,'mm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  scale_color_manual(values = c('PF'=darken('#66C2A5',0.2),
                                'HS'=darken('#FC8D62',0.2)))+
  scale_fill_manual(values = c('PF'=darken('#66C2A5',-0.5),
                               'HS'=darken('#FC8D62',-0.5)))+
  scale_x_continuous(limits=c(1, 9), breaks=c(1,3,5,7,9), expand = c(0.02, 0.02))

ggsave('./Figures/PHATE1_动态.pdf',width = 5, height = 4)


# 分析获得p值
ggplot(phate_coords, aes(group, PHATE1, color = group)) +
  geom_boxplot()+
  geom_point(shape=21, size=4, fill="white") +
  theme_bw()+
  theme(
    axis.text = element_blank(),
    legend.position = c(0.1,0.8),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  scale_color_manual(values = c('PF'=darken('#66C2A5',0.2),
                                'HS'=darken('#FC8D62',0.2)))+
  scale_fill_manual(values = c('PF'=darken('#66C2A5',-0.5),
                               'HS'=darken('#FC8D62',-0.5)))+
  facet_grid(~day)+
  geom_signif(comparisons = list(c('PF','HS')))

ggsave('./Figures/PHATE1_动态.pdf',width = 6, height = 4)


# 计算不同时间HS和PF之间的BC距离

bc_dist_df <- as.matrix(bc_dist) %>% as.data.frame() %>% rownames_to_column(var = 'SampleID1')

bc_dist_df_long <- pivot_longer(bc_dist_df, cols = 2:40, values_to = 'bc_dist', names_to = 'SampleID') %>% unite(col = 'samplepair',1:2,remove = F)

bc_dist_df_long <- bc_dist_df_long[!duplicated(bc_dist_df_long$samplepair) & bc_dist_df_long$bc_dist != 0, ]

bc_dist_df_long <- left_join(bc_dist_df_long, metadata,by='SampleID')

bc_dist_df_long <- left_join(bc_dist_df_long, metadata,by=c('SampleID1'='SampleID'))

bc_dist_df_long <- bc_dist_df_long[bc_dist_df_long$group.x != bc_dist_df_long$group.y, ]

bc_dist_df_long <- bc_dist_df_long[bc_dist_df_long$day.x == bc_dist_df_long$day.y, ]


ggplot(bc_dist_df_long, aes(x=day.x,y=bc_dist))+
  lvplot::geom_lv(k=5, outlier.shape = NA, aes(x=day.x,y=bc_dist, color=day.x, fill=day.x, group = day.x))+
  geom_boxplot(outlier.shape=NA, coef=0, fill="#00000000", size=2,aes(x=day.x,y=bc_dist, color=day.x, fill=day.x, group = day.x)) +
  geom_jitter(shape=21, size=4, fill="white", width=0.1, height=0,aes(x=day.x,y=bc_dist, color=day.x, fill=day.x, group = day.x))+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background =element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        # axis.text = element_blank(),
        # axis.title.x = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  scale_color_gradient(low = darken('#66C2A5',0.5), high = darken('#FC8D62',0.5))+
  scale_fill_gradient(low = darken('#66C2A5',0), high = darken('#FC8D62',0))+
  scale_x_continuous(limits=c(0, 10), breaks=c(1,3,5,7,9))

ggsave('./Figures/组间BC动态变化.pdf',width = 7, height = 4)


######时间序列分析############

library(Mfuzz)
library(maSigPro)

library(tidyverse)

# 定义通用分析函数
analyze_group <- function(group_name, target_days) {
  # 筛选组别元数据
  metadata_group <- metadata %>%
    filter(group == group | day == 1) %>%
    mutate(day_group = ifelse(day == 1, "baseline", "target"))
  
  # 确保数据格式正确
  meta_species_group <- meta_species_remain %>%
    as.data.frame() %>%
    remove_rownames() %>%
    column_to_rownames("species") %>%
    dplyr::select(all_of(metadata_group$SampleID))
  
  # 定义结果存储列表
  result_list <- list()
  
  # 循环处理每个目标天数（重命名循环变量避免冲突）
  for (target_day in target_days) {  # 关键修改点：day -> target_day
    # 筛选目标天数数据
    metadata_temp <- metadata_group %>%
      filter(day %in% c(1, target_day)) %>%  # 明确使用target_day
      mutate(
        day = factor(
          day, 
          levels = c(1, target_day),  # 明确指定两个水平
          labels = c("baseline", "target")  # 添加可读标签
        )
      )
    
    # 确保样本存在
    valid_samples <- intersect(colnames(meta_species_group), metadata_temp$SampleID)
    species_temp <- meta_species_group[, valid_samples, drop = FALSE]
    
    # 计算统计指标
    stats_df <- species_temp %>%
      tibble::rownames_to_column("species") %>%
      pivot_longer(
        cols = -species,
        names_to = "SampleID",
        values_to = "abundance"
      ) %>%
      left_join(metadata_temp, by = "SampleID") %>%
      group_by(species) %>%
      summarise(
        wilcox_p = wilcox.test(
          abundance ~ day,
          exact = FALSE  # 避免小样本警告
        )$p.value,
        mean_baseline = mean(abundance[day == "baseline"], na.rm = TRUE) + 1e-6,
        mean_target = mean(abundance[day == "target"], na.rm = TRUE) + 1e-6,
        log2fc = log2(mean_target / mean_baseline)
      ) %>%
      rename_with(~ paste0(., "_", target_day, "d"), -species)
    
    result_list[[as.character(target_day)]] <- stats_df
  }
  
  # 合并所有天数结果
  purrr::reduce(result_list, full_join, by = "species") %>%
    mutate(group = group_name) %>%
    relocate(group, species)
}

# 执行分析
target_days <- c(3, 5, 7, 9)
hs_results <- analyze_group("HS", target_days)
pf_results <- analyze_group("PF", target_days)

# 合并最终结果
final_results <- bind_rows(hs_results, pf_results) %>%
  arrange(group, species)

# 结果输出（示例）
final_results %>%
  dplyr::select(group, species, contains("wilcox"), contains("log2fc")) %>%
  head()

write.table(hs_results, './Figures/hs_results-1122.txt',sep = '\t',row.names = F)
write.table(pf_results, './Figures/pf_results-1122.txt',sep = '\t',row.names = F)



# 读取数据（示例数据为表达矩阵，行为基因，列为时间点）
data_matrix <- hs_results %>% dplyr::select(species,contains('mean_baseline_3d')|contains('mean_target')) 

data_matrix <- data_matrix %>% column_to_rownames(var = 'species')

# HS样本时间序列分析

eset_HS <- new("ExpressionSet", exprs = as.matrix(data_matrix))

# 数据预处理
eset_std_HS <- standardise(eset_HS)  # 标准化

# 确定最佳聚类数（以 6 类为例）
m_HS <- mestimate(eset_std_HS)  # 计算模糊参数

# 最佳聚类数计算

set.seed(18)

# tmp <- cselection(eset_std_HS, m=m_HS, crange=seq(5,20,2), repeats=5, visu=TRUE)

tmp <- Dmin(eset_std_HS, m=m_HS, crange=seq(2,20,1), repeats=5, visu=TRUE)

cl_HS <- mfuzz(eset_std_HS, c = 8, m = m_HS)  # 执行模糊 C 均值聚类

# 可视化聚类结果
# pdf("mfuzz_clusters.pdf", width = 8, height = 12)
# 
# mfuzz.plot2(
#   eset_std_HS, cl_HS,
#   mfrow = c(1, 1),  # 分面布局
#   time.labels = colnames(data_matrix),  # 时间标签
#   centre = TRUE,  # 显示聚类中心
#   centre.lwd = 2,
#   cluster
# )
# 
# dev.off()

# 计算Pearson相关系数及p值
library(Hmisc)
library(pheatmap)

cor_matrix <- rcorr(as.matrix(t(centroids)), type = "spearman")
p_values <- cor_matrix$P  # 获取p值矩阵

# 生成显著性星号矩阵（不显示相关系数）
sig_stars <- matrix("", nrow = nrow(p_values), ncol = ncol(p_values))
rownames(sig_stars) <- rownames(p_values)
colnames(sig_stars) <- colnames(p_values)

# 根据p值填充星号
sig_stars[p_values < 0.001] <- "***"
sig_stars[p_values >= 0.001 & p_values < 0.01] <- "**"
sig_stars[p_values >= 0.01 & p_values < 0.05] <- "*"

# 绘制热图（仅显示显著性星号）
pheatmap(
  mat = cor_matrix$r, # 创建全零矩阵作为背景
  display_numbers = sig_stars,  # 仅显示星号
  cluster_rows = TRUE,          # 行聚类
  cluster_cols = TRUE,          # 列聚类
  border_color = "NA",      # 单元格边框颜色
  fontsize_number = 12,         # 星号大小
  number_color = "black",       # 星号颜色
)

# 提取分析过程中标准化后的表达值（绘制曲线图用的那个值，而非原始蛋白表达值）
gene_cluster <- cl_HS$cluster
gene_standard <- eset_std_HS@assayData$exprs
gene_standard_cluster <- cbind(gene_standard[names(gene_cluster), ], gene_cluster)
head(gene_standard_cluster)

gene_standard_cluster <- gene_standard_cluster %>% as.data.frame()

colnames(gene_standard_cluster)[1:5] <- seq(1,9,2)

# 将新的簇进行替换
gene_standard_cluster$merged_clusters <- ifelse(gene_standard_cluster$gene_cluster %in% c(7), 1,
                                                ifelse(gene_standard_cluster$gene_cluster %in% c(1,2), 2,
                                                       ifelse(gene_standard_cluster$gene_cluster %in% c(6,8), 3,
                                                              ifelse(gene_standard_cluster$gene_cluster %in% c(3,4), 5,4)
                                                       )
                                                )
)

# openxlsx::write.xlsx(gene_standard_cluster, "./Figures/mfuzz_clusters_HS.xlsx")

# gene_standard_cluster <- openxlsx::read.xlsx('./Figures/mfuzz_clusters_HS.xlsx')

gene_standard_cluster <- gene_standard_cluster %>% rownames_to_column(var = 'species')

gene_standard_cluster_long <- pivot_longer(gene_standard_cluster, cols = 2:6, names_to = 'time', values_to = 'abun_change')

gene_standard_cluster_long$time <- as.numeric(gene_standard_cluster_long$time)

# ggplot可视化

ggplot(gene_standard_cluster_long, aes(x=time,y=abun_change))+
  geom_line(aes(group = species),alpha=0.5)+
  geom_smooth(method = 'loess',color='#B20000',fill='#C18484')+
  facet_wrap(~merged_clusters,scales = "free", nrow=2,ncol=3)+
  ylim(-2,2)+
  theme_bw()+
  theme(axis.title = element_blank(),
        # axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )+
  scale_x_continuous(limits=c(0.8, 9.2), breaks=c(1,3,5,7,9))

ggsave('./Figures/HS_clusters.pdf', width = 9, height = 6.2)  


means_pf <- pf_results %>% dplyr::select(species,contains('mean_baseline_3d')|contains('mean_target')) 

# 假设原始数据框名为df，物种列为species
# 提取数值列（排除物种名称列）
numeric_cols <- means_pf[, -1]  

# 按行标准化（均值为0，标准差为1）
normalized_data <- t(apply(numeric_cols, 1, function(x) {
  (x - mean(x)) / sd(x)
}))

# 合并物种名称和标准化后的数据
normalized_pf <- data.frame(
  species = means_pf$species,
  normalized_data,
  check.names = FALSE
)

# 查看结果
head(normalized_pf)

normalized_pf <- left_join(normalized_pf,gene_standard_cluster[,c(1,8)],by='species')

colnames(normalized_pf)[2:6] <- seq(1,9,2)

normalized_pf_long <- normalized_pf %>% pivot_longer(cols = 2:6, names_to = 'time', values_to = 'abun_change')

normalized_pf_long$time <- as.numeric(normalized_pf_long$time)

# ggplot可视化
ggplot(normalized_pf_long, aes(x=time,y=abun_change))+
  geom_line(aes(group = species),alpha=0.5)+
  geom_smooth(method = 'loess',color='#B20000',fill='#C18484')+
  facet_wrap(~merged_clusters,scales = "free", nrow=2,ncol=3)+
  ylim(-2,2)+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )+
  scale_x_continuous(limits=c(0.8, 9.2), breaks=c(1,3,5,7,9))

ggsave('./Figures/PF_clusters.pdf', width = 9, height = 6.2)  


# HS heatmap

hs_heatmap_data <- gene_standard_cluster %>% column_to_rownames(var = 'species') %>% arrange(by=merged_clusters)

hs_heatmap_data$merged_clusters <- factor(hs_heatmap_data$merged_clusters, levels = c(1,2,3,4,5))

color1 <- colorRampPalette(
  rev(c("#FFFFFF", "#3399FF", "#003366"))
)(42)

color2 <- colorRampPalette(
  rev(c("#990000", "#FF3333", "#FFFFFF"))
)(58)

pheatmap::pheatmap(hs_heatmap_data[,1:5], 
                   color = c(color1,color2),
                   annotation_row = hs_heatmap_data[,6:7], 
                   cluster_cols = F, 
                   cluster_rows = F,
                   filename = "./Figures/HSheatmap.pdf",
                   width = 8,
                   height = 8,
                   border_color=NA
)

# PF heatmap

pf_heatmap_data <- normalized_pf %>% column_to_rownames(var = 'species') %>% arrange(by=merged_clusters)

pf_heatmap_data$merged_clusters <- factor(pf_heatmap_data$merged_clusters, levels = c(1,2,3,4,5))

color1 <- colorRampPalette(
  rev(c("#FFFFFF", "#3399FF", "#003366"))
)(42)

color2 <- colorRampPalette(
  rev(c("#990000", "#FF3333", "#FFFFFF"))
)(58)


pheatmap::pheatmap(pf_heatmap_data[,1:5], 
                   color = c(color1,color2),
                   annotation_row = hs_heatmap_data[,6:7], 
                   cluster_cols = F, 
                   cluster_rows = F,
                   filename = "./Figures/PFheatmap.pdf",
                   width = 8,
                   height = 8,
                   border_color=NA)




# 计算merged_cluster质心
# 计算质心（各列均值）
centroids_hs <- hs_heatmap_data %>%
  group_by(merged_clusters) %>%
  summarise(
    centroid_1 = mean(`1`),
    centroid_3 = mean(`3`),
    centroid_5 = mean(`5`),
    centroid_7 = mean(`7`),
    centroid_9 = mean(`9`)
  )

# 结果输出
print(centroids_hs)

# pf
centroids_pf <- pf_heatmap_data %>%
  group_by(merged_clusters) %>%
  summarise(
    centroid_1 = mean(`1`),
    centroid_3 = mean(`3`),
    centroid_5 = mean(`5`),
    centroid_7 = mean(`7`),
    centroid_9 = mean(`9`)
  )

# 结果输出
print(centroids_pf)

#计算相关性

library(openxlsx)
library(tidyverse)
library(pheatmap)
library(Hmisc)
library(corrplot)

centriods <- read.xlsx('clusters质心.xlsx', rowNames = T)

centroids_hs <- t(centroids_hs[,2:6])
colnames(centroids_hs) <- paste0('HS',seq(1,5))


centroids_pf <- t(centroids_pf[,2:6])
colnames(centroids_pf) <- paste0('PF',seq(1,5))


centriods <- cbind(centroids_hs,centroids_pf) %>% as.data.frame()

# 相关性

correlation <- data.frame(group=paste0('HS',seq(1,5,1),'_PF',seq(1,5,1)),
                          correlation=c(cor.test(centriods$PF1,centriods$HS1)$estimate,
                                        cor.test(centriods$PF2,centriods$HS2)$estimate,
                                        cor.test(centriods$PF3,centriods$HS3)$estimate,
                                        cor.test(centriods$PF4,centriods$HS4)$estimate,
                                        cor.test(centriods$PF5,centriods$HS5)$estimate),
                          pvalue=c(cor.test(centriods$PF1,centriods$HS1)$p.value,
                                   cor.test(centriods$PF2,centriods$HS2)$p.value,
                                   cor.test(centriods$PF3,centriods$HS3)$p.value,
                                   cor.test(centriods$PF4,centriods$HS4)$p.value,
                                   cor.test(centriods$PF5,centriods$HS5)$p.value)
)


library(scales)

# 原始数据
values <- correlation$correlation
abs_values <- abs(values)

# 颜色映射 ---------------------------------------------------------------
# 创建双色渐变（负值蓝→白，正值白→红）
blue_white <- colorRampPalette(c("#0018A8", "white"))  # 深蓝到白
white_red <- colorRampPalette(c("white", "#B30000"))    # 白到深红

# 分离负值和正值
neg_mask <- values < 0
pos_mask <- values >= 0

# 处理负值颜色（标准化到0-1区间）
neg_values <- values[neg_mask]
if(length(neg_values) > 0){
  neg_norm <- rescale(neg_values, to = c(0, 1), from = c(min(values), 0))
  neg_colors <- blue_white(100)[as.integer(neg_norm * 99) + 1]
} else {
  neg_colors <- character(0)
}

# 处理正值颜色（标准化到0-1区间）
pos_values <- values[pos_mask]
if(length(pos_values) > 0){
  pos_norm <- rescale(pos_values, to = c(0, 1), from = c(0, max(values)))
  pos_colors <- white_red(100)[as.integer(pos_norm * 99) + 1]
} else {
  pos_colors <- character(0)
}

# 合并颜色
colors <- character(length(values))
colors[neg_mask] <- neg_colors
colors[pos_mask] <- pos_colors

# 线宽映射 ---------------------------------------------------------------
line_widths <- rescale(abs_values, to = c(1, 6))  # 线宽范围1到6

# 可视化呈现 -------------------------------------------------------------
par(mar = c(4, 4, 3, 6), bg = "gray95")
plot(NA, xlim = c(0.5, length(values)+0.5), ylim = c(-1, 1),
     xlab = "Index", ylab = "Value", main = "双维度数据映射",
     axes = FALSE, frame.plot = FALSE)

# 绘制带样式的线段
segments(x0 = 1:5-0.4, y0 = 0, y1 = values, 
         col = colors, lwd = line_widths*2, 
         lend = "round", lty = "solid")

# 添加参考线
abline(h = 0, col = "gray50", lty = 2)
axis(2, at = seq(-1, 1, 0.5), las = 1)

# 添加数值标签
text(x = 1:5, y = rep(-0.9, 5), labels = sprintf("%.2f", values),
     col = colors, font = 2, cex = 0.9)

# 添加双图例
legend(x = 5.8, y = 0.9, 
       legend = c("最小值(-0.18)", "0", "最大值(0.98)"),
       col = c("#0018A8", "white", "#B30000"), 
       pch = 15, pt.cex = 2.5, bty = "n",
       title = "颜色映射", xpd = TRUE)

legend(x = 5.8, y = 0.4, 
       legend = c("最小绝对值", "最大绝对值"),
       lwd = c(1, 6)*2, col = "gray30",
       title = "线宽映射", bty = "n", xpd = TRUE)

# HS和PF网络分析


hs_samples <- metadata %>% subset( group == 'HS')
pf_samples <- metadata %>% subset( group == 'PF')


select_species <- gene_standard_cluster[gene_standard_cluster$merged_clusters %in% c(3,5),]$species

# 提取cluster 1,2,5的菌种

network_species_hs <- meta_species_remain[,c('species',hs_samples$SampleID)] %>% subset(species %in% select_species) %>% remove_rownames() %>% column_to_rownames(var = 'species')

network_species_pf <- meta_species_remain[,c('species',pf_samples$SampleID)] %>% subset(species %in% select_species) %>% remove_rownames() %>% column_to_rownames(var = 'species')


# sparcc网络分析

# 环境准备
# 加载包
library(SpiecEasi)
library(parallel)
library(igraph)
library(tidyverse)
library(corrplot)
library(ggraph)

set.seed(123)
# 计算观测相关矩阵

spcorr <- cor(t(network_species_hs),method = 'spearman')
sp_pvalue <- cor.mtest(t(network_species_hs),method='spearman')$p

table(abs(spcorr)>0.7)
table(sp_pvalue<=0.5)

# 网络构建与筛选
# 设置阈值筛选
spcorr[lower.tri(spcorr,diag=TRUE)]=NA
sp_pvalue[lower.tri(sp_pvalue,diag=TRUE)]=NA

spcorr_df <- as.data.frame(as.table(spcorr)) %>% na.exclude()
sp_pvalue_df <- as.data.frame(as.table(sp_pvalue)) %>% na.exclude()

net_hs <- cbind(spcorr_df,sp_pvalue_df$Freq)

colnames(net_hs) <- c('species1','species2','corr','pvalue')

net_hs_filter <- net_hs %>% subset(abs(corr)>=0.7 & pvalue < 0.05)

net_hs[abs(net_hs$corr)>=0.7,]$species1
net_hs[abs(net_hs$pvalue)<=0.05,]$species1

net_hs_filter <- left_join(net_hs_filter,gene_standard_cluster[,c('merged_clusters','species')],by=c('species1'='species'))

write.table(net_hs_filter,'./Figures/network_source_target_hs.txt', sep = '\t',row.names = F)

species_cluster <- gene_standard_cluster[gene_standard_cluster$species %in% net_hs_filter$species1 | gene_standard_cluster$species %in% net_hs_filter$species2,7:8]

write.table(species_cluster,'./Figures/species_clusters-hs.txt', sep = '\t',row.names = F)

# 基础igraph可视化
net <- graph_from_data_frame(net_hs_filter[,1:2], directed=FALSE)

# 计算度和介
degree<-degree(net)
betweenness<-betweenness(net)
Node_nw_st <- data.frame(degree, betweenness)
Rank_stat <- rowMeans(cbind(rank(Node_nw_st[,1]), rank(Node_nw_st[,2])))
Node_nw_st <- cbind(Node_nw_st, Rank_stat)

write.table(Node_nw_st,file="./Figures/Node_nw_st-hs.txt", sep="\t",col.names = NA, quote=F)

#  导出Gephi格式
write_graph(net, "./Figures/network-hs.graphml", format="graphml")

meta_species_remain

# PF

set.seed(123)
# 计算观测相关矩阵

spcorr <- cor(t(network_species_pf),method = 'spearman')
sp_pvalue <- cor.mtest(t(network_species_pf),method='spearman')$p

table(spcorr>0.7)
table(sp_pvalue<=0.5)

# 网络构建与筛选
# 设置阈值筛选
spcorr[lower.tri(spcorr,diag=TRUE)]=NA
sp_pvalue[lower.tri(sp_pvalue,diag=TRUE)]=NA

spcorr_df <- as.data.frame(as.table(spcorr)) %>% na.exclude()
sp_pvalue_df <- as.data.frame(as.table(sp_pvalue)) %>% na.exclude()

net_pf <- cbind(spcorr_df,sp_pvalue_df$Freq)

colnames(net_pf) <- c('species1','species2','corr','pvalue')

net_pf_filter <- net_pf %>% subset(abs(corr)>=0.7 & pvalue <= 0.05)

net_pf[abs(net_pf$corr)>=0.7,]$species1
net_pf[abs(net_pf$pvalue)<=0.05,]$species1

write.table(net_pf_filter,'./Figures/network_source_target_PF.txt', sep = '\t',row.names = F)

species_cluster <- merged_clusters[merged_clusters$species %in% net_pf_filter$species1 | merged_clusters$species %in% net_pf_filter$species2,7:8]

write.table(species_cluster,'./Figures/species_clusters_PF.txt', sep = '\t',row.names = F)

# 基础igraph可视化
net <- graph_from_data_frame(net_pf_filter[,1:2], directed=FALSE)

# 计算度和介
degree<-degree(net)
betweenness<-betweenness(net)
Node_nw_st <- data.frame(degree, betweenness)
Rank_stat <- rowMeans(cbind(rank(Node_nw_st[,1]), rank(Node_nw_st[,2])))
Node_nw_st <- cbind(Node_nw_st, Rank_stat)

write.table(Node_nw_st,file="./Figures/Node_nw_st_pf.txt", sep="\t",col.names = NA, quote=F)

#  导出Gephi格式
write.graph(net, "./Figures/network-pf.graphml", format="graphml")


# cluster 1 2 5与VFDB基因的相关性

tpm_vfdb_top20 <- openxlsx::read.xlsx('./Figures/tpm_vfdb_top20.xlsx') %>% column_to_rownames(var = 'vfdbid') %>% t() %>% 
  as.data.frame() %>% rownames_to_column(var = 'meta_SampleID')

vfdb_species <- left_join(tpm_vfdb_top20,hs_species_t[,c(1,17:34)], by='meta_SampleID') %>% column_to_rownames(var = 'meta_SampleID')

vfdb_species_diff <- vfdb_species[,c(colnames(vfdb_species)[1:20], c('s__Klebsiella_pneumoniae','s__Enterococcus_faecalis','s__Citrobacter_braakii','s__Pseudomonas_putida','s__Serratia_marcescens',        
                                                                     's__Acinetobacter_parvus','s__Acinetobacter_bereziniae','s__Macrococcus_caseolyticus', 's__Pseudomonas_anguilliseptica'))] 



vfdb_species_corr <- cor(vfdb_species_diff, method = 'spearman')

vfdb_species_pvalue <- cor.mtest(vfdb_species_diff,method='spearman')$p

library(pheatmap)

pheatmap(vfdb_species_corr)

# 导出相关性数据，绘制cytoscape网络

vfdb_species_corr[lower.tri(vfdb_species_corr,diag=TRUE)]=NA
vfdb_species_pvalue[lower.tri(vfdb_species_pvalue,diag=TRUE)]=NA

vfdb_species_corr_df <- as.data.frame(as.table(vfdb_species_corr)) %>% na.exclude()
vfdb_species_pvalue_df <- as.data.frame(as.table(vfdb_species_pvalue)) %>% na.exclude()

vfdb_species_net <- cbind(vfdb_species_corr_df,vfdb_species_pvalue_df)

colnames(vfdb_species_net)[c(3,6)] <- c('Correlation','pvalue')

vfdb_species_net_filter <- vfdb_species_net %>% subset(pvalue<0.05&abs(Correlation)>0.5)

write.table(vfdb_species_net_filter,file="./Figures/vfdb_species_net_filter.txt", sep="\t",col.names = NA, quote=F)



vfdb_species_net_annot <- data.frame(identity=rownames(vfdb_species_corr),
                                     type=c(rep('vfdb',20),rep('species',9)))

# 基础igraph可视化
net <- graph_from_data_frame(vfdb_species_net_filter[,1:2], directed=FALSE)

# 计算度和介
degree<-degree(net)
betweenness<-betweenness(net)
Node_nw_st <- data.frame(degree, betweenness)
Rank_stat <- rowMeans(cbind(rank(Node_nw_st[,1]), rank(Node_nw_st[,2])))
Node_nw_st <- cbind(Node_nw_st, Rank_stat)
Node_nw_st <- rownames_to_column(Node_nw_st,var = 'identity')


vfdb_species_net_annot <- left_join(vfdb_species_net_annot,Node_nw_st,by='identity')


write.table(vfdb_species_net_annot,file="./Figures/vfdb_species_net_annot.txt", sep="\t",col.names = NA, quote=F)


################ 菌丰度差异分析##########################

meta_species_select <- meta_species_remain %>% subset(species %in% select_species & kindom != 'k__Eukaryota') %>% remove_rownames()

species_select <- meta_species_select[,8:48]

library(tidyverse)

# 定义更健壮的计算函数
calculate_level_sums <- function(data, level_col) {
  # 检查列是否存在
  if (!level_col %in% colnames(data)) {
    stop(paste("Column", level_col, "not found in data"))
  }
  
  data %>%
    group_by(!!sym(level_col)) %>%  # 使用sym+!!直接引用
    summarise(across(where(is.numeric), sum), .groups = "drop") %>%
    mutate(taxonomic_level = level_col) %>%
    rename_with(~"taxon_name", all_of(level_col))  # 更安全的rename方式
}

# 确保数据列名正确（检查实际列名）
print("Available columns:")
print(colnames(meta_species_select))

# 确认要计算的层级列确实存在
target_levels <- c("phylum", "class", "order", "family", "genus", "species")
valid_levels <- target_levels[target_levels %in% colnames(meta_species_select)]
print(paste("Valid levels to process:", paste(valid_levels, collapse = ", ")))

# 计算各层级的丰度总和（只处理存在的列）
level_sums <- map_dfr(
  .x = valid_levels,
  .f = ~ {
    print(paste("Processing level:", .x))  # 调试输出
    calculate_level_sums(meta_species_select, .x)
  }
)

# 合并结果
flat_table <- level_sums %>%
  select(taxonomic_level, taxon_name, everything()) %>%
  arrange(taxonomic_level, taxon_name)

# 查看结果
print(flat_table, n = 20, width = Inf)

# 保存结果
write_tsv(flat_table, "taxonomic_abundance_summary.tsv")


# 加载必要包
library(tidyverse)
library(broom)

# 假设数据已加载：
# species_select: 菌种丰度矩阵（行为菌种，列为样本）
# metadata: 包含SampleID和group(HS/PF)信息

# Step 1: 数据预处理
# 转换species_select为长格式（如果尚未转换）
species_long <- flat_table[,-1] %>%
  pivot_longer(
    cols = -taxon_name,
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>%
  left_join(metadata %>% select(SampleID, group), by = "SampleID") %>%
  filter(group %in% c("HS", "PF"))  # 确保只包含HS和PF组

# Step 2: 对每个菌种进行差异分析
results <- species_long %>%
  group_by(taxon_name) %>%
  summarise(
    # 计算组间均值
    mean_HS = mean(Abundance[group == "HS"], na.rm = TRUE),
    mean_PF = mean(Abundance[group == "PF"], na.rm = TRUE),
    # 计算log2FC（HS/PF），添加伪计数避免除零
    log2FC = log2((mean_HS + 1e-10) / (mean_PF + 1e-10)),
    # Wilcoxon检验
    wilcox_p = wilcox.test(Abundance ~ group, exact = FALSE)$p.value,
    # 计算显著性标记
    significance = case_when(
      wilcox_p < 0.001 ~ "***",
      wilcox_p < 0.01 ~ "**",
      wilcox_p < 0.05 ~ "*",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(wilcox_p)  # 按p值排序

# Step 3: 输出结果
print(results, n = Inf)  # 打印所有结果

# 可选：保存到CSV
write_csv(results, "./Figures/microbiome_HS_vs_PF_stats.csv")

results$average <- (results$mean_HS+results$mean_PF)/2

results$color <- ifelse(results$wilcox_p < 0.05 & results$log2FC > 0, 'HS',
                        ifelse(results$wilcox_p < 0.05 & results$log2FC < 0, 'PF', 'NS.'))


library(data.tree)
library(tidyverse)
library(ggtree)
library(ape)
library(phytools)



# ===== 第一步：读取表格 =====
taxonomy <- meta_species_select %>% select(kindom,phylum, class,order,family,genus,species)

# 构造边表
edges <- taxonomy %>%
  select(kindom,phylum, class,order,family,genus,species) %>%
  pivot_longer(cols = everything(), names_to = "level", values_to = "node") %>%
  group_by(row = rep(1:nrow(taxonomy), each = 7)) %>%
  mutate(next_node = lead(node)) %>%
  ungroup() %>%
  filter(!is.na(next_node)) %>%
  select(from = node, to = next_node) %>%
  distinct()

# 创建igraph对象
g <- graph_from_data_frame(edges, directed = TRUE)

# ==== 添加丰度信息到节点 ====

# 添加到图中
V(g)$Abundance <- results$average[match(V(g)$name, results$taxon_name)]

V(g)$color <- results$color[match(V(g)$name, results$taxon_name)]

# ==== ggraph 画图 ====
library(ggraph)
library(ggplot2)

ggraph(g, layout = 'dendrogram', circular = T) + 
  geom_edge_diagonal(color = "gray60", edge_width = 1) + 
  geom_node_point(aes(size =Abundance, fill = color, color= color), shape = 21) + 
  geom_node_text(aes(label = name), hjust = -0.1, size = 3) + 
  scale_size_continuous(range = c(3,10)) + 
  coord_flip() +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5),
        legend.position = 'none',
        panel.background = element_rect(fill = "white"))+
  scale_color_manual(values = c('PF'='black',
                                'HS'='black',
                                'NS.'=NA,
                                'NA'=NA))+
  scale_fill_manual(values = c('PF'=darken('#66C2A5',0),
                               'HS'=darken('#FC8D62',0),
                               'NS.'='grey',
                               'NA'='grey'))


ggsave('./Figures/different-species.pdf', width = 8, height = 8)


#################### 分析 肺炎克雷伯菌的具体时间变化 ##################

time_microbe <- merge(metadata_diff, t(meta_species_select[,8:48] %>% column_to_rownames(var = 'species')) , by='row.names')

time_microbe$day <- as.numeric(time_microbe$day)

ggplot(data = time_microbe, aes(x=day,y=s__Klebsiella_pneumoniae,colour = group)) +
  geom_smooth(aes(fill = group))+
  geom_point(shape=21, size=4, fill="white") +
  theme_bw()+
  theme(
    # axis.text = element_blank(),
    legend.position = c(0.1,0.8),
    legend.title = element_blank(),
    axis.ticks.length = unit(2,'mm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  scale_color_manual(values = c('PF'=darken('#66C2A5',0.2),
                                'HS'=darken('#FC8D62',0.2)))+
  scale_fill_manual(values = c('PF'=darken('#66C2A5',-0.5),
                               'HS'=darken('#FC8D62',-0.5)))

ggsave('./Figures/s__Klebsiella_pneumoniae_动态.pdf',width = 5, height = 4)


time_microbe$SCC <- log10(time_microbe$`Somatic_Cells(thousands/mL)`/1000)+3


lm(s__Klebsiella_pneumoniae~SCC, data = time_microbe[-5,]) %>% summary()


ggplot(data = time_microbe[-5,], aes(x=s__Klebsiella_pneumoniae,y=SCC)) +
  geom_smooth(method = 'lm')+
  geom_point(shape=21, size=4, fill="white") +
  theme_bw()+
  theme(
    # axis.text = element_blank(),
    legend.position = c(0.1,0.8),
    legend.title = element_blank(),
    axis.ticks.length = unit(2,'mm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

ggsave('./Figures/s__Klebsiella_pneumoniae_SCC_动态.pdf',width = 5, height = 4)


####### 肺炎克雷伯菌中VFDB基因拷贝数 #############

VFDB_matrix <- read_tsv("./tables/vf_gene_copy_matrix.tsv")

VFDB_matrix <- VFDB_matrix %>% column_to_rownames(var = 'VF_Gene')

# 将拷贝数大于10的值用10代替
VFDB_matrix_processed <- as.matrix(VFDB_matrix)
VFDB_matrix_processed[VFDB_matrix_processed > 3] <- 3

library(pheatmap)

# 绘制热图
# 设置高分辨率输出
png("./Figures/VFDB_heatmap.png", 
    width = 12, height = 8, 
    units = "in", res = 600)  # 600 dpi高分辨率

# 绘制热图
pheatmap(VFDB_matrix_processed,
         color = colorRampPalette(c("white", "red"))(10),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "none",
         show_colnames = FALSE,
         show_rownames = FALSE,
         border_color = NA,
         fontsize = 8, 
         legend = TRUE,
         display_numbers = FALSE,
         treeheight_row = 50,
         treeheight_col = 50,
         annotation_legend = FALSE)

dev.off()  # 关闭图形设备

## 提取VF基因


# 基因组总数
n_genomes <- ncol(VFDB_matrix)

# 计算每个基因在多少个基因组中 copy number > 1
counts <- apply(VFDB_matrix, 1, function(x) sum(x >= 1))

# 计算占比
ratio <- counts / n_genomes

# 筛选超过 80% 基因组 copy number > 1 的基因
selected <- data.frame(
  Gene = names(counts),
  NumGenomes = counts,
  Ratio = ratio
) %>% filter(Ratio > 0.5)

# 保存结果
write_tsv(selected, "vfdb_high_copy_genes.tsv")


##### 统计各个基因组中VFgene的copy数########

genome_sum <- read_tsv("./tables/genome_vf_summary.tsv")

# 计算均值，避免除零错误
mean_primary <- mean(genome_sum$Total_Copy_Count, na.rm = TRUE)
mean_secondary <- mean(genome_sum$Unique_VF_Genes, na.rm = TRUE)

# 计算缩放因子
scale_factor <- mean_primary / mean_secondary

# 创建缩放后的次Y轴数据，以便与主Y轴使用同一尺度
genome_sum$Unique_VF_Genes_scaled <- genome_sum$Unique_VF_Genes * scale_factor


# 创建ggplot对象，x轴为Assembly_Accession
ggplot(genome_sum, aes(x = reorder(Assembly_Accession, desc(Total_Copy_Count)))) +
  
  # 绘制主Y轴数据 (Total_Copy_Count)，例如用柱状图表示总量
  geom_segment(aes(x = reorder(Assembly_Accession, desc(Total_Copy_Count)),
                   y = 0, 
                   xend = reorder(Assembly_Accession, desc(Total_Copy_Count)), 
                   yend = Total_Copy_Count), alpha = 0.7, color='blue') + 
  
  # 绘制次Y轴数据 (缩放后的Unique_VF_Genes)，例如用折线或点表示趋势
  geom_point(aes(y = Unique_VF_Genes_scaled), color = "salmon", size = 1) +
  
  # 美化主题
  theme_classic() +
  theme(
    axis.text.x = element_blank(), # 旋转x轴标签避免重叠
    legend.position = "inside",
    legend.position.inside = c(0.8,0.8)
    # 图例放在底部
  )


ggsave('./Figures/genome_vfdb_summary.pdf', height = 6, width = 15)


##### VF基因差异分析 #############

vfdb_gene_count <- read_tsv("./tables/vfdb_abundance_matrix.tsv") %>% column_to_rownames(var = 'GeneID')

# 加载必要的包
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(apeglm)

# --------------------------
# 步骤1: 创建模拟数据（实际使用时替换为你的数据）
# --------------------------
set.seed(123)

metadata <- openxlsx::read.xlsx('./metadata.xlsx')

rownames(metadata) <- metadata$SampleID

metadata <- metadata[colnames(vfdb_gene_count),]

vfdb_gene_count <- vfdb_gene_count*1e9

# --------------------------
# 步骤2: 构建DESeq2对象并运行模型
# --------------------------
# 创建DESeqDataSet对象
dds <- DESeqDataSetFromMatrix(
  countData = vfdb_gene_count %>% round(),
  colData = metadata,
  design = ~ group + day + group:day  # 包含交互项的模型
)

# 查看group因子的水平（第一个为参考组）
levels(colData(dds)$group)

colData(dds)$group <- factor(colData(dds)$group, levels = c("PF", "HS"))

# 过滤低表达基因（至少5个样本的计数>10）
keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep, ]

# 运行DESeq2分析（使用似然比检验LRT检验交互效应）
dds <- DESeq(
  dds,
  test = "LRT",  # 似然比检验
  reduced = ~ group + day  # 简化模型（无交互项）
)

resultsNames(dds)


# 提取处理组主效应结果（Group B vs A）
res_group <- results(
  dds,
  contrast = c("group", "HS", "PF"),
  alpha = 0.05
)

res_group_view <- as.data.frame(res_group)


res_group_view$color <- ifelse(res_group_view$pvalue < 0.05 & res_group_view$log2FoldChange > 0, 'HS',
                        ifelse(res_group_view$pvalue < 0.05 & res_group_view$log2FoldChange < 0, 'PF', 'NS.'))

res_group_view$log2FoldChange[res_group_view$log2FoldChange>10] <- 10

res_group_view$log2FoldChange[res_group_view$log2FoldChange < -10] <- -10

ggplot(res_group_view,
       aes(x=log2FoldChange, y= -log10(pvalue), fill=color, color=color))+
  geom_point(shape=21, size=3) +
  theme_bw()+
  theme(legend.position = 'none',
        legend.background =element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        # axis.text = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  scale_color_manual(values = c('PF'='black',
                                'HS'='black',
                                'NS.'=NA,
                                'NA'=NA))+
  scale_fill_manual(values = c('PF'=darken('#66C2A5',0),
                               'HS'=darken('#FC8D62',0),
                               'NS.'='grey',
                               'NA'='grey'))+
  xlim(-10,10)


ggsave('./Figures/different-vfdb gene.pdf', width = 5, height = 5)


########### VFG048632(gb|WP_004176434) # 差异性时间分析 ##############

vfdb_gene_count_select <- t(vfdb_gene_count) %>% as.data.frame() %>% select(`VFG048632(gb|WP_004176434)`)

vfdb_gene_count_select <- merge(metadata, vfdb_gene_count_select, by='row.names')


# shannon
ggplot(data = vfdb_gene_count_select, aes(x=group,y=`VFG048632(gb|WP_004176434)`, fill=group,colour = group)) +
  lvplot::geom_lv(k=5, outlier.shape = NA)+
  geom_boxplot(outlier.shape=NA, coef=0, fill="#00000000", size=2,) +
  geom_jitter(shape=21, size=4, fill="white", width=0.1, height=0)+
  theme_bw()+
  theme(legend.position = 'none',
        legend.background =element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        axis.text = element_blank(),
        # axis.title.x = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  scale_color_manual(values = c('PF'=darken('#66C2A5',0.5),
                                'HS'=darken('#FC8D62',0.5)))+
  scale_fill_manual(values = c('PF'=darken('#66C2A5',0),
                               'HS'=darken('#FC8D62',0)))+
  facet_grid(~day)+
  geom_signif(comparisons = list(c('PF','HS')), color='black',tip_length = 0)


ggsave("./Figures/VFG048632_连续.pdf", width = 10, height = 4)

vfdb_gene_count_select$SCC <- log10(vfdb_gene_count_select$`Somatic_Cells(thousands/mL)`/1000)+3

lm(`VFG048632(gb|WP_004176434)`~SCC, data = vfdb_gene_count_select[-39,]) %>% summary()

ggplot(data = vfdb_gene_count_select[-39,], aes(x=`VFG048632(gb|WP_004176434)`,y=SCC)) +
  geom_smooth(method = 'lm',color='red',fill='salmon')+
  geom_point(shape=21, size=4, fill="white") +
  theme_bw()+
  theme(
    # axis.text = element_blank(),
    legend.position = c(0.1,0.8),
    legend.title = element_blank(),
    axis.ticks.length = unit(2,'mm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

ggsave('./Figures/VFG048632_SCC.pdf',width = 5, height = 4)



################## CAZyme 基因丰度 ####################


tpm <- read.csv("D:/BaiduSyncdisk/内师大/在研课题/国自然-热应激/昌平牛奶微生物组/tables/merged_tpm.csv")

apply(tpm[,-1], 2, sum)


# low count filter
rmn_feat <- nrow(tpm)
minLen <- 0.2 * (ncol(tpm)-1)
kept.inx <- apply(tpm[,-1], MARGIN = 1, function(x) {
  sum(x > 10) >= minLen 
}) 

tpm <- tpm[kept.inx,]

tpm <- tpm %>%
  mutate_at(vars(2:ncol(.)), as.numeric)

apply(tpm[,2:41], 2, sum)


library(edgeR)

library(edgeR)
library(dplyr)
library(tibble)

# 1. 创建 DGEList 并计算 TMM
dge <- tpm %>% remove_rownames() %>% 
  column_to_rownames("target_id") %>%
  DGEList() %>%
  calcNormFactors(method = "TMM")

# 2. 转 CPM
tpm_norm <- dge %>%
  cpm(log = FALSE) %>%
  as.data.frame() %>%
  rownames_to_column("target_id")


# 查看前5行
apply(tpm_norm[,2:41], 2, sum)

# 保存结果
write_tsv(tpm_norm, "normalized_tpm.tsv")


vfdb_annota <- openxlsx::read.xlsx("./tables/metagG.vfdb.annotation.unique.xlsx", colNames = F)

vfdb_annota <- vfdb_annota %>%
  mutate(
    temp = str_split_fixed(X2, pattern = "\\(gb\\|", n = 2),
    gene_id = temp[,1],
    accession = str_remove(temp[,2], "\\)$")  # 移除末尾的右括号
  ) %>%
  select(-temp)

tpm_nor <- left_join(tpm_norm, vfdb_annota[,c(1,13)], by=c('target_id'='X1'))

tpm_vfdb <- tpm_nor %>% na.exclude()


vfdb_name <- openxlsx::read.xlsx("D:/BaiduSyncdisk/内师大/在研课题/国自然-热应激/昌平牛奶微生物组/tables/VFDB_list.xlsx", colNames = F)


vfdb_name <- vfdb_name %>% 
  extract(
    col = X1,
    into = c("gene_id", "accession", "gene_name", "protein_function", "pathway", "strain"),
    regex = "^([^\\s]+)\\(gb\\|([^)]+)\\)\\s+\\(([^)]+)\\)\\s+([^\\[]+)\\[([^\\]]+)\\]\\s+\\[([^\\]]+)\\]$",
    remove = TRUE
  ) %>%
  mutate_all(str_trim)  # 去除多余空格


tpm_vfdb <- left_join(tpm_vfdb,vfdb_name,by='gene_id')


openxlsx::write.xlsx(tpm_vfdb, './Figures/vfdb_tpm.xlsx')



# Step 1: 数据预处理
# 转换species_select为长格式（如果尚未转换）
vfdb_long <- tpm_vfdb[,c(2:42)] %>%
  pivot_longer(
    cols = -gene_id,
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>%
  left_join(metadata %>% select(SampleID, group), by = "SampleID") %>%
  filter(group %in% c("HS", "PF"))  # 确保只包含HS和PF组

# Step 2: 对每个菌种进行差异分析
results_vfdb <- vfdb_long %>%
  group_by(gene_id) %>%
  summarise(
    # 计算组间均值
    mean_HS = mean(Abundance[group == "HS"], na.rm = TRUE),
    mean_PF = mean(Abundance[group == "PF"], na.rm = TRUE),
    # 计算log2FC（HS/PF），添加伪计数避免除零
    log2FC = log2((mean_HS + 1e-10) / (mean_PF + 1e-10)),
    # Wilcoxon检验
    wilcox_p = wilcox.test(Abundance ~ group, exact = FALSE)$p.value,
    # 计算显著性标记
    significance = case_when(
      wilcox_p < 0.001 ~ "***",
      wilcox_p < 0.01 ~ "**",
      wilcox_p < 0.05 ~ "*",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(wilcox_p)  # 按p值排序

results_vfdb <- left_join(results_vfdb,tpm_vfdb,by='gene_id')

# 可选：保存到CSV
write_csv(results_vfdb, "./Figures/vfdb_HS_vs_PF_stats.csv")

results_huamnn$log2FoldChange <- log2(results_huamnn$mean_HS/results_huamnn$mean_PF)
results_huamnn$color <- ifelse(results_huamnn$wilcox_p < 0.05 & results_huamnn$log2FC > 0, 'HS',
                               ifelse(results_huamnn$wilcox_p < 0.05 & results_huamnn$log2FC < 0, 'PF', 'NS.'))

ggplot(results_huamnn,
       aes(x=log2FoldChange, y= -log10(wilcox_p), fill=color, color=color))+
  geom_point(shape=21, size=3) +
  theme_bw()+
  theme(legend.position = 'none',
        legend.background =element_blank(),
        legend.title = element_blank(),
        legend.key = element_blank(),
        # axis.text = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  scale_color_manual(values = c('PF'='black',
                                'HS'='black',
                                'NS.'=NA,
                                'NA'=NA))+
  scale_fill_manual(values = c('PF'=darken('#66C2A5',0),
                               'HS'=darken('#FC8D62',0),
                               'NS.'='grey',
                               'NA'='grey'))+
  xlim(-5,5)


ggsave('./Figures/different-pathway.pdf', width = 5, height = 5)




# 分析disease通路的蛋白与HS cluster菌的相关性
library(readxl)
meta_species_remain <- read_excel("Figures/meta_species_remain.xlsx")
View(meta_species_remain)                                                                                                                                          

mfuzz_clusters_HS <- read_excel("Figures/mfuzz_clusters_HS.xlsx")
View(mfuzz_clusters_HS)   

hs_species <- meta_species_remain %>% subset(species %in% mfuzz_clusters_HS[mfuzz_clusters_HS$merged_clusters %in% c(1,2,5),]$species)

hs_species_t <- hs_species[,c(7,9:48)] %>% remove_rownames() %>% column_to_rownames(var = 'species') %>% t() %>% as.data.frame() %>% rownames_to_column(var = 'sampleid')

hs_species_t <- left_join(metadata, hs_species_t, by=c('meta_SampleID'='sampleid'))

openxlsx::write.xlsx(hs_species_t,'./Figures/目标菌丰度.xlsx')

# 蛋白数据整理













