##### 让报错变成英文 方便google
Sys.setenv(LANGUAGE = "en")
##### 禁止chr转成factor
options(stringsAsFactors = FALSE)
##### 清空环境
rm(list=ls())

getwd()
setwd("D:/bioinformation analysis/covid swab")

# data clean requirement:
# 1.Add the group column for each sample.
# 2.the column names of data are the samples
# 3.the row names of data are the symbols, sampleid and group

library(readxl) 

# 1.数据导入与预处理
# 1.1 从Excel文件导入原始数据
data_swab <- read_excel("01-data/swab_proteins_imarker_with_sampleID.xlsx")
data_plasma <- read_excel("01-data/plasma_proteins_imarker_with_sampleID.xlsx")
data_platelet <- read_excel("01-data/platelet_proteins_imarker_with_sampleID.xlsx")

# 1.2 数据聚合（按基因符号分组求均值）
# 按'Gene Symbol'列分组，计算所有数值列的平均值，结果生成新列Group.1存放基因符号
data_swab <- aggregate(data_swab, by = list(data_swab$`Gene Symbol`), mean)
data_plasma <- aggregate(data_plasma, by = list(data_plasma$`Gene Symbol`), mean)
data_platelet <- aggregate(data_platelet, by = list(data_platelet$`Gene Symbol`), mean)

# 1.3 设置行名与列处理
# 将分组后的基因符号设置为行名
rownames(data_swab) <- data_swab$Group.1
rownames(data_plasma) <- data_plasma$Group.1
rownames(data_platelet) <- data_platelet$Group.1

# 1.4 删除前9列（包含分组列Group.1和其他元数据列）
data_swab <- data_swab[,-1:-9]
data_plasma <- data_plasma[,-1:-9]
data_platelet <- data_platelet[,-1:-10]

# 1.5 数据结构转换
# 转置数据框：将行（基因）转为列，列（样本）转为行
data_swab <- as.data.frame(t(data_swab))
data_plasma <- as.data.frame(t(data_plasma))
data_platelet <- as.data.frame(t(data_platelet))

# 1.6 添加新列sampleid存放原始行名（即样本ID）
data_swab$sampleid <- rownames(data_swab)
data_plasma$sampleid <- rownames(data_plasma)
data_platelet$sampleid <- rownames(data_platelet)

# 1.7 添加实验分组信息
# 创建group列：前9个样本标记为"HD"，后7个样本标记为"OXG"
data_swab$group <- c(rep("HD",9),rep("OXG",7))
data_plasma$group <- c(rep("HD",9),rep("OXG",7))
data_platelet$group <- c(rep("HD",9),rep("OXG",7))

# data QC:
# 1. no abnormal sample data
# 2. no non-biological differences

library(ggplot2)
library(tidyr)

## 2.1 boxplot 
# 1. 创建数据副本(保留原始数据完整性)
data <- data_swab
data <- data_plasma
data <- data_platelet

# 2. 将分组列转换为有序因子
data$group <- factor( data$group, levels = c("OXG", "HD") )

# 3. 数据重塑：宽表变长表
data_ggplot <- tidyr::gather(data, key = "key", value = "value",-c("sampleid", "group")) 

# 4. 创建分组颜色映射
value_colour <- c("HD" = "#00A087FF", "OXG" = "#E64B35FF")        

# 5. 构建箱线图
ggplot(data_ggplot, aes(x = sampleid, y = log2(value + 1), fill = group)) + 
  geom_boxplot() +      # 核心箱线图图层
  scale_fill_manual(values = value_colour) +   # 应用自定义颜色方案
  theme_classic() +  # 无背景网格
  theme(axis.text.x = element_text( angle = 45,  hjust = 1, colour = "black", size = 10 ),
        axis.text.y = element_text( hjust = 1,  colour = "black", size = 10 )) + 
  labs(x = "") +
  ggtitle("swab") +  # 创建主标题
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold" ))

# 6. 保存
ggsave(filename = "QC_boxplot.pdf", height = 5, width = 5,plot = last_plot())

