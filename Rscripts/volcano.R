your_data <- data.frame(
  Symbol = c("HNRNPC", "SLC2A3", "SCD", "RCN2", "DEGS1", "PCSK6", "CD44", "STAG2", "KLK5", "MAL2", "GGT5"),
  log2.fc. = c(-3.848, 3.436, -3.098, -2.636, -2.426, -2.147, -2.131, -2.011, -2.408, -2.225, -2.002),
  PValue = c(6.07e-187, 3.65e-165, 5.42e-100, 1.54e-82, 3.95e-69, 2.35e-67, 4.54e-64, 2.21e-59, 6.09e-59, 6.03e-57, 7.57e-56),
  FDR = c(1.05e-182, 3.17e-161, 3.13e-96, 6.69e-79, 1.37e-65, 6.80e-64, 1.12e-60, 4.78e-56, 1.17e-55, 1.05e-53, 1.19e-52)
)
library(ggVolcano)

data<-shNC_vs_shHNRNPC2t
data$log2.fc. <- data$`log2(fc)`
# 查看列名和类型
names(data)
class(data$log2.fc.)  # 应为 "numeric"
class(data$FDR)       # 应为 "numeric"
data$log2.fc.<-data$`log2(fc)`


# 转换数据类型（若FDR为字符型）

data$FDR <- as.numeric(data$FDR)
data$log2.fc. <- as.numeric(data$log2.fc.)
data$FDR[data$FDR == 0] <- 1e-300

data <- add_regulate(data,
                     log2FC_name = "log2.fc.",  # 对应你的log2(fc)列
                     fdr_name = "FDR",          # 对应你的FDR列
                     log2FC = 1,                # log2倍变化阈值
                     fdr = 0.05)                # 显著性阈值

data <- add_regulate(data,
                     log2FC_name = "log2FoldChange",  # 对应你的log2(fc)列
                     fdr_name = "padj",          # 对应你的FDR列
                     log2FC = 1,                # log2倍变化阈值
                     fdr = 0.05)                # 显著性阈值

names(data)
p <- ggvolcano(
  data,
  x = "log2FoldChange",     # 指定log2FC列
  y = "padj",          # 指定FDR列
  label = "Symbol",   # 基因标签列
  label_number = 20,  # 显示前10个显著基因（按FDR排序）
  fills = c("#44B7BE", "grey90", "#FF5044"),  # 自定义颜色（上调/无/下调）
  colors = c("#44B7BE", "grey90", "#FF5044"),
  log2FC_cut = 1,     # 虚线标记log2FC阈值
  FDR_cut = 0.05      # 虚线标记FDR阈值
)

ggsave("volcano.png", p, width = 8, height = 6, dpi = 300)
gradual_volcano(data, 
                x = "log2.fc.", 
                y = "FDR",
                fills = RColorBrewer::brewer.pal(9, "RdYlBu"),  # 红黄蓝渐变
                colors = RColorBrewer::brewer.pal(9, "RdYlBu"),
                label = "Symbol",
                add_line = TRUE,
                log2FC_cut = 2,
                FDR_cut = 0.05,
                pointSizeRange = c(1),
                add_label = TRUE,                   # 确保标签显示
                label_number = 20                 # 显示更多标签以观察效果
                )  # 点大小根据显著性动态调整[3](@ref)
term_data <- data.frame(
  Symbol = c("HNRNPC", "SLC2A3", "SCD"),
  Term = c("RNA splicing", "Glucose transport", "Lipid metabolism")
)
term_volcano(data, term_data, 
             x = "log2.fc.", 
             y = "FDR",
             label = "Symbol",
             term_col = "Term")  # 标记通路基因[3](@ref)






