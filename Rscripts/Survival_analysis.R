library(readr)
library(tidyverse)
library(dplyr)
library(data.table)
library(survminer)
library(survival)
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "limma"))
library(TCGAbiolinks)
library(SummarizedExperiment)
###这里我提前下载了临床数据
clin <- fread("clinical.tsv",data.table = F)

clin_time <- clin %>% 
  dplyr::select(case_submitter_id,vital_status,days_to_death,days_to_last_follow_up,age_at_index,gender,ajcc_pathologic_t,ajcc_pathologic_m,ajcc_pathologic_n,ajcc_pathologic_stage) %>%
  dplyr::filter(!duplicated(case_submitter_id))

#创建新的一列OS.time，如果vital_status列是Alive，就用days_to_last_follow_up数据填充，如果vital_status列是Dead，就用days_to_death数据填充。
#创建一个OS列，生存用0表示，死亡用1表示。
clin_merge <- clin_time %>% 
  dplyr::mutate(OS.time = case_when(vital_status == "Alive" ~ days_to_last_follow_up,
                                    vital_status == "Dead" ~ days_to_death)) %>%
  dplyr::mutate(OS = case_when(vital_status == "Alive" ~ 0,
                               vital_status == "Dead" ~ 1))
head(clin_merge)
#使用rename函数修改列名
clin_merge1 <- clin_merge %>% 
  rename(
    age = age_at_index,  # 将age_at_index列名修改为age
    T = ajcc_pathologic_t,  # ajcc_pathologic_t列名修改为T
    M = ajcc_pathologic_m, 
    N = ajcc_pathologic_n, 
    stage = ajcc_pathologic_stage 
  )

clin<-clin_merge1[,-(3:4)]
clin<-na.omit(clin)
write.table(clin, file = "clinical.txt",sep = "\t",row.names = F,quote = F)
write.csv(clin, file = "clinical.csv", row.names = F, quote = FALSE)
#临床数据输出


# 查询并下载HNRNPC表达数据
query <- GDCquery(
  project = "TCGA-ESCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
data <- GDCprepare(query)

# 提取基因表达矩阵（HNRNPC的Ensembl ID: ENSG00000111731）先查找好对应的ID
exp_matrix <- assay(data, "unstranded")
rownames(exp_matrix) <- sub("\\..*", "", rownames(exp_matrix))
hnrnpc_exp <- exp_matrix["ENSG00000111731", ]


# 转换为数据框
hnrnpc_df <- data.frame(
  patient_id = colnames(exp_matrix),
  HNRNPC_expression = as.numeric(hnrnpc_exp)
)
head(hnrnpc_df)

survival_data<-clin
head(survival_data)

hnrnpc_df$case_submitter_id <- substr(hnrnpc_df$patient_id, 1, 12)
head(hnrnpc_df)

merged_data <- inner_join(hnrnpc_df, survival_data, by = "case_submitter_id") %>%
  dplyr::select(HNRNPC_expression, OS.time,OS,case_submitter_id)
merged_data$OS.time[is.na(merged_data$OS.time) & merged_data$OS == 0] <- max(merged_data$OS.time, na.rm = TRUE)
head(merged_data)
merged_data$OS.time<-as.numeric(merged_data$OS.time)
head(merged_data)
#输出merged_data
write.csv(merged_data, file = "OS_HNRNPC.csv", row.names = F, quote = FALSE)

set.seed(123)
cutpoint <- surv_cutpoint(
  merged_data,
  time = "OS.time",
  event = "OS",
  variables = "HNRNPC_expression"
)
?surv_cutpoint

print(cutpoint)
merged_data$group <- surv_categorize(cutpoint)$HNRNPC
write.csv(merged_data, file = "svrv_data_group.csv", row.names = F, quote = FALSE)



fit <- survfit(Surv(OS.time, OS) ~ group, data = merged_data)
head(fit)
ggsurvplot(
  fit,
  data = merged_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  palette = c("#E7B800", "#2E9FDF"),
  xlab = "Time (days)",
  legend.title = "HNRNPC Expression",
  legend.labs = c("High", "Low"),
  title = "TCGA-ESCA HNRNPC Survival Analysis"
)

ggsurvplot(fit, 
           data = merged_data, 
           size = 1, # change line size 
           cex.lab=2,
           break.time.by = 250, # break X axis in time intervals by 500.
           xlim = c(0,2134),
           axis.title.x =element_text(size=5), 
           axis.title.y = element_text(size=5),
           palette = c("#FF5044","#44B7BE"),# custom color palettes 
           #conf.int = TRUE, # Add confidence interval 
           pval = TRUE, # Add p-value
           pval.coord = c(1500, 0.9),
           risk.table = TRUE, # Add risk table 
           xlab = "Time(days)", # customize X axis label. 
           ylab="Survival probablity ",
           risk.table.col = "strata",# Risk table color by groups 
           legend.labs =  c("High","Low"), # Change legend labels
           risk.table.height = 0.3, # Useful to change when you have multiple groups 
           ggtheme = theme_bw(), # Change ggplot2 theme
           title = "TCGA-ESCA HNRNPC Survival Analysis",
           legend.title = "",
           ylim=c(0.4,1)
)
?ggsurvplot
###注意legend.labs和fit的映射关系！
cox_model <- coxph(Surv(OS.time, OS) ~ group, data = merged_data)
summary(cox_model)
?ggsurvplot

# 将生存时间（OS.time，天数）转换为月数（假设每月30天）
merged_data$time_month <- merged_data$OS.time / 30

# 查看转换后的数据示例
head(merged_data[, c("OS.time", "time_month")])
ggsurvplot(fit_month, 
           data = merged_data, 
           size = 1, 
           break.time.by = 6,        # 每6个月分割X轴
           xlim = c(0, 72),          # 原2134天≈71.13月，取整为72
           xlab = "Time (months)",   # X轴标签改为月
           ylab = "Survival probability",
           palette = c("#FF5044", "#44B7BE"),
           pval = TRUE, 
           risk.table = TRUE,
           risk.table.col = "strata",
           legend.labs = c("High", "Low"),
           risk.table.height = 0.3,
           ggtheme = theme_bw(),
           title = "TCGA-ESCA HNRNPC Survival Analysis (Monthly Time)",
           ylim = c(0.4, 1),
           legend.title = "",
           pval.coord = c(54, 0.85)
)
merged_data$time_year <- merged_data$OS.time / 365.25

# 重新拟合生存分析模型
fit_year <- survfit(Surv(time_year, OS) ~ group, data = merged_data)
head(merged_data)
ggsurvplot(fit_year,
           break.time.by = 1,          # 每年显示一个刻度
           xlab = "Time (Years)",
           risk.table = TRUE,
)
ggsurvplot(fit_year,
           data = merged_data, 
           size = 1, 
           break.time.by = 1,        # 每6个月分割X轴
                     # 原2134天≈71.13月，取整为72
           xlab = "Time (Years)",   # X轴标签改为月
           ylab = "Survival probability",
           palette = c("#FF5044", "#44B7BE"),
           pval = TRUE, 
           risk.table = TRUE,
           risk.table.col = "strata",
           legend.labs = c("High", "Low"),
           risk.table.height = 0.3,
           ggtheme = theme_bw(),
           title = "TCGA-ESCA HNRNPC Survival Analysis",
           ylim = c(0.4, 1),
           legend.title = "",
           pval.coord = c(4, 0.85)
)
