# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
library(dplyr)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

## 转录本cDNA长度
if(T){
  # 拿到所有基因的所有转录本,CDS,5UTR,3UTR等长度信息
  t_l=transcriptLengths(txdb)
  t_l=transcriptLengths(txdb,with.cds_len=TRUE,
                        with.utr5_len=TRUE, 
                        with.utr3_len=TRUE)
  head(t_l)
  t_l=na.omit(t_l)
  
  #增加一列，mRNA长度
  t_l <- mutate(t_l,mRNA_len = cds_len + utr5_len + utr3_len)
  
  # 按基因名，mRNA长度排序
  t_l=t_l[order(t_l$gene_id,t_l$tx_len,decreasing = T),]
  str(t_l)
  library(data.table)
  t_l_1 <- setDT(t_l)  # 转换为data.table对象
  
  t_l_max <- t_l_1[
    , .SD[which.max(tx_len)], 
    by = gene_id
  ]
  
  
  ## 取最长的mRNA代表转录本加工后的长度
  # 已经排序，第一位最长，去除后面冗余
  t_l_max=t_l[!duplicated(t_l$gene_id),] 
  head(t_l_max)
  t_l_max=t_l_max[,c(1,5)]
  t_l_max=t_l_max[order(t_l_max$gene_id),]
  head(t_l_max)
  ## 取cds平均长度位基因长度
  library(dplyr)
  
  t_l_mean <- t_l %>% 
    group_by(gene_id) %>% 
    summarise(
      avg_len = mean(tx_len, na.rm = TRUE)) 
  head(t_l_mean)
  
  ## 取mRNA长度中位数作为基因长度
  t_l_median <- t_l %>% 
    group_by(gene_id) %>% 
    summarise(
      median_len = median(tx_len, na.rm = TRUE)) 
  
  head(t_l_median)
  
}

head(t_l_max)

## id转基因symbol
library(org.Hs.eg.db)
ls("package:org.Hs.eg.db")  
s2g=toTable(org.Hs.egSYMBOL)
head(s2g)
t_l_max=merge(t_l_max,s2g,by='gene_id')
head(t_l_max)
t_l_mean=merge(t_l_mean,s2g,by='gene_id')
t_l_median=merge(t_l_median,s2g,by='gene_id')
head(t_l_mean)
print(genes[2,"Symbol"])
class(genes[2,"Symbol"])
print(t_l_median[which(t_l_max$symbol=="SLC2A3"),])

t_l_max=t_l[!duplicated(t_l$gene_id),] 
head(t_l_max)
t_l_max=t_l_max[,c(3,9)]
t_l_max=t_l_max[order(t_l_max$gene_id),]
head(t_l_max)
t_l_max[which(t_l_max$symbol=='BRAF'),]



for (i in 1:nrow(genes)){
  genes$transcDNA_MAX <- t_l_max[which(t_l_max$symbol==genes[i,"SYMBOL"]),]
  
}
# 预分配列表存储匹配结果
genes$transcDNA <- vector("list", nrow(genes))
match_idx <- which(t_l_max[which(t_l_max$symbol==genes[i,"Symbol"]),])
print(t_l_max[match_idx, ])
current_symbol <- genes$Symbol[1]
match_idx <- which(t_l_max$symbol == current_symbol)



for (i in 1:nrow(genes)) {
  current_symbol <- genes$Symbol[i]
  match_idx <- which(t_l_max$symbol == current_symbol)
  if (length(match_idx) > 0) {
    genes$transcDNA_MAX[[i]] <- t_l_max[match_idx, ]
    genes$transcDNA_mean[[i]] <- t_l_mean[match_idx, ]
    genes$transcDNA_median[[i]] <- t_l_median[match_idx, ]
  } else {
    genes$transcDNA_MAX[[i]] <- NA  # 无匹配时填充NA
    genes$transcDNA_mean[[i]] <- NA
    genes$transcDNA_median[[i]] <- NA
  }
}

for (i in 1:nrow(genes)) {
  current_symbol <- genes$Symbol[i]
  match_idx <- which(t_l_max$symbol == current_symbol)
  if (length(match_idx) > 0) {
    genes$transcDNA_MAX[[i]] <- t_l_max[match_idx, ]
    
  } else {
    genes$transcDNA_MAX[[i]] <- NA  # 无匹配时填充NA
   
  }
}


library(tidyr)
genes_expanded <- genes %>%
  unnest(transcDNA_MAX, names_sep = "_MAX_", keep_empty = TRUE) %>%
  unnest(transcDNA_mean, names_sep = "_MEAN_", keep_empty = TRUE) %>% 
  unnest(transcDNA_median, names_sep = "_MEDIAN_", keep_empty = TRUE) %>%
  distinct()

# 新版tidyr语法（推荐）
genes_expanded <- genes %>% 
  unnest(cols = c(transcDNA_MAX, transcDNA_mean))
# 写入 CSV
write.csv(genes_expanded, "Trans_cDNA_expanded003.csv", row.names = FALSE)
write.csv(genes,"Trans_cDNA.csv")
getwd()
head(genes)
genes$transcDNA <- NULL


