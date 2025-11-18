# 基因Id转换方法一 详见：https://mp.weixin.qq.com/s/dSWmMBGqtpMMNh4L49decA?forceh5=1
library(AnnotationDbi)
# mapIds函数基本用法mapIds(x, keys, column, keytype, ..., multiVals)# x：AnnotationDb对象# keys：从数据库中选择记录的key# column：要搜索的列（针对mapId），与columns的不同之处在于，值只能为单个元素# keytype：与所使用的key相匹配的key类型
library(org.Hs.eg.db)# 使用keytypes函数查看org.Hs.eg.db包支持的id类型> 
keytypes(org.Hs.eg.db) 
#[1] "ACCNUM"  "ALIAS" "ENSEMBL" "ENSEMBLPROT"  
#[5] "ENSEMBLTRANS" "ENTREZID" "ENZYME" "EVIDENCE" 
#[9] "EVIDENCEALL" "GENENAME" "GO"  "GOALL" 
#[13] "IPI" "MAP" "OMIM"  "ONTOLOGY" 
#[17] "ONTOLOGYALL" "PATH" "PFAM" "PMID" 
#[21] "PROSITE" "REFSEQ" "SYMBOL" "UCSCKG" 
#[25] "UNIGENE" "UNIPROT"
# ENSEMBL表示Ensembl Gene ID，ENSEMBLPROT表示Ensembl Protein ID，ENSEMBLTRANS表示Ensembl Transcript ID
# 将Gene Symbol转换为Gene Name，key为TP53，属于Gene Symbol，所以对应的keytype为SYMBOL> 
ENTREZID<-mapIds(x = org.Hs.eg.db,keys = genes$id,column = "ENTREZID",keytype = "ENSEMBL")
#'select()' returned 1:1 mapping between keys and columns  TP53 "tumor protein p53"
# 此处对人类基因进行ID转换，所以使用人类基因组注释包org.Hs.eg.db，如果想要对小鼠（org.Mm.eg.db）或其他物种进行基因ID转换，请加载对应物种的注释R包。# 同时对多个基因进行ID转换，将基因id放到一个向量中即可，转换后可使用na.omit函数删除转换失败的结果。
output_ENTREZID<-file("output_ENTREZID.csv","w")
writeLines(ENTREZID,output_ENTREZID)
close(output_ENTREZID)



library(RCurl)
library(stringr)
library(XML)
library(clusterProfiler)

rm(list=ls())
# 读入基因列表：
genes <- read.table("test_genes",header = T,stringsAsFactors = F)

# 基因Id转换方法二
# 将gene symbol转为entrze ID:
genes <- bitr(genes$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# 网址数据框：
genes$NCBI_url <- paste("https://www.ncbi.nlm.nih.gov/gene/",genes$ENTREZID,sep="")
head(genes)
class(genes)

# Official Symbol: //*[@id="summaryDl"]/dd[1]
# Official Full name的xpath：//*[@id="summaryDl"]/dd[2]
# HGNC ID的xpath：//*[@id="summaryDl"]/dd[3]/a
# Gene type的xpath：//*[@id="summaryDl"]/dd[5]/text()
# Summary的xpath：//*[@id="summaryDl"]/dd[10]/text()


# 根据xpath获取节点内容：
getNodesTxt <- function(html_txt1,xpath_p){
  els1 = getNodeSet(html_txt1, xpath_p)
  # 获得Node的内容，并且去除空字符：
  els1_txt <- sapply(els1,xmlValue)[!(sapply(els1,xmlValue)=="")]
  # 去除\n：
  str_replace_all(els1_txt,"(\\n )+","")
}

# 处理节点格式，为character且长度为0的赋值为NA：
dealNodeTxt <- function(NodeTxt){
  ifelse(is.character(NodeTxt)==T && length(NodeTxt)!=0 , NodeTxt , NA)
}


for(i in 1:nrow(genes)){
  # 获得网址：
  doc <- getURL(genes[i,"NCBI_url"])
  cat("成功获得网页！\t")
  # 获得网页内容
  html_txt1 = htmlParse(doc, asText = TRUE)
  
  # 获得Full Name:
  genes[i,"FullName"] <- dealNodeTxt(getNodesTxt(html_txt1,'//*[@id="summaryDl"]/dd[2]/text()'))
  cat("写入基因\t")
  
  genes[i,"Exon Count"] <- dealNodeTxt(getNodesTxt(html_txt1,'//*[@id="padded_content"]/div[6]/div[2]/div[2]/div/div/dl/dd'))
  cat("写入Exon Count\t")
  
  print(paste("完成第",i,"个了！"))
  
}

write.csv(genes,"output_EXON0.csv")

