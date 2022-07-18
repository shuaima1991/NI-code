setwd("")
options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")

library(clusterProfiler)
library(limma)
library(ggplot2)
library(data.table)
library(ggpubr)
library(GSVA)
source("twoclasslimma.R")
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE)

samAnno <- read.table("easy_input_sample_annotation.txt", sep = "\t",header = T, check.names = F)

expr <- fread("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",sep = "\t",stringsAsFactors = F,check.names = F,header = T)
expr[1:3,1:3]
expr <- as.data.frame(expr); rownames(expr) <- expr[,1]; expr <- expr[,-1]
gene <- sapply(strsplit(rownames(expr),"|",fixed = T), "[",1)
expr$gene <- gene
expr <- expr[!duplicated(expr$gene),]
rownames(expr) <- expr$gene; expr <- expr[,-ncol(expr)]
expr[expr < 0] <- 0 # 对于这份泛癌数据，将略小于0的数值拉到0，否则不能取log（其他途径下载的泛癌数据可能不需要此操作）
colnames(expr) <- substr(colnames(expr),1,15) # 非TCGA数据可忽略这行
gc()

#genelist <- read.table("easy_input_genelist.txt", header = T)$GeneSymbol

genelist <- c("FADD","FAS","FASLG","MLKL","RIPK1","RIPK3","TLR3","TNF")


es <- gsva(expr = as.matrix(log2(expr + 1)),
           gset.idx.list = list("genelist" = genelist),
           method = "ssgsea",
           parallel.sz = 0) # 若采用LINUX或者MacOS可以设置为0，Windows请设置为1
write.table(es, file = "output_ssgsea enrichment score of interested gene list in pancancer.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


#es <- log2(expr[rownames(expr) == "TP53",] + 1) 

es <- read.table("output_ssgsea enrichment score of interested gene list in pancancer.txt", sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

msigdb.hallmark <- read.gmt("h.all.v7.2.symbols.gmt") 

pct <- 0.3 

gseaTab <- NULL


tumors <- c("BLCA","BRCA","CESC","CHOL","COAD",
            "ESCA","GBM","HNSC","KICH","KIRC",
            "KIRP","LIHC","LUAD","LUSC","PAAD",
            "PRAD","READ","STAD","THCA","UCEC") 
for (i in tumors) {
  message("--",i,"...")
  sam <- samAnno[which(samAnno$`cancer type` == i),"simple_barcode"]
  comsam <- intersect(colnames(es), sam) 
  
  tumsam <- comsam[substr(comsam,14,14) == "0"] 
  # tumsam <- comsam # 
  es_subset <- as.data.frame(t(es[,tumsam]))
  es_subset <- es_subset[order(es_subset$genelist,decreasing = T),,drop = F] 
  
  high.es <- rownames(es_subset[1:(nrow(es_subset) * pct),,drop = F]) # 取前30%为高组
  low.es <- rownames(es_subset[nrow(es_subset):(nrow(es_subset) - nrow(es_subset) * pct + 1),,drop = F]) # 取后30%为低组
  
  # 高低两组样本的limma差异表达分析
  subt <- data.frame(condition = rep(c("high","low"),c(length(high.es),length(low.es))),
                     row.names = c(high.es, low.es),
                     stringsAsFactors = F)
  gset <- log2(na.omit(expr[,rownames(subt)]) + 1)
  twoclasslimma(subtype  = subt, # 亚型信息（必须包含一列为condition）
                featmat  = gset, # 表达谱（会自动检测数据标准化与否）
                treatVar = "high", # “治疗组”的名字
                ctrlVar  = "low", # “对照组”的名字
                prefix   = paste0("TCGA_enrichment_score_",i), # 差异表达文件的前缀
                overwt   = T, # 决定是否覆盖已经存在的差异表达结果
                sort.p   = F, # 决定是否对结果按照FDR排序
                verbose  = TRUE, # 决定是否显示一些提示
                res.path = ".") # 确定结果路径
  
  # 加载差异表达分析结果
  res <- read.table(paste0("TCGA_enrichment_score_",i,"_limma_test_result.high_vs_low.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  
  # 产生pre-ranked基因列表
  res <- res[order(res$log2fc, decreasing = T),]
  glist <- res$log2fc; names(glist) <- rownames(res)
  
  # 运行GSEA富集分析
  set.seed(20211114)
  gsea <- GSEA(geneList = glist,
               pvalueCutoff = 1, # 得到所有结果
               seed = TRUE,
               TERM2GENE = msigdb.hallmark)
  gc()
  gsea.df <- as.data.frame(gsea) # 数据转换为数据框
  write.table(gsea.df,file = "output_gsea between high and low group of enrichment score in pancancer.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
  
  gseaTab <- rbind.data.frame(gseaTab,
                              data.frame(term = gsea.df$ID,
                                         NES = gsea.df$NES,
                                         FDR = gsea.df$p.adjust,
                                         tumor = i,
                                         stringsAsFactors = F),
                              stringsAsFactors = F)
}

# 将所有肿瘤的富集分析结果输出到文件
write.table(gseaTab, "output_TCGA_pancan_gsea_regarding_es_group.txt",sep = "\t",row.names = F,col.names = T,quote = F)


# 设置颜色
darkblue <- "#303B7F"
darkred <- "#D51113"
yellow <- "#EECA1F"

# 生成泡泡图
tmp <- gseaTab
tmp$term <- gsub("HALLMARK_","",tmp$term)
my_palette <- colorRampPalette(c(darkblue,yellow,darkred), alpha=TRUE)(n=128)
ggplot(tmp, aes(x=tumor,y=term)) +
  geom_point(aes(size=-log10(FDR),color=NES)) +
  scale_color_gradientn('NES', 
                        colors=my_palette) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, size = 12, hjust = 0.3, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,1,1), "lines"))
