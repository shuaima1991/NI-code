rm(list = ls())
setwd("D:/Rcode/文章思路/细胞焦亡与抗肿瘤免疫/泛癌与riskscore新/干性与风险评分泛癌")
gene="riskScore"           
geneRT=read.table("gene.txt",sep="\t",header=F)          
files='输入文件33.txt'
data=read.table(files, header=T,sep="\t",check.names=F,row.names = 1)

##错误点就在Cancertype没有转化问因子，索引没有32癌症列表
CancerType=levels(as.factor(data$CancerType))
sameGenes=sameGenes=intersect(as.vector(geneRT[,1]),colnames(data))
outTab=matrix(1,32,6)

corTab=matrix(1,32,6)

for(i in 1:32){
  data1=data[data$CancerType==CancerType[i],]
  x=as.numeric(data1[,gene])
  
  for(j in 1:6){
    y=as.numeric(data1[,j])
    corT=cor.test(x,y)
    cor=corT$estimate
    pValue=corT$p.value
    
    outTab[i,j]=pValue
    corTab[i,j]=cor
  }}

colnames(outTab)=sameGenes
row.names(outTab)=CancerType
write.table(outTab,file="geneCor.pvalue2.txt",sep="\t",row.names=T,quote=F)
colnames(corTab)=sameGenes
row.names(corTab)=CancerType
write.table(corTab,file="geneCor.cor2.txt",sep="\t",row.names=T,quote=F)

