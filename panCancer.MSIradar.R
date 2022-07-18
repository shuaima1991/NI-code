
setwd("")             #设置工作目录

library(fmsb) 
#setwd("D:\\biowolf\\panCancer\\18.MSIradar")         #设置工作目录
data=read.table("fmsbInput.txt",header=T,sep="\t",row.names=1,check.names=F)   #读取输入文件
maxValue=ceiling(max(abs(data))*10)/10
data=rbind(rep(maxValue,ncol(data)),rep(-maxValue,ncol(data)),data)
#定义颜色
colors="blue"
#定义显著性
corStat=read.table("corStat.txt",header=T,sep="\t",row.names=1,check.names=F)
colnames(data)=paste0(colnames(data),corStat$sig)

#输出结果
pdf(file="radar.pdf",height=7,width=7)
radarchart( data, axistype=1 , 
    pcol=colors,                
    plwd=2 ,                     
    plty=1,                     
    cglcol="grey",             
    cglty=1,                    
    caxislabels=seq(-maxValue,maxValue,maxValue/2),   
    cglwd=1.2,                   
    axislabcol="green",          
    vlcex=0.8                   
)
dev.off()

