
setwd("D:/Rcode/文章思路/坏死性凋亡/泛癌分析/NI指数与MSI")             #设置工作目录

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
    pcol=colors,                 #设置颜色
    plwd=2 ,                     #线条粗线
    plty=1,                      #虚线，实线
    cglcol="grey",               #背景线条颜色
    cglty=1,                     #背景线条虚线，实线 
    caxislabels=seq(-maxValue,maxValue,maxValue/2),    #坐标刻度
    cglwd=1.2,                   #背景线条粗细
    axislabcol="green",           #刻度颜色
    vlcex=0.8                    #字体大小
)
dev.off()

