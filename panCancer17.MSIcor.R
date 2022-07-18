

setwd("D:/Rcode/文章思路/坏死性凋亡/泛癌分析/NI指数与MSI")             #设置工作目录

#读取表达文件
exp=read.table("输入文件33.txt", header=T,sep="\t",row.names=1,check.names=F)
#读取MSI文件
MSI=read.table("MSI.txt", header=T,sep="\t",row.names=1,check.names=F)
#去除正常样品
group=sapply(strsplit(row.names(exp),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
exp=exp[group==0,]
#样品取交集
sameSample=intersect(row.names(MSI),row.names(exp))
MSI=MSI[sameSample,]
exp=exp[sameSample,]

#相关性检验
outTab=data.frame()
fmsbTab=data.frame()
for(i in levels(as.factor(exp[,"CancerType"]))){
     exp1=exp[(exp[,"CancerType"]==i),]
     MSI1=MSI[(MSI[,"CancerType"]==i),]
	 x=as.numeric(MSI1[,1])
	 y=as.numeric(exp1[,1])
	 corT=cor.test(x,y,method="spearman")
	 cor=corT$estimate
	 pValue=corT$p.value
	 sig=ifelse(pValue<0.001,"***",ifelse(pValue<0.01,"**",ifelse(pValue<0.05,"*"," ")))
	 outTab=rbind(outTab,cbind(CancerType=i,cor=cor,pValue=pValue,sig))
	 fmsbTab=rbind(fmsbTab,cbind(CancerType=i,cor=cor))
}
write.table(outTab,file="corStat.txt",sep="\t",row.names=F,quote=F)
write.table(t(fmsbTab),file="fmsbInput.txt",sep="\t",col.names=F,quote=F)

