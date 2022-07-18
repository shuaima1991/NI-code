

library(limma)
library(survival)
library(survminer)
library(forestplot)

pFilter=0.05                                                      #km方法pvalue过滤条件
cliFile="Survival_SupplementalTable_S1_20171025_xena_sp"          #临床数据文件名称
setwd("D:/Rcode/文章思路/坏死性凋亡/泛癌分析/泛癌dFI")             #设置工作目录

#读取基因表达文件
rt=read.table("输入文件33.txt",header=T,sep="\t",check.names=F,row.names=1)
gene=colnames(rt)[1]

#读取临床数据文件
cli=read.table(cliFile,header=T,sep="\t",check.names=F,row.names=1)
cli=cli[,c("DFI.time","DFI")]
cli=na.omit(cli)
colnames(cli)=c("futime","fustat")

#对肿瘤类型进行循环
outTab=data.frame()
for(i in levels(as.factor(rt[,"CancerType"]))){
	rt1=rt[(rt[,"CancerType"]==i),]
	#交集基因
	data=cbind(rt1,gene=rt1[,gene])
	data=as.matrix(data[,c(gene,"gene")])
	if(nchar(row.names(data)[1])!=nchar(row.names(cli)[1])){
		row.names(data)=gsub(".$","",row.names(data))}
	data=avereps(data)
	sameSample=intersect(row.names(data),row.names(cli))
	sameData=data[sameSample,]
	sameCli=cli[sameSample,]
	rt1=cbind(sameCli,sameData)
	rt1$futime=rt1$futime/365
	if(nrow(rt1)<3){next}

	group=ifelse(rt1[,gene]>median(rt1[,gene]),"high","low")
	diff=survdiff(Surv(futime, fustat) ~group,data = rt1)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<pFilter){
		if(pValue<0.001){
			pValue="p<0.001"
		}else{
			pValue=paste0("p=",sprintf("%.03f",pValue))
		}
		fit <- survfit(Surv(futime, fustat) ~ group, data = rt1)
		#绘制生存曲线
		surPlot=ggsurvplot(fit, 
				    data=rt1,
				    title=paste0("Cancer: ",i),
				    pval=pValue,
				    pval.size=6,
				    legend.labs=c("high","low"),
				    legend.title=paste0(gene," levels"),
				    font.legend=12,
				    xlab="Time(years)",
				    ylab="Disease-free interval",
				    break.time.by = 1,
				    palette=c("red","blue"),
				    conf.int=F,
				    fontsize=4,
				    risk.table=TRUE,
				    risk.table.title="",
				    risk.table.height=.25)
		pdf(file=paste0("DFI.",i,".pdf"),onefile = FALSE, width = 6, height =5)
		print(surPlot)
		dev.off()
	}
	#cox分析
	cox=coxph(Surv(futime, fustat) ~ gene, data = rt1)
	coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	outTab=rbind(outTab,
	             cbind(cancer=i,
	                   HR=coxSummary$conf.int[,"exp(coef)"],
	                   HR.95L=coxSummary$conf.int[,"lower .95"],
	                   HR.95H=coxSummary$conf.int[,"upper .95"],
			           pvalue=coxP) )
}
write.table(outTab,file="cox.result.txt",sep="\t",row.names=F,quote=F)    #输出肿瘤类型和p值表格文件

############绘制森林图函数############
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
    #读取输入文件
	rt=read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
	data=as.matrix(rt)
	HR=data[,1:3]
	hr=sprintf("%.3f",HR[,"HR"])
	hrLow=sprintf("%.3f",HR[,"HR.95L"])
	hrHigh=sprintf("%.3f",HR[,"HR.95H"])
	pVal=data[,"pvalue"]
	pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
	clrs <- fpColors(box=forestCol,line="darkblue", summary="royalblue")      #定义颜色
	tabletext <- 
	  list(c(NA, rownames(HR)),
	       append("pvalue", pVal),
	       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )   #定义图片文字
	pdf(file=forestFile,width = 9,height = 6,onefile = FALSE)
	forestplot(tabletext, 
	           rbind(rep(NA, 3), HR),
	           col=clrs,
	           graphwidth=unit(50, "mm"),
	           xlog=T,
	           lwd.ci=4,
	           boxsize=0.6,
	           xlab="Hazard ratio",
	           txt_gp=fpTxtGp(ticks=gpar(cex=1.1),xlab=gpar(cex = 1.25))
	           )
	dev.off()
}
############绘制森林图函数############

bioForest(coxFile="cox.result.txt",forestFile="forest.pdf",forestCol="red")


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
