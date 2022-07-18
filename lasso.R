
library("glmnet")
library("survival")

setwd("D:/Rcode/文章思路/坏死性凋亡/单因素和多因素cox")                #设置工作目录
rt=read.table("输入文件2.txt",header=T,sep="\t",row.names=1)       #读取文件
#rt$Survival_time=rt$Survival_time/365
rt$Survival_time=rt$Survival_time+0.04
rt$Survival_time=rt$Survival_time/365

gene=read.table("gene.txt",header=F)
#rt=rt[,c("Survival_time","Vital_status",as.vector(gene[,1]))]

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$Survival_time,rt$Vital_status))

fit <- glmnet(x, y, family = "cox", maxit = 1000)
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)

riskScore=predict(cvfit, newx = x, s = "lambda.min",type="response")
outCol=c("Survival_time","Vital_status",lassoGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(riskScore),risk)
write.table(cbind(id=rownames(outTab),outTab),
    file="lassoRisk.txt",
    sep="\t",
    quote=F,
    row.names=F)

###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######生信自学网: http://www.biowolf.cn/
######合作邮箱：2749657388@qq.com
######答疑微信: 18520221056
