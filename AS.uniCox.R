
pFilter=0.05                                                      #���嵥����������

setwd("D:/Rcode/����˼·/ϸ�������뿹��������/DEGs ͻ�䣬�׻����������������� 2") #���ù���Ŀ¼
library(survival)                                                 #���ð�
library(UpSetR)
rt=read.table("�����������ļ���3.txt",header=T,sep="\t",check.names=F,row.names=1)       #��ȡ�����ļ�

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
 coxSummary = summary(cox)
 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
 outTab=rbind(outTab,
              cbind(id=i,
                    z=coxSummary$coefficients[,"z"],
                    HR=coxSummary$conf.int[,"exp(coef)"],
                    HR.95L=coxSummary$conf.int[,"lower .95"],
                    HR.95H=coxSummary$conf.int[,"upper .95"],
                    pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
              )
}
#������е����صĽ��
outTab = outTab[is.na(outTab$pvalue)==FALSE,]
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab,file="uniCoxResult.txt",sep="\t",row.names=F,quote=F)
#��������������Ľ��
sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<pFilter,]
write.table(sigTab,file="uniCoxResult.Sig.txt",sep="\t",row.names=F,quote=F)
#�������������AS��PSIֵ�����ں�����ģ
sigGenes=c("futime","fustat")
sigGenes=c(sigGenes,as.vector(sigTab[,1]))
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)


#����upsetͼ
gene=sapply(strsplit(sigGenes,"\\|"),"[",1)
asType=sapply(strsplit(sigGenes,"\\|"),"[",3)
upsetList=list(AA=unique(gene[asType=="AA"]),
               AD=unique(gene[asType=="AD"]),
               AP=unique(gene[asType=="AP"]),
               AT=unique(gene[asType=="AT"]),
               ES=unique(gene[asType=="ES"]),
               ME=unique(gene[asType=="ME"]),
               RI=unique(gene[asType=="RI"]) )
upsetData=fromList(upsetList)

pdf(file="uniCoxUpset.pdf",onefile = FALSE,width=8,height=5)              #����ͼƬ
upset(upsetData,
      nsets = 7,                                    #չʾ�ɱ�������͸���
      order.by = "freq",                            #������Ŀ����
      show.numbers = "yes",                         #��״ͼ�Ϸ��Ƿ���ʾ��ֵ
      number.angles = 20,                           #����Ƕ�
      point.size = 1.5,                             #��Ĵ�С
      matrix.color="red",                           #��������ɫ
      line.size = 0.8,                              #��������
      mainbar.y.label = "Gene Intersections",
      sets.x.label = "Set Size")
dev.off()

######Video source: http://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺2749657388@qq.com
######����΢��: 18520221056