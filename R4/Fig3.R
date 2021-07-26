library(visNetwork)
library(rpart.plot)
library(rpart)
library(circlize)
library(EpiDISH)
library(pheatmap)
library(ggpubr)
library(ggthemes)
library(ggplot2)
library(ggplotify)
library(pROC)


########### Figure 3B
m<-read.csv("TCGA.MeImm.Meth.txt",stringsAsFactors = F,row.names = 1)
IPMS<-read.delim("IPMS.txt",stringsAsFactors = F,header=F)
imc<-read.delim("sample_info.txt",stringsAsFactors = F)
imc<-imc[,c("sample_id","IMC")]
imc$IMC=gsub("IMC1|IMC2","Cold",imc$IMC)
imc$IMC=gsub("IMC3|IMC4|IMC5","Hot",imc$IMC)
sample<-split(imc$sample_id,imc$IMC)
sample$normal<-colnames(m)[grep("[.]1..$",colnames(m))]
bed<-data.frame(iDMC=rownames(m),
                hot=rowMeans(m[,sample$Hot]),
                cold=rowMeans(m[,sample$Cold]),
                normal=rowMeans(m[,sample$normal]))
IPMS$V8=IPMS$V3+1
bed<-merge(IPMS[,c(2,3,8,5,1)],bed,by.x="V1",by.y="iDMC")
bed<-bed[,c(2:4,1,5:8)]
bed=bed[order(bed$V3),]
bed$V5=gsub("Hot","#F7B731",bed$V5)
bed$V5=gsub("Cold","#8854D0",bed$V5)

methdiff<-read.delim("DMC_methdiff.txt",stringsAsFactors = F)
methdiff<-methdiff[which(methdiff$iDMC %in% bed$V1),]
methdiff<-merge(methdiff[,c(1,4,6)],bed[,1:4],by.x="iDMC",by.y="V1")
methdiff<-methdiff[,c(4:6,1:3)]
bed_list<-split(methdiff[,c(1:3,6)],methdiff$Type)

col=c("#1B3463","#4E9CC7","#D1E8F2","white","#FCE3D2","#F4A784","#C03339","#720A21")
col_fun = colorRamp2(c(0,0.5, 1), c("#1B3463","white","#720A21"))
pdf("circos.pdf",height = 10,width = 10)
circos.par("start.degree" = 260)
circos.initializeWithIdeogram(chromosome.index = paste0("chr", c(1:8,10:12,15:17,20,"X")))
circos.genomicLabels(bed, labels.column = 4, side = "inside",cex = 0.6,col = bed$V5)
circos.genomicHeatmap(bed[,-5:-4], col = col_fun, side = "inside",heatmap_height = 0.15)
#circos.genomicRainfall(bed_list, pch = 16, cex = 0.8, col = c("#FF000080", "#0000FF80"))
circos.genomicDensity(bed_list$Hypermethylated, col = c("#FF000080"), track.height = 0.1)
circos.genomicDensity(bed_list$Hypomethylated, col = c("#0000FF80"), track.height = 0.1)
circos.clear()
dev.off()

########### Figure 3C
m<-read.csv("TCGA.MeImm.Meth.txt",stringsAsFactors = F,header=F,row.names = 1)
m<-data.frame(t(m))
colnames(m)[1]<-"sample_id"
imc<-read.delim("sample_info.txt",stringsAsFactors = F)
imc<-imc[,c("sample_id","IMC")]
imc$IMC=gsub("IMC1|IMC2","Cold",imc$IMC)
imc$IMC=gsub("IMC3|IMC4|IMC5","Hot",imc$IMC)
mm<-merge(imc,m,by="sample_id")
nm<-m[grep("[.]1..$",m$sample_id),]
nm<-data.frame(sample_id=nm$sample_id,IMC="normal",nm[,-1])
m<-rbind(mm,nm)
write.table(m,"Tree.data.txt",quote=F,col.names = T,row.names = F,sep="\t")
m<-read.delim("Tree.data.txt",stringsAsFactors = F)
m<-m[,-1]

m.rpart <- rpart(IMC ~ .-IMC, data = m,minsplit=10)
rpart.plot(m.rpart, digits = 3,)
visTree(m.rpart,main = "title",
        fallenLeaves = T,
        simplifyRules=F,
        colorEdges="#d1d8e0",nodesFontSize = 20,
        edgesFontSize = 18,
        colorY = c("blue","red","green"))
        

########### Figure 3D
load("iPMS_sigMatrix.RData") ##sigm
ip.infer <- epidish(beta.m = rexp, 
				  ref.m = sigm,
				  method = "RPC")$estF 
ipmsc<-ip.infer[,2]-ip.infer[,1]
names(ipmsc)<-rownames(ip.infer)
ip.infer <- c("Cold","Hot")[apply(ip.infer, 1, which.max)]
names(ip.infer)<-colnames(rexp)
imc$ip.infer=ip.infer[imc$sample_id]
roc.data<-data.frame(imc=imc$IMC,ipmsc)
ggplot(roc.data,aes(x=imc,y=ipmsc,fill=imc))+
  geom_violin(alpha=0.5)+
  geom_boxplot(width=0.1,position=position_dodge(0.8),
               outlier.alpha = 0,alpha=0.5)+
  theme_few()+
  scale_fill_manual(values=c("#6AB5DD","#EA614A"))+
  guides(fill=FALSE)+
  xlab(NULL)+
  ylab("iPMS Score")+
  geom_hline(yintercept= -0.253,linetype="dashed",color="#C30D23",size=0.5)+
  scale_y_continuous(breaks = c(-1,-0.253,1))  
ggsave("ipmsScore.pdf",height=3,width=3)

########### Figure 3E
pdf("ROC.pdf",height=4,width=4)
roc(roc.data$imc,roc.data$ipmsc, print.thres=TRUE, percent=F,
    ci=TRUE, boot.n=100, ci.alpha=0.95, stratified=FALSE,
    plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE)
dev.off()

########### Figure 3F
load("TCGA_imc.test.RData")
sv<-imc.test[,c("IMC","sample_id",
				"TCGA.PanCanAtlas.Cancer.Type.Acronym",
              "Overall.Survival..Months.","Overall.Survival.Status")]
sv$IMC=ip.infer[sv$sample_id]
lapply(c("Overall.Survival"), function(x){
	ssv<-sv[,c(1,2,grep(x,colnames(sv)),11)]
	loc<-union(union(which(is.na(ssv[,3])),
					which(is.na(ssv[,4]))),
			  which(is.na(ssv[,5])))

	if(length(loc)>0){ssv<-ssv[-loc,]}
	colnames(ssv)[2:4]<-c("cancer","time","status")
	ssv$status=as.numeric(gsub(":.*","",ssv$status))
	ssv$time=as.numeric(as.character(ssv$time))

	library("survival")
	library("survminer")
	library(ggplot2)
	library(ggthemes)
	library(ggsci)
	fit <- survfit(Surv(time, status) ~ IMC, data=ssv)   
			   ggsurvplot(fit,data=ssv,
						  pval = T, conf.int = T, 
						  palette=pal_d3()(5), 
						  risk.table.col = "strata", # Change risk table color by groups
						  linetype = "strata", # Change line type by groups
						  surv.median.line = "hv", # Specify median survival
						  ggtheme = theme_few()
			   )
	ggsave(paste0(x,"_IMC_survival.pdf"),height = 5,width = 4.5,dpi=600)
})

########### Figure 3G
pheat<-imc.test[,c(33:41,43:48)]
rownames(pheat)<-imc.test$sample_id
index<-data.frame(V1=colnames(pheat),
                  type=c("Immunoactivity","Immunoactivity",
                         rep("Immune_signature",5),
                         "Immunoactivity","Immunoactivity",
                         rep("TIMER",6)))
index$V1=as.character(index$V1)
index<-split(index$V1,index$type)
pheat=t(pheat)[,c(is$Cold,is$Hot)]
heatplot<-lapply(index, function(x){
  ann_col=data.frame(row.names=names(ip.infer[imc.test$sample_id]),
                     immunophenotype=ip.infer[imc.test$sample_id])	
  is<-split(names(ip.infer[imc.test$sample_id]),ann_col$immunophenotype)
  ann_colors = list(immunophenotype=c(Cold="#6AB5DD",Hot="#EA614A"))
  
  p<-pheatmap(pheat[x,],
              filename=paste(x[1],".pdf"),
              scale="row",
              annotation_col = ann_col,
              cluster_cols = F,
              cellheight = 9,
              #cellwidth = 5,
              fontsize = 7,
              show_colnames = F,
              annotation_colors = ann_colors, 
              gaps_col = length(is$Cold),
              color = colorRampPalette(c("#1B3463","#4E9CC7","#D1E8F2","white","#FCE3D2","#F4A784","#C03339","#720A21"))(150)
  )
  p<-as.ggplot(p)
  p
})
cowplot::plot_grid(heatplot$Immune_signature,
                   heatplot$Immunoactivity,
                   heatplot$TIMER, nrow = 3)

########### Figure 3H
ipmsc<-ipmsc[colnames(pheat)]
gd.cor<-reshape2::melt(pheat[index$TIMER,])
gd.cor$Var2=as.character(gd.cor$Var2)
gd.cor$ipmsc=ipmsc[gd.cor$Var2]
gd.cor$immunophenotype=ip.infer[gd.cor$Var2]
gd.cor$Var1=gsub("[(].*[)]","",gd.cor$Var1)
p2<-ggplot(gd.cor,aes(x=ipmsc,y=value,color=immunophenotype))+
  geom_point()+
  geom_smooth(color="black",size=0.5)+
  stat_cor(color="black")+
  theme_few()+
  scale_color_manual(values=c(Cold="#6AB5DD",Hot="#EA614A"))+
  ylab("Cell abundance estimated by TIMER")+
  xlab("iPMS score")+
  theme(legend.position = "none")+
  facet_wrap(Var1~.,scales = "free",ncol=3)


