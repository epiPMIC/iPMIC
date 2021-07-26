library(survival)
library(survminer)
library(ggplot2)
library(ggthemes)
library(ggsci)
library(ggpubr)


################### Figure 4A-E
pheatX<-function(rexp,filename){
	library(EpiDISH)
	library(pheatmap)
	library(ggpubr)
	library(ggthemes)
	library(ggplot2)
	library(ggplotify)
	
	load("iPMS_sigMatrix.RData") ##sigm
	ip.infer <- epidish(beta.m = rexp, 
                      ref.m = sigm,
                      method = "RPC")$estF 
	ipmsc<-ip.infer[,2]-ip.infer[,1]
	names(ipmsc)<-rownames(ip.infer)
	ip.infer <- c("Cold","Hot")[apply(ip.infer, 1, which.max)]
	names(ip.infer)<-colnames(rexp)
		
	sigm<-as.matrix(read.delim("lung_NSCLC_not.specified_v2_Signature.txt.txt",row.names=1))
	
	
	BloodFrac.m <- epidish(beta.m = rexp, 
						   ref.m = sigm,
						   method = "RPC")$estF 
	pheat<-t(BloodFrac.m[,c(1:5)])
	rownames(pheat)<-c("CD14+ monocyte lineage", "CD19+ B-lymphocytes", 
						"CD4+ T-cells","CD56+ NK cells", "CD8+ cytotoxic T-lymphocytes")
	pheat<-pheat[,names(sort(ip.infer))]
	ann_col=data.frame(row.names=names(ip.infer),immunophenotype=ip.infer)	
	is<-split(names(ip.infer),ann_col$immunophenotype)
	sig<-apply(pheat, 1, function(x){
		w<-wilcox.test(x[is$Hot],x[is$Cold],alternative = "greater")
		lab<-c(0.0001,0.001,0.05,1,10)
		names(lab)<-c("***","**","*","NS","NS")
		return(names(which(w$p.value<lab))[1])
	  })
	rownames(pheat)=paste0(gsub("_TIMER","",rownames(pheat)),"(",sig,")")
	
	ann_colors = list(immunophenotype=c(Cold="#6AB5DD",Hot="#EA614A"))
	
	p<-pheatmap(scale(pheat),
			   filename=filename,
			   scale="row",
			   annotation_col = ann_col,
			   cluster_cols = F,
			   #cellheight = 9,
			   #cellwidth = 5,
			   #fontsize = 7,
			   show_colnames = F,
			   annotation_colors = ann_colors, 
			   gaps_col = length(is$Cold),
			   color = colorRampPalette(c("#1B3463","#4E9CC7","#D1E8F2","white","#FCE3D2","#F4A784","#C03339","#720A21"))(150)
			   )
	p<-as.ggplot(p)
	ipmsc<-ipmsc[colnames(pheat)]
	gd.cor<-reshape2::melt(pheat)
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
		  ylab("Cell abundance estimated\nby MethylCIBERSORT")+
		  xlab("iPMS score")+
		  theme(legend.position = "none")+
		  facet_wrap(Var1~.,scales = "free",ncol=5)+
		  ggtitle(paste0(unlist(strsplit(filename,"[.]"))[1]," (n = ",ncol(pheat),")"))
	cowplot::plot_grid(p2, p, ncol=2)
	ggsave(filename,height = 3,width = 25,dpi=600)
	return(list(ip.infer,ipmsc))
}

load("GSE119144.RData") ###rexp
ip.infer<-pheatX(rexp,filename="GSE119144.png")
save(ip.infer,file="GSE119144.ip.infer.RData")

load("GSE115246.RData") ###rexp
pheatX(rexp=rexp,filename="GSE115246.png")

load("GSE56044.RData") ###rexp
pheatX(rexp=rexp,filename="GSE56044.png")

load("GSE66836.RData") ###rexp
pheatX(rexp=rexp,filename="GSE66836.png")

load("GSE39279.RData") ##rexp
ip.infer<-pheatX(rexp=rexp,filename="GSE39279.png")
save(ip.infer,file="GSE39279.ip.infer.RData")

################### Figure 4F
load("GSE39279.ip.infer.RData")
ipmsc<-ip.infer[[2]]
ip.infer<-ip.infer[[1]]
sample<-read.delim("sample.txt",stringsAsFactors = F)
sample$ipmsc=ipmsc[sample$sample]
sample$ip.infer=ip.infer[sample$sample]
sample$recurrence=gsub("^T.*",NA,sample$recurrence)
sample$recurrence=gsub("recurrence: ","",sample$recurrence)
sample<-sample[-which(is.na(sample$recurrence)),]
sample<-sample[-which(sample$recurrence==""),]
sample$recurrence<-as.numeric(as.factor(sample$recurrence))-1
sample$Stage=as.numeric(as.factor(sample$Stage))

my_comparisons=list(c("Cold","Hot"))
tz<-ggplot(sample,aes(y=tumor.size_cm,x=ip.infer,fill=ip.infer))+
  geom_violin(alpha=0.5)+
  geom_boxplot(width=0.1,position=position_dodge(0.8),
               outlier.alpha = 0,alpha=0.5)+
  theme_few()+
  scale_fill_manual(values=c("#6AB5DD","#EA614A"))+
  guides(fill=FALSE)+
  xlab(NULL)+
  ylab("Tumor size (cm)")+
  stat_compare_means(comparisons=my_comparisons,
                     label = "p.signif",
                     method = "wilcox.test",
                     method.args = list(alternative = "greater"))
################### Figure 4G
ts<-ggplot(sample,aes(y=Stage,x=ip.infer,fill=ip.infer))+
  geom_violin(alpha=0.5)+
  geom_boxplot(width=0.1,position=position_dodge(0.8),
               outlier.alpha = 0,alpha=0.5)+
  theme_few()+
  scale_fill_manual(values=c("#6AB5DD","#EA614A"))+
  guides(fill=FALSE)+
  xlab(NULL)+
  ylab("Tumor stage")+
  stat_compare_means(comparisons=my_comparisons,
                     label = "p.signif",
                     method = "wilcox.test",
                     method.args = list(alternative = "greater"))
library(patchwork)
tz+ts
ggsave("tumor_size.pdf",height=3,width=6)

################### Figure 4H
fit <- survfit(Surv(adjuvant.chemotherapy, recurrence) ~ ip.infer, data=sample)   
res_cox<-coxph(Surv(adjuvant.chemotherapy, recurrence) ~ ip.infer, data=sample)
summary(res_cox)
ggsurvplot(fit,data=sample,
           pval = T, conf.int = T, 
           palette=c("#6AB5DD","#EA614A"), 
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_few()
)
ggsave("adjuvant.chemotherapy_IMC_survival.pdf",height = 5,width = 4.5,dpi=600)









