###############Fig4A
pheatX<-function(rexp,filename){
	library(EpiDISH)
	library(pheatmap)
	library(ggpubr)
	library(ggthemes)
	library(ggplot2)
	library(ggplotify)
	library(patchwork)
	
	load("iPMS_sigMatrix.RData") ##sigm
	ip.infer <- epidish(beta.m = rexp, 
                      ref.m = sigm,
                      method = "RPC")$estF 
	ip.infer=data.frame(sample=rownames(ip.infer),ip.infer)
	ip.infer$ipms.s=ip.infer$Hot-ip.infer$Cold
	ip.infer$type=NA
	ip.infer$type[which(ip.infer$ipms.s>0)]<-"Hot"
	ip.infer$type[which(ip.infer$ipms.s<0)]<-"Cold"
	if(length(is.na(ip.infer$type))>0){ip.infer<-ip.infer[-which(is.na(ip.infer$type)),]}

					  
	sigm<-as.matrix(read.delim("lung_NSCLC_not.specified_v2_Signature.txt.txt",row.names=1))
		
	BloodFrac.m <- epidish(beta.m = rexp, 
						   ref.m = sigm,
						   method = "RPC")$estF 
	pheat<-t(BloodFrac.m[,c(1:5)])
	rownames(pheat)<-c("CD14+ monocyte lineage", "CD19+ B-lymphocytes", 
						"CD4+ T-cells","CD56+ NK cells", "CD8+ cytotoxic T-lymphocytes")
	gd<-reshape2::melt(pheat)
	gd<-merge(gd,ip.infer,by.x="Var2",by.y="sample")
	save(gd,file=gsub(".png",".RData",filename))
	comp<-list(c("Cold","Hot"))
	
	gd.plot<-split(gd,gd$Var1)
	
	pp<-lapply(gd.plot,function(gd){
		max.v=max(gd$value)+sd(gd$value)*0.75
		p1<-ggplot(gd,aes(x=type,y=value,fill=type,color=type))+
		  geom_jitter(width = 0.15,size=0.5)+
		  geom_violin(width=0.35,color="black",alpha=.5)+
		  geom_boxplot(width=0.35,color="black",fill=NA,
					   outlier.alpha = 0,size=0.5)+
		  theme_few()+
		  theme(legend.position = "none")+
		  #ylab("Cell abundance estimated by MethylCIBERSORT")+
		  xlab("")+ylab(NULL)+
		  #facet_wrap(Var1~.,scales = "free",ncol=5)+
		  ylim(0,max.v)+
		  stat_compare_means(comparisons = comp,method = "wilcox.test",
							label="p.signif",size=5,
							method.args = list(alternative="less"))+
		  scale_color_manual(values = c(Cold="#6AB5DD",Hot="#EA614A"))+
		  scale_fill_manual(values = c(Cold="#6AB5DD",Hot="#EA614A"))
		  
		  
		p2<-ggplot(gd,aes(x=ipms.s,y=value,color=type))+
			  geom_point()+
			  geom_smooth(color="black",size=0.5)+
			  stat_cor(color="black",method="spearman",size=4.8)+
			  #facet_wrap(Var1~.,scales = "free",ncol=5)+
			  theme_few()+
			  scale_color_manual(values=c(Cold="#6AB5DD",Hot="#EA614A"))+
			  ylab(NULL)+
			  xlab("")+
			  theme(legend.position = "none")	
		layout <- 'AAABB'
		pp<-p2+p1+ plot_layout(design = layout)
	#	+ plot_annotation(title = unique(gd$Var1))
		ggsave(paste0(filename,unique(gd$Var1),".png"),
		height = 3,width = 5,dpi=600)	

	})

}


###############Fig 4B
load("GSE39279.non-knn.RData")
gd<-unique(gd[,c("Var2","ipms.s","type")])
sample<-read.delim("sample.txt",stringsAsFactors = F)
sample<-merge(sample,gd,by.x="sample",by.y="Var2")

my_comparisons=list(c("Cold","Hot"))
tz<-ggplot(sample,aes(y=tumor.size_cm,x=type,fill=type,color=type))+
  geom_jitter(width = 0.15,size=0.5)+
  geom_violin(width=0.35,color="black",alpha=.5)+
  geom_boxplot(width=0.35,color="black",fill=NA,
               outlier.alpha = 0,size=0.5)+
  theme_few()+
  scale_color_manual(values=c("#6AB5DD","#EA614A"))+
  scale_fill_manual(values=c("#6AB5DD","#EA614A"))+
  guides(fill=FALSE)+
  xlab(NULL)+
  ylab("Tumor size (cm)")+
  ylim(min(sample$tumor.size_cm,na.rm = T)-0.1,
       max(sample$tumor.size_cm,na.rm = T)+sd(sample$tumor.size_cm,na.rm = T)*0.75)+
  stat_compare_means(comparisons=my_comparisons,
                     label = "p.signif",
                     method = "wilcox.test",
                     method.args = list(alternative = "greater"))+
  theme(legend.position = "none")
p2<-ggplot(sample,aes(y=tumor.size_cm,x=ipms.s,color=type))+
  geom_point()+
  geom_smooth(color="black",size=0.5)+
  stat_cor(color="black",method="spearman",size=4.8)+
  #facet_wrap(Var1~.,scales = "free",ncol=5)+
  theme_few()+
  scale_color_manual(values=c(Cold="#6AB5DD",Hot="#EA614A"))+
  xlab("iPMS score")+
  ylab(NULL)+
  theme(legend.position = "none")	

library(patchwork)
layout<-'AABBB'
tz + p2+ plot_layout(design = layout)
ggsave("tumor_size.png",height=3.5,width=6,dpi = 600)


###############Fig 4C
sample$recurrence=gsub("recurrence: ","",sample$recurrence)
sample<-sample[-which(is.na(sample$adjuvant.chemotherapy)),]
sample$recurrence<-as.numeric(as.factor(sample$recurrence))-1
sample$Stage=as.numeric(factor(sample$Stage))
sample$gender=factor(sample$gender)
sample$Age=as.numeric(sample$Age)

res.cox <- coxph(Surv(adjuvant.chemotherapy, recurrence) ~ type+Stage+Age+gender, data=sample)
df <- with(sample,
           data.frame(type = sort(unique(type)), 
                      Age = rep(mean(Age),2),
                      gender = rep("male",2),
                      Stage = rep(mean(Stage, na.rm = TRUE),2))
           )


fit <- survfit(res.cox, newdata = df)
ggsurvplot(fit, data = df,
           pval = TRUE, conf.int = T, 
           risk.table = F,
           palette=c("#6AB5DD","#EA614A"), 
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_base()
)

ggsave("COX_cheomT_stage_gender_age_type.pdf",height = 5,width = 4.5,dpi=600)

ggforest(res.cox, data = sample)
ggsave("COX_cheomT_stage_gender_age_type_HR.pdf",height = 6,width = 6,dpi=600)

