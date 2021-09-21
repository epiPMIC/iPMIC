#####Fig6B-D
library(patchwork)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(viridis)
library(survival)
library(survminer)

load("exp/TCGA.exp.RData")
comp<-list(c("Cold","Hot"))

p2<-ggplot(gd.exp,aes(x=type,y=value,color=type,fill=type))+
  geom_jitter(width = 0.2,size=0.5,alpha=0.5)+
  geom_violin(width=0.4,color="black",alpha=.5)+
  geom_boxplot(width=0.4,color="black",fill=NA,
               outlier.alpha = 0,size=0.5)+
  theme_few()+
  scale_fill_manual(values=c(Cold="#6AB5DD",Hot="#EA614A"))+
  scale_color_manual(values=c(Cold="#6AB5DD",Hot="#EA614A"))+
  stat_compare_means(comparisons = comp,
                     label = "p.signif",
                     method.args = list(alternative="less"))+
  facet_wrap(Var2~.,nrow=1)+
  xlab(NULL)+
  ylab("Gene expression levels")+
  theme(legend.position = "none")

load("TCGA_LUAD.knn.RData")
load("exp/TCGA.exp.RData")

load("TCGA_Cancer_sample.RDa")
sample<-sample[which(sample$type=="LUAD"),]
IMC<-read.delim("IMC/kmcluster.txt",stringsAsFactors = F,header = T)
IMC=IMC[sample$sample_id,]
IMC$sub5=gsub("C","IMC",IMC$sub5)
sample$IMC=IMC$sub5
sample$IMC=gsub("IMC1|IMC2","Cold",sample$IMC)
sample$IMC=gsub("IMC3|IMC4|IMC5","Hot",sample$IMC)
sample=sample[-which(is.na(sample$IMC)),]

gd<-merge(gd,gd.exp,by.x="Var2",by.y="Var1")
gd<-gd[grep("CD8",gd$Var1),]
gd<-merge(gd,sample[,c("sample_id","cyt")],by.x="Var2",by.y="sample_id")

p1<-ggplot(gd,aes(x=type.x,y=value.x,color=type.x,fill=type.x))+
  geom_jitter(width = 0.2,size=0.5,alpha=0.5)+
  geom_violin(width=0.4,color="black",alpha=.5)+
  geom_boxplot(width=0.4,color="black",fill=NA,
               outlier.alpha = 0,size=0.5)+
  theme_few()+
  scale_fill_manual(values=c(Cold="#6AB5DD",Hot="#EA614A"))+
  scale_color_manual(values=c(Cold="#6AB5DD",Hot="#EA614A"))+
  stat_compare_means(comparisons = comp,
                     label = "p.signif",
                     method.args = list(alternative="less"))+
  xlab(NULL)+
  ylab("Infiltration abundance by MethylCIBERSORT")+
  theme(legend.position = "none")

p4<-ggplot(gd,aes(x=type.x,y=cyt,color=type.x,fill=type.x))+
  geom_jitter(width = 0.2,size=0.5,alpha=0.5)+
  geom_violin(width=0.4,color="black",alpha=.5)+
  geom_boxplot(width=0.4,color="black",fill=NA,
               outlier.alpha = 0,size=0.5)+
  theme_few()+
  scale_fill_manual(values=c(Cold="#6AB5DD",Hot="#EA614A"))+
  scale_color_manual(values=c(Cold="#6AB5DD",Hot="#EA614A"))+
  stat_compare_means(comparisons = comp,
                     label = "p.signif",
                     method.args = list(alternative="less"))+
  xlab(NULL)+
  ylab("CYT score")+
  theme(legend.position = "none")


p3<-ggplot(gd,aes(x=value.x,y=value.y))+
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level..,fill = ..level..))+
  facet_grid(Var2.y~type.y,scales = "free")+
  scale_color_manual(values=c(Cold="#6AB5DD",Hot="#EA614A"))+
  geom_point(color="black",size=0.1,alpha=0.5)+
  theme_few()+
  xlim(-0.03,0.3)+
  ylim(-0.04,6)+
  stat_cor(method = "spearman")+
  xlab("CD8+ cytotoxic T-lymphocytes abundance")+
  ylab("Gene expression levels")+
  theme(legend.position = "none")+
  scale_fill_viridis(option ="D")

#######Fig 5E
data <- data.frame(
  group=c("Cold","Hot"),
  value=c(3,11)
)
ggplot(data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1,size=2, color="white") +
  coord_polar("y", start=0) +
  theme_void() 

#######Fig 5F
sample<-read.delim("methSample.txt",stringsAsFactors = F,header=F)
sv<-read.delim("sv.txt",stringsAsFactors = F,sep = "\t")
sample<-merge(sample,sv,by.x="V1",by.y="geo_accession")

load("GSE119144.ip.infer.RData") ##ip.infer
ipmsc<-ip.infer[[2]]
ip.infer<-ip.infer[[1]]
sample$ipmsc=ipmsc[sample$V1]
sample$ip.infer=ip.infer[sample$V1]
sample=sample[-which(sample$Clinical.benefit=="na"),]
sample$PFS=as.numeric(as.character(sample$PFS))
sample$PD_Event.1_Censoring.0=as.numeric(as.character(sample$PD_Event.1_Censoring.0))
fit <- survfit(Surv(PFS, PD_Event.1_Censoring.0) ~ ip.infer,
               data=sample[which(sample$Clinical.benefit=="NDB"),])  
res.cox <- coxph(Surv(PFS, PD_Event.1_Censoring.0) ~ ip.infer, 
                 data=sample[which(sample$Clinical.benefit=="NDB"),])
round(summary(res.cox)$conf.int[1],2)
ggsurvplot(fit,
           data=sample[which(sample$Clinical.benefit=="NDB"),],
           pval = T, conf.int = T, 
           palette=c("#6AB5DD","#EA614A"), 
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_few()
)
ggsave("NDB_PFS_survival.pdf", height = 5,width = 4.5,dpi=600)

res.cox <- coxph(Surv(PFS, PD_Event.1_Censoring.0) ~ ip.infer + Clinical.benefit, data=sample)
df <- with(sample,data.frame(ip.infer = c("Cold","Cold","Hot","Hot"),
                             Clinical.benefit=c("DCB","NDB","DCB","NDB")))
fit <- survfit(res.cox, newdata = df)

ggsurvplot(fit, data = df,
           pval = TRUE, conf.int = F, 
           risk.table = F, 
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_base()
)
ggsave("All_CB_PFS_survival.pdf", height = 5,width = 4.5,dpi=600)

#######Fig6G
sample<-read.delim("methSample.txt",stringsAsFactors = F,header=F)
sv<-read.delim("sv.txt",stringsAsFactors = F,sep = "\t")
sample<-merge(sample,sv,by.x="V1",by.y="geo_accession")
gd<-sample[,c("Clinical.benefit","type","riskscore")]
comp<-list(c("Cold","Hot"))
p1<-ggplot(gd,aes(x=type,y=riskscore,color=type,fill=type))+
  geom_jitter(width = 0.15,size=0.5)+
  geom_violin(width=0.35,color="black",alpha=.5)+
  geom_boxplot(width=0.35,color="black",fill=NA,
               outlier.alpha = 0,size=0.5)+
  theme_few()+
  xlab(NULL)+ylab("iRS")+
  theme(legend.position = "none")+
  stat_compare_means(comparisons = comp,method = "wilcox.test",
                     method.args = list(alternative="greater"))+
  scale_color_manual(values = c(Cold="#6AB5DD",Hot="#EA614A"))+
  scale_fill_manual(values = c(Cold="#6AB5DD",Hot="#EA614A"))

comp<-list(c("NDB","DCB"))
gd$Clinical.benefit=factor(gd$Clinical.benefit,levels = c("NDB","DCB"))
p2<-ggplot(gd,aes(x=Clinical.benefit,y=riskscore,color=Clinical.benefit,fill=Clinical.benefit))+
  geom_jitter(width = 0.15,size=0.5)+
  geom_violin(width=0.35,color="black",alpha=.5)+
  geom_boxplot(width=0.35,color="black",fill=NA,
               outlier.alpha = 0,size=0.5)+
  theme_few()+
  xlab(NULL)+ylab("")+
  theme(legend.position = "none")+
  stat_compare_means(comparisons = comp,method = "wilcox.test",
                     method.args = list(alternative="greater"))+
  scale_color_manual(values = c("#6AB5DD","#EA614A"))+
  scale_fill_manual(values = c("#6AB5DD","#EA614A"))
library(patchwork)
p1+p2
ggsave("iRS.pdf",height = 4,width = 6)

#####Fig6H
sample$PFS=as.numeric(as.character(sample$PFS))
sample$PD_Event.1_Censoring.0=as.numeric(as.character(sample$PD_Event.1_Censoring.0))
sample<-merge(sample,RFM,by.x="V1",by.y="sample")
res.cox <- coxph(Surv(PFS, PD_Event.1_Censoring.0) ~ cg10122865+cg15138339+cg00249511+cg18103859,data=sample)
sample$riskscore<-4.2235*sample$cg10122865+ 1.3652*sample$cg15138339-0.6366*sample$cg00249511-2.5145*sample$cg18103859
save(sample,file="Sample_risk.RData")

setwd("GSE119144_GSE135222/")
load("Sample_risk.RData")
col<-pal_lancet()(6)
names(col)=c("Mutation_burden","Aneuploidy","cg10122865","cg15138339","cg00249511","cg18103859")
pdf("roc.immbenefit.pdf",height = 4,width = 4)
roc(sample$Clinical.benefit,sample$riskscore,
    print.thres=F, percent=F,
    ci=TRUE, boot.n=100, ci.alpha=0.95, stratified=FALSE,
    plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE)
for (i in c("Mutation_burden","Aneuploidy","cg10122865","cg15138339","cg00249511","cg18103859")) {
  plot.roc(sample$Clinical.benefit,sample[,i], 
           col=col[i], ci=TRUE,
           lwd=1, print.auc=T, print.auc.cex=0.5, 
           add=TRUE, print.auc.y=0.45)  
}
dev.off()

















