library(ggplot2)
library(ggsci)
library(ggthemes)
library(ggpubr)
library(umap)
library(ggradar)
library(dplyr)
library(scales)
library(tibble)

##########Figure 1A
load("TCGA_Cancer_sample.RDa")
sample<-sample[which(sample$type=="LUAD"),]
IMC<-read.delim("IMC/kmcluster.txt",stringsAsFactors = F,header = T)
IMC=IMC[sample$sample_id,]
IMC$sub5=gsub("C","IMC",IMC$sub5)
sample$IMC=IMC$sub5
sample=sample[-which(is.na(sample$IMC)),]
write.table(sample,"sample_info.txt",quote = F,col.names = T,row.names = F,sep = "\t")

udat<-sample[,35:39]
rownames(udat)=sample$sample_id
umap<-umap(udat,n_neighbors = 60, min_dist = 0.001, target_weight = 0.5)
umap<-umap$layout
umap=data.frame(sample=rownames(udat),umap)
colnames(umap)<-c("sample","UMAP1","UMAP2")
sample$UMAP1=umap$UMAP1
sample$UMAP2=umap$UMAP2
sample$pheno=gsub("IMC1|IMC2","Cold",sample$IMC)
sample$pheno=gsub("IMC3|IMC4|IMC5","Hot",sample$pheno)

gd<-data.frame(sample[,c("UMAP1","UMAP2")],
               IMC=sample$IMC,
               pheno=sample$pheno)
p1<-ggplot(gd,aes(x=UMAP1,y=UMAP2,color=IMC))+
  geom_point(size=1.2)+
  theme_few()+
  scale_color_d3(palette = "category20")
p2<-ggplot(gd,aes(x=UMAP1,y=UMAP2,color=pheno))+
    geom_point(size=1.2)+
    theme_few()+
    scale_color_manual(values = c(Cold="#6AB5DD",Hot="#EA614A"))
p1+p2+ plot_layout(ncol = 1)


##########Figure 1B
ic<-read.csv("infiltration_estimation_for_tcga.csv",stringsAsFactors = F)
ic<-ic[,c(1,grep("TIMER",colnames(ic)))]
ic<-merge(sample,ic,by.x="sample",by.y="cell_type")
gd<-reshape2::melt(ic[,c("IMC",colnames(ic)[grep("TIMER",colnames(ic))])])
ggplot(gd,aes(x=variable,y=value,color=variable,fill=variable))+
  geom_jitter(width = 0.3,size=0.5)+
  geom_boxplot(width=0.3,color="black",
               outlier.alpha = 0,size=0.5)+
  theme_hc()+
  scale_fill_lancet()+
  scale_color_lancet()+
  ylab("Infiltration levels")+xlab(NULL)+
  facet_wrap(IMC~.,nrow=1)
ggsave("immune.pdf",height = 4,width = 14,dpi=600)

aa<-aggregate(sample[,31:41],by=list(IMC=sample$IMC),mean)
d1<-aa[,c(1,4,5,11,12)]
colnames(d1)=c("IMC","ESTIMATE_score","Immune_score","CYT_score","HLA_score")
radar <- d1[,-1] %>% as_tibble(rownames = "IMC") %>% 
  mutate_at(vars(-IMC), rescale) 
ggradar(radar,group.colours = pal_d3()(5),legend.position = "none")
ggsave("E:/MeImm/LUAD/0.Comb/R1/Immune_score.pdf",width = 5,height = 5)

d2<-aa[,c(1,6:10)]
radar <- d2[,-1] %>% as_tibble(rownames = "IMC") %>% 
  mutate_at(vars(-IMC), rescale) 
ggradar(radar,group.colours = pal_d3()(5),legend.position = "none")
ggsave("Pathway_score.pdf",width = 5,height = 5)

##########Figure 1C-F
meth<-read.delim("LUAD_MethLevel.txt",header = T,stringsAsFactors = F)
meth$sample=gsub("-",".",meth$sample)
meth$meth=as.numeric(meth$meth)
colnames(meth)=c("sample_id","GMeth")
sample<-merge(sample,meth,by="sample_id")
dmeth<-read.delim("LUAD_DMethLevel.txt",header = T)
dmeth$sample=gsub("[.]","-",dmeth$sample)

p1<-ggplot(sample,aes(x=IMC,y=TMB,fill=IMC,color=IMC))+
  geom_jitter(width = 0.2,size=0.5)+
  geom_boxplot(width=0.2,color="black",
               outlier.alpha = 0,size=0.5)+
  theme_few()+
  scale_fill_d3()+
  scale_color_d3()+
  theme(legend.position = "none")+
  ylab("TMB (/Mb)")+xlab(NULL)+
  stat_compare_means(method = "kruskal.test")+
  geom_hline(yintercept = 10,linetype="dashed",color="gray40")

p2<-ggplot(sample,aes(x=IMC,y=Aneuploidy.Score,fill=IMC,color=IMC))+
  geom_jitter(width = 0.2,size=0.5)+
  geom_boxplot(width=0.2,color="black",
               outlier.alpha = 0,size=0.5)+
  theme_few()+
  scale_fill_d3()+
  scale_color_d3()+
  theme(legend.position = "none")+
  ylab("Aneuploidy")+xlab(NULL)+
  stat_compare_means(method = "kruskal.test")

p3<-ggplot(sample,aes(x=IMC,y=GMeth,fill=IMC,color=IMC))+
  geom_jitter(width = 0.2,size=0.5)+
  geom_boxplot(width=0.2,color="black",
               outlier.alpha = 0,size=0.5)+
  theme_few()+
  scale_fill_d3()+
  scale_color_d3()+
  theme(legend.position = "none")+
  ylab("Methylation levels (β-value)")+xlab(NULL)+
  stat_compare_means(method = "kruskal.test")
  
library(patchwork)
p1+p2+p3

##########Figure 1G-H, FigureS2B-D
load("TCGA_Cancer_sample.RDa")
sample<-sample[which(sample$type=="LUAD"),]
IMC<-read.delim("IMC/kmcluster.txt",stringsAsFactors = F,header = T)
IMC=IMC[sample$sample_id,]
IMC$sub5=gsub("C","IMC",IMC$sub5)
sample$IMC=IMC$sub5
sample=sample[-which(is.na(sample$IMC)),]
colnames(sample)[11]="stage"
sv<-sample[,c("IMC","stage","Diagnosis.Age","Sex",
              "Overall.Survival..Months.","Overall.Survival.Status",
              "Progress.Free.Survival..Months.",
              "Progression.Free.Status","Radiation.Therapy")]
lapply(c("Overall.Survival","Progress"), function(x){
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
           
           res.cox <- coxph(Surv(time, status) ~ IMC+stage+Sex+Diagnosis.Age, data=ssv)
           df <- with(ssv,
                      data.frame(IMC = sort(unique(ssv$IMC)), 
                                 Diagnosis.Age = rep(mean(ssv$Diagnosis.Age, na.rm = TRUE),length(unique(ssv$IMC))),
                                 Sex = rep("Male",length(unique(ssv$IMC))),
                                 stage = rep("STAGE I",length(unique(ssv$IMC)))
                      )
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

           })
})



##########Figure S2A
gd<-data.frame(stage=as.factor(sample$Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code),
               IMC=sample$IMC)
p1<-ggplot(gd,aes(x=IMC,y=as.numeric(stage),fill=IMC,color=IMC))+
  geom_jitter(width = 0.2,size=0.5)+
  geom_boxplot(width=0.2,color="black",
               outlier.alpha = 0,size=0.5)+
  theme_few()+
  scale_fill_d3()+
  scale_color_d3()+
  theme(legend.position = "none")+
  ylab("Disease Stage")+xlab(NULL)+
  stat_compare_means(method = "kruskal.test")+
  scale_y_continuous(breaks=1:9,labels = gsub("STAGE","",levels(gd$stage)))

ggsave("Stage.png",height = 2.5,width = 4,dpi=600)


