library(ggplot2)
library(ggthemes)
library(viridis)
library(pheatmap)
library(maftools)
library(ggpubr)
library(EpiDISH)
library(rpart)
library(visNetwork)
library(rpart.plot)
library(pROC)
library(survival)
library(survminer)

###########Figure 5A
imc<-read.delim("sample_info.txt",stringsAsFactors = F)
imc$IMC=gsub("IMC1|IMC2","Cold",imc$IMC)
imc$IMC=gsub("IMC3|IMC4|IMC5","Hot",imc$IMC)
m<-read.csv("TCGA.MeImm.Meth.txt",stringsAsFactors = F,header=T,row.names = 1)
rexp<-m[,imc$sample_id]
load("iPMS_sigMatrix.RData") ##sigm
ip.infer <- epidish(beta.m = rexp, 
                    ref.m = sigm,
                    method = "RPC")$estF 
ipmsc<-ip.infer[,2]-ip.infer[,1]
names(ipmsc)<-rownames(ip.infer)
ip.infer <- c("Cold","Hot")[apply(ip.infer, 1, which.max)]
names(ip.infer)<-colnames(rexp)
imc$ipmsc<-ipmsc[imc$sample_id]
imc$ip.infer<-ip.infer[imc$sample_id]
exp<-read.delim("TCGA-LUAD.htseq_fpkm.tsv",stringsAsFactors = F)
exp$Ensembl_ID<-gsub("[.].*$","",exp$Ensembl_ID)
icp<-read.delim("final_ICP.txt",stringsAsFactors = F,header = F)
icp$V3=gsub(" ","",icp$V3)
icp.exp<-exp[which(exp$Ensembl_ID %in% icp$V3),]
icp.exp<-merge(icp[,c("V1","V3")],icp.exp,by.x="V3",by.y="Ensembl_ID")
rownames(icp.exp)<-icp.exp$V1
icp.exp<-icp.exp[,-2:-1]

pheat=aggregate(icp.exp,by=list(gene=rownames(icp.exp)),mean)
rownames(pheat)<-pheat$gene
pheat<-pheat[,-1]
sig<-apply(pheat, 1, function(x){
  w<-wilcox.test(x[sample$Hot],x[sample$Cold])
  lab<-c(0.0001,0.001,0.05,1,10)
  names(lab)<-c("***","**","*","NS","NS")
  return(names(which(w$p.value<lab))[1])
})
rownames(pheat)=paste0(rownames(pheat)," (",sig,")")

icp.cor<-apply(icp.exp, 1, function(x){
  c.hot<-cor.test(x[sample$Hot],imc[sample$Hot,"ipmsc"],metohd="spearman")
  c.cold<-cor.test(x[sample$Cold],imc[sample$Cold,"ipmsc"],metohd="spearman")
  rbind(data.frame(rho=c.hot$estimate,p=c.hot$p.value,type="Hot"),
        data.frame(rho=c.cold$estimate,p=c.cold$p.value,type="Cold"))
})
gd<-do.call(rbind,icp.cor)
gd$icp<-gsub("[.].*$","",rownames(gd))

pheat.m<-as.matrix(data.frame(Hot=rowMeans(pheat[,sample$Hot]),
                              Cold=rowMeans(pheat[,sample$Cold]),
                              Normal=rowMeans(pheat[,sample$normal])))
p1<-pheatmap(scale(pheat.m),
            filename="pheat.icp.png",
            scale="row",
            cluster_cols = F,
            clustering_method = "average",
            cellheight = 8,
            cellwidth = 8,
            fontsize = 6,
            show_colnames = T,
            annotation_colors = ann_colors, 
            gaps_col = 1:3,
            color = colorRampPalette(c("#1B3463","#4E9CC7","#D1E8F2","white","#FCE3D2","#F4A784","#C03339","#720A21"))(150)
)
ord<-gsub(" [(].*","",p1$tree_row$labels[p1$tree_row$order])
gd$icp=factor(gd$icp,levels = rev(ord))
p<-ggplot(gd,aes(y=icp,x=type,color=rho,size=-log10(p)))+
  geom_point(alpha=0.5)+
  scale_size_area(max_size = 5)+
  theme_pander()+
  xlab(NULL)+
  ylab(NULL)+
  #theme(legend.position = "bottom")+
  scale_color_viridis(option ="A")

p1<-as.ggplot(p1)
cowplot::plot_grid(p1, p, ncol=2)
ggsave("TCGA_icp.png",height = 20,width = 6,dpi=600)

###########Figure 5B
icp_maf = read.maf(maf ="LUAD_icp_merge.maf",clinicalData = imc)
plotmafSummary(maf = icp_maf, rmOutlier = TRUE, addStat = 'median')

icp.hot.maf = read.maf(maf ="LUAD_icp_Hot.maf",clinicalData = imc)
icp.cold.maf = read.maf(maf ="LUAD_icp_Cold.maf",clinicalData = imc)
p.hot<-oncoplot(maf = icp.hot.maf, top = 30,
         clinicalFeatures="IMC",sortByAnnotation = TRUE)
p.cold<-oncoplot(maf = icp.cold.maf, top = 30,
         clinicalFeatures="IMC",sortByAnnotation = TRUE)

###########Figure 5C and Figure S5

sp.icp<-c("LILRB5","CD200R1","CD96","EZH2","BTNL2",
  "CD44","HAVCR1","TLR8","CD226","BTN3A2",
  "SPAG9","MERTK","CD229","CD244","ITGB2",
  "BTK","PP2A","EDNRB")

gd<-lapply(sp.icp, function(x){
  data.frame(sample=colnames(icp.exp),
             value=as.numeric(as.matrix(icp.exp[x,])),
             sp.icp=x,
             pdl1=as.numeric(as.matrix(icp.exp["PD-L1",])))
})
gd<-do.call(rbind,gd)
gd<-merge(gd,imc[,c("sample_id","IMC")],by.x="sample",by.y="sample_id")
library(ggplot2)
library(ggpubr)
library(ggthemes)

gd$sp.icp=gsub("PP2A","PPP2R1A",gd$sp.icp)
gd$sp.icp=gsub("CD229","LY9",gd$sp.icp)
gd$sp.icp=factor(gd$sp.icp,levels = c("LILRB5","CD200R1","CD96","EZH2","BTNL2",
                                      "CD44","HAVCR1","TLR8","CD226","BTN3A2",
                                      "SPAG9","MERTK","LY9","CD244","ITGB2",
                                      "BTK","PPP2R1A","EDNRB"))

p<-ggplot(gd,aes(x=pdl1,y=value,color=IMC))+
  geom_point(alpha=0.5)+
  geom_smooth(size=0.5)+
  stat_cor(method = "spearman")+
  theme_few()+
  scale_color_manual(values=c(Cold="#6AB5DD",Hot="#EA614A"))+
  ylab("Expression levels (FPKM)")+
  xlab("PD-L1 expression levels (FPKM)")+
  theme(legend.position = "none")+
  facet_wrap(sp.icp~.,scales = "free",ncol=9)
ggsave("PDL1.co.exp.pdf",height = 5,width = 19,dpi = 600)


gd<-lapply(sp.icp, function(x){
  data.frame(sample=colnames(icp.exp),
             value=as.numeric(as.matrix(icp.exp[x,])),
             sp.icp=x,
             pdl1=as.numeric(as.matrix(icp.exp["PD-1",])))
})
gd<-do.call(rbind,gd)
gd<-merge(gd,imc[,c("sample_id","IMC")],by.x="sample",by.y="sample_id")

gd$sp.icp=gsub("PP2A","PPP2R1A",gd$sp.icp)
gd$sp.icp=gsub("CD229","LY9",gd$sp.icp)
gd$sp.icp=factor(gd$sp.icp,levels = c("LILRB5","CD200R1","CD96","EZH2","BTNL2",
                                      "CD44","HAVCR1","TLR8","CD226","BTN3A2",
                                      "SPAG9","MERTK","LY9","CD244","ITGB2",
                                      "BTK","PPP2R1A","EDNRB"))

p<-ggplot(gd,aes(x=pdl1,y=value,color=IMC))+
  geom_point(alpha=0.5)+
  geom_smooth(size=0.5)+
  stat_cor(method = "spearman")+
  theme_few()+
  scale_color_manual(values=c(Cold="#6AB5DD",Hot="#EA614A"))+
  ylab("Expression levels (FPKM)")+
  xlab("PD-1 expression levels (FPKM)")+
  theme(legend.position = "none")+
  facet_wrap(sp.icp~.,scales = "free",ncol=4)
ggsave("PD1.co.exp.png",height = 14,width = 10,dpi = 600)

###########Figure 5D
sample<-read.delim("methSample.txt",stringsAsFactors = F,header=F)
sv<-read.delim("sample.txt",stringsAsFactors = F,sep = "\t")
sample<-merge(sample,sv,by.x="V1",by.y="geo_accession")

load("GSE126045.ip.infer.RData") ##ip.infer
ipmsc<-ip.infer[[2]]
ip.infer<-ip.infer[[1]]
sample$ipmsc=ipmsc[sample$V1]
sample$ip.infer=ip.infer[sample$V1]

load("GSE126045.knn.RData")
sigm<-as.matrix(read.delim("lung_NSCLC_not.specified_v2_Signature.txt.txt",row.names=1))
BloodFrac.m <- epidish(beta.m = rexp, 
                       ref.m = sigm,
                       method = "RPC")$estF 
pheat<-t(BloodFrac.m)
rownames(pheat)<-c("CD14+ monocyte lineage", "CD19+ B-lymphocytes", 
                   "CD4+ T-cells","CD56+ NK cells", "CD8+ cytotoxic T-lymphocytes")
pheat<-pheat[,names(sort(ip.infer))]
rownames(sample)<-sample$V1
ann_col=data.frame(row.names=names(ip.infer),
                   Responsiveness = sample[names(ip.infer),"Responsiveness"],
                   progression = sample[names(ip.infer),"Best.response"],
                   immunophenotype=ip.infer)
is<-split(names(ip.infer),ann_col$immunophenotype)
sig<-apply(pheat, 1, function(x){
  w<-wilcox.test(x[is$Hot],x[is$Cold],alternative = "greater")
  lab<-c(0.0001,0.001,0.05,1,10)
  names(lab)<-c("***","**","*","NS","NS")
  return(names(which(w$p.value<lab))[1])
})
rownames(pheat)=paste0(gsub("_TIMER","",rownames(pheat)),"(",sig,")")
ann_colors = list(immunophenotype=c(Cold="#6AB5DD",Hot="#EA614A"))
pheatmap(scale(pheat),
         filename="immune.pdf",
         scale="row",
         annotation_col = ann_col,
         cluster_cols = F,
         cellheight = 9,
         cellwidth = 5,
         fontsize = 7,
         show_colnames = F,
         annotation_colors = ann_colors, 
         gaps_col = length(is$Cold),
         color = colorRampPalette(c("#1B3463","#4E9CC7","#D1E8F2","white","#FCE3D2","#F4A784","#C03339","#720A21"))(150)
)

###########Figure 5E
load("iPMS_sigMatrix.RData")
ipm<-rexp[rownames(sigm),sample$V1]
ipm<-ipm[-which(is.na(ipm$GSM3589655)),]
ipm<-data.frame(t(ipm))
ipm$V1=rownames(ipm)
pm<-merge(sample,ipm,by="V1")
gexp<-read.delim("GSE126044_counts.txt",stringsAsFactors = F)
gexp<-gexp[which(gexp$X %in% c("LILRB5","CD200R1","CD96","EZH2","BTNL2",
                               "CD44","HAVCR1","TLR8","CD226","BTN3A2",
                               "SPAG9","MERTK","LY9","CD244","ITGB2",
                               "BTK","PP2R1A","EDNRB")),]

rownames(gexp)<-gexp$X
gexp<-gexp[,-1]
ipm<-data.frame(t(gexp))
ipm$V1=rownames(ipm)
pm<-merge(pm,ipm,by.x="V2",by.y="V1")
pm<-pm[,c(8,4,5,6,7,12,13:71)]

m.rpart <- rpart(Best.response ~ .-Best.response,
                 data = pm[,c(-4:-1,-7)], 
                 minsplit=1,maxdepth=30)

rpart.plot(m.rpart)

visTree(m.rpart,main = "title",
        fallenLeaves = T,
        simplifyRules=F,
        colorEdges="#d1d8e0",nodesFontSize = 20,
        edgesFontSize = 18,
        colorY = c("blue","red","green"))

###########Figure 5F
sample<-read.delim("methSample.txt",stringsAsFactors = F,header=F)
sv<-read.delim("sv.txt",stringsAsFactors = F,sep = "\t")
sample<-merge(sample,sv,by.x="V1",by.y="geo_accession")
load("GSE119144.knn.RData")
load("E:/MeImm/LUAD/0.Comb/R4/iPMS_sigMatrix.RData")
ipm<-rexp[rownames(sigm),sample$V1]
ipm<-ipm[-which(is.na(ipm$GSM3902399)),]
ipm<-data.frame(t(ipm))
ipm$V1=rownames(ipm)
pm<-merge(sample,ipm,by="V1")
gexp<-read.csv("GSE135222_exp.csv",row.names = 1)
icp<-read.delim("../ICInet/ICP_gene_symbol.txt",header = F,stringsAsFactors = F)
gexp<-gexp[c("LILRB5","CD200R1","CD96","EZH2","BTNL2",
             "CD44","HAVCR1","TLR8","CD226","BTN3A2",
             "SPAG9","MERTK","LY9","CD244","ITGB2",
             "BTK","PPP2R1A","EDNRB","CD274","PDCD1"),]
ipm<-data.frame(t(gexp))
ipm$V1=rownames(ipm)
pm<-merge(pm,ipm,by.x="V2",by.y="V1")
load("GSE119144.ip.infer.RData") ##ip.infer
ipmsc<-ip.infer[[2]]
ip.infer<-ip.infer[[1]]
sample$ipmsc=ipmsc[sample$V1]
sample$ip.infer=ip.infer[sample$V1]
sample=sample[-which(sample$Clinical.benefit=="na"),]

col<-pal_npg()(4)
names(col)=c("CD244","cg06436185","cg06772221","CD274")

pdf("roc.immbenefit.pdf",height = 4,width = 4)
roc(pm$Clinical.benefit,pm$CD244+pm$cg06436185+pm$cg06772221+pm$CD274, 
    print.thres=F, percent=F,
    ci=TRUE, boot.n=100, ci.alpha=0.95, stratified=FALSE,
    plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
    print.auc=TRUE, show.thres=TRUE)
for (i in names(col)) {
  plot.roc(pm$Clinical.benefit,pm[,i], 
           col=col[i], 
           lwd=1, print.auc=T, print.auc.cex=0.5, 
           add=TRUE, print.auc.y=0.45)  
}
dev.off()

###########Figure 5G
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




