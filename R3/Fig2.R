library(ggplot2)
library(ggthemes)
library(ggsci)

##########Figure 2A
bed<-read.delim("hg19.HM450K.bed",stringsAsFactors = F,header=F)
dmc_hot<-read.delim("DMC_Hot.txt",stringsAsFactors = F)
dmc_cold<-read.delim("DMC_Cold.txt",stringsAsFactors = F)
const<-intersect(dmc_cold$probe,dmc_hot$probe)
dmc<-unique(rbind(dmc_cold,dmc_hot))
rownames(dmc)<-dmc$probe
group=rep(NA,length(dmc$probe))
names(group)=dmc$probe
group[dmc_hot$probe]<-"Hot"
group[dmc_cold$probe]<-"Cold"
group[const]<-"Constitutive"
dmc$group=group
bed<-merge(bed,dmc[,c("probe","derta_ac","derta_at")],by.x = "V4",by.y="probe")

load("Density.RData")
get.gd<-function(x,y){
  x<-data.frame(t(x))
  x$probe=rownames(x)
  colnames(x)=as.character(as.matrix(x[1,]))
  x<-x[-1,]
  x$type=y
  return(x)
}
gd.hot<-get.gd(x=AcM,y="Hot")
gd.cold<-get.gd(x=AtM,y="Cold")
gd<-rbind(gd.hot,gd.cold)
gd$type[which(gd$Group.1 %in% intersect(gd.cold$Group.1,gd.hot$Group.1))]="Constitutive"
gd<-unique(gd)
gd<-reshape::melt(gd,id=c("Group.1","type"))
gd$value=as.numeric(as.character(gd$value))
ggplot(gd,aes(x=value,color=variable))+
  geom_line(stat="density")+
  facet_grid(type~.)
  
##########Figure 2B
exp<-read.delim("TCGA-LUAD.htseq_fpkm.tsv",header = T,row.names = 1)
cs<-grep("[.]0..$",colnames(exp))
ns<-grep("[.]1..$",colnames(exp))
imc<-read.delim("sample_info.txt",stringsAsFactors = F)
imp<-split(imc$sample_id,imc$IMC)
deg<-function(cs,ns){
  p<-apply(exp,1,function(x){
    w<-wilcox.test(x[cs],x[ns])
    w$p.value
  })
  q=p.adjust(p,method = "BH")
  FC=log2(rowMeans(exp[,cs])/rowMeans(exp[,ns]))
  DEG<-data.frame(gene_id=rownames(exp),MC=rowMeans(exp[,cs]),MN=rowMeans(exp[,ns]),p=p,q=q,log2FC=FC)
  return(DEG)
}
deg.hot<-deg(cs=c(imp$IMC3,imp$IMC4,imp$IMC5),ns)
colnames(deg.hot)<-paste0(colnames(deg.hot),".hot")
deg.cold<-deg(cs=c(imp$IMC1,imp$IMC2),ns)
colnames(deg.cold)<-paste0(colnames(deg.cold),".cold")
DEG<-deg(cs,ns)
DEG<-cbind(DEG,cbind(deg.hot,deg.cold))
DEG$gene_id=gsub("[.].*","",DEG$gene_id)

DMG<-merge(gd,merge(DMC,DEG,by="gene_id"),by.x="Group.1",by.y = "probe")
DMG$delta.hot=as.numeric(as.character(DMG$Active))-as.numeric(as.character(DMG$normal))
DMG$delta.cold=as.numeric(as.character(DMG$Attacked))-as.numeric(as.character(DMG$normal))
DMG$delta=(as.numeric(as.character(DMG$Active))+as.numeric(as.character(DMG$Attacked)))/2-as.numeric(as.character(DMG$normal))

ggplot(DMG,aes(x=delta,y=log2FC,color=type))+
  geom_point(alpha=0.5,size=1)+
  theme_few()+
  xlim(-0.4,0.4)+
  ylim(-5,8)+
  xlab("Methylation delta, Mean(cancer)-Mean(normal)")+
  scale_color_manual(values = c(Cold="#8854D0",Hot="#F7B731",Constitutive="#FA8231"))+
  geom_hline(yintercept = 0,linetype="dashed",color="gray60")+
  geom_vline(xintercept = 0,linetype="dashed",color="gray60")+
  theme(legend.position = "none")
ggsave("DMC_DEG.png",height = 5,width = 5)

######Figure 2C
cg<-read.delim("cancermine_collated.tsv",header=T,stringsAsFactors = F)
cg<-cg[grep("lung",cg$cancer_normalized),]
pcg<-read.delim("Pcg_id.txt",header=F,stringsAsFactors = F)
dg<-read.delim("DMC_TargetGene.txt",header=T,stringsAsFactors = F)
dg<-merge(dg,pcg,by.x="gene_id",by.y="V1")
gd<-data.frame(type=factor(c("CancerMine","DisGeNET","CGC Hallmark"),
                           levels = c("CancerMine","DisGeNET","CGC Hallmark")),
               p=c(phyper(55,1403,21000-1403,403,lower.tail = F),
                   phyper(14,362,21000-362,403,lower.tail = F),
                   phyper(16,294,21000-294,403,lower.tail = F)),
               num=c(55,14,16))
gd$logP=-log10(gd$p)
ggplot(gd,aes(x=type,y=logP,size=num))+
  geom_point(alpha=0.5)+
  scale_size_area(max_size = 18)+
  ylim(0,8)+
  xlab(NULL)+
  theme_classic()+
  geom_hline(yintercept = -log10(0.05),linetype="dashed",color="gray80")
ggsave("TumorGene_enrich.pdf",height = 3,width = 5)

######Figure 2D
hallmark<-read.csv("seTG_Hallmark.csv",stringsAsFactors = F)
dt<-data.frame(table(hallmark$hallmark))
colnames(dt)<-c("hallmark","Count")
dt$hallmark=factor(dt$hallmark,levels = levels(dt$hallmark)[c(9,8,3,5,6,1,7,10,2,4)])
library(ggplot2)
ggplot(dt, aes(x = hallmark, y = Count, fill = hallmark)) + 
  geom_bar(stat = "identity") + 
  coord_polar() + 
  theme_bw() + 
  labs(x = "", y = "") + 
  theme(axis.text.y = element_blank()) +  
  theme(axis.ticks = element_blank()) +  
  theme(panel.border = element_blank()) +   
  theme(legend.position = "none") +
  theme(title = element_text(vjust = -56, face = "bold")) +
  scale_fill_manual(values = c(`suppression of growth`="#563236",
                               `escaping immune response to cancer`="#a2002f",
                               `cell replicative immortality`="#04656c",
                               `tumour promoting inflammation`="#c95202",
                               `invasion and metastasis`="#262626",
                               `angiogenesis`="#ff0202",
                               `genome instability and mutations`="#022b59",
                               `escaping programmed cell death`="#676767",
                               `change of cellular energetics`="#47055c",
                               `proliferative signalling`="#114009"))
ggsave("hallmark.pdf",height = 3,width = 3)

######Figure 2E
gdd<-read.delim("FuncEnrich.txt",stringsAsFactors = F)
p<-ggplot(gdd,aes(x=Description,y=NodeA,size=-LogP))+
  geom_point(shape=15,alpha=0.5)+
  geom_point(shape=0,color="black")+
  theme_pander()+
  scale_size_area(max_size = 10)+
  scale_color_gradient(low = "#ffffff",high = "#326280")+
  xlab(NULL)+ylab(NULL)+
  theme(axis.text.x = element_text(vjust = 0.5, 
                                   hjust = 1, angle = 45))
ggsave("function_enrich.pdf",height = 4,width = 8.5)
