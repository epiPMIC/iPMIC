######Fig 3B
load("RFMatrix.hot.cold.RData")
load("iPMS_sigMatrix.RData")
rownames(RFMatrix)<-RFMatrix$sample
RFMatrix<-cbind(sample_id=RFMatrix[,1],RFMatrix[,rownames(sigm)])

imc<-read.delim("R1/sample_info.txt",stringsAsFactors = F)
imc$IMC=gsub("IMC1|IMC2","Cold",imc$IMC)
imc$IMC=gsub("IMC3|IMC4|IMC5","Hot",imc$IMC)
ilf<-read.csv("R1/infiltration_estimation_for_tcga.csv",stringsAsFactors = F)
ilf<-ilf[,c(1,grep("TIMER",colnames(ilf)))]
imc<-merge(imc,ilf,by.x="sample",by.y="cell_type")

mm<-merge(imc,RFMatrix,by="sample_id")

ann_col=data.frame(row.names=mm$sample_id,
                   immunophenotype=mm$IMC)	

is<-split(mm$sample_id,ann_col$immunophenotype)
ann_colors = list(immunophenotype=c(Cold="#6AB5DD",Hot="#EA614A"))
rownames(mm)<-mm$sample_id
pheat<-mm[c(is$Cold,is$Hot),tail(colnames(mm),5)]
library(pheatmap)
pheatmap(t(pheat),
        filename="TCGA.iPMS.pdf",
        scale="none",
        annotation_col = ann_col,
        cluster_cols = F,
        cellheight = 15,
        cellwidth = 0.5,
        fontsize = 7,
        show_colnames = F,
        annotation_colors = ann_colors, 
        gaps_col = length(is$Cold),
        legend_breaks = c(0,0.3,0.5,0.7,1),
        color = colorRampPalette(c("#1B3463","#4E9CC7","#D1E8F2","white","#FCE3D2","#F4A784","#C03339","#720A21"))(150)
)

#####Fig 3C 
###python program
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
df=pd.read_csv('array.txt',sep="\t",index_col=0)
data = np.array(df)
ypos, xpos  = np.indices(data.shape) 
xpos = xpos.flatten()
ypos = ypos.flatten()
zpos = np.zeros(xpos.shape)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x_scale=30
y_scale=10
z_scale=12

scale=np.diag([x_scale, y_scale, z_scale, 1.0])
scale=scale*(1.0/scale.max())
scale[3,3]=1.0

def short_proj():
  return np.dot(Axes3D.get_proj(ax), scale)

ax.get_proj=short_proj
colors = pd.read_csv('color.txt',sep="\t",index_col=0)
colors = np.array(colors)
ax.bar3d(xpos,ypos,zpos, .3,.3,data.flatten(),color=colors.flatten())
ax.w_xaxis.set_ticks(xpos)
ax.w_yaxis.set_ticks(ypos)
plt.show()

########Fig3E-F
library(pROC)
roc.data<-data.frame(imc=imc.test$IMC,
                     hot=imc.test$Hot,
                     ipmsc=imc.test$ipms.s)
comp<-list(c("Cold","Hot"))
p1<-ggplot(roc.data,aes(x=imc,y=ipmsc,fill=imc,color=imc))+
  geom_jitter(width = 0.2,size=0.5)+
  geom_violin(width=0.4,color="black",alpha=.5)+
  geom_boxplot(width=0.4,color="black",fill=NA,
               outlier.alpha = 0,size=0.5)+
  theme_few()+
  theme(legend.position = "none")+
  ylab("iPMS score")+xlab(NULL)+
  stat_compare_means(comparisons = comp,method = "wilcox.test",label="p.signif")+
  scale_color_manual(values = c(Cold="#6AB5DD",Hot="#EA614A"))+
  scale_fill_manual(values = c(Cold="#6AB5DD",Hot="#EA614A"))+
  geom_hline(yintercept = 0,linetype="dashed",color="gray40")
ggsave("TCGA.test.iPMS.pdf",height = 4,width = 3,dpi=600)


