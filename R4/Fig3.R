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





