##########Fig5A

library(maftools)

load("DE_ICP.RData")
icp<-icp[which(abs(icp$logFC.normal) > 0.5 &icp$q.normal<0.01),] 
icp<-icp[which(abs(icp$logFC) > 0.5 & icp$q < 0.01),] 

icp.cold.maf = read.maf(maf ="LUAD_icp_Cold.maf")
oncoplot(maf = icp.cold.maf, top = 20,sortByAnnotation = TRUE)

icp.hot.maf = read.maf(maf ="LUAD_icp_Hot.maf")
oncoplot(maf = icp.hot.maf, top = 20,sortByAnnotation = TRUE)

