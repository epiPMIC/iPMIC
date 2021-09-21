##########Fig5A
library(survival)
library(survminer)
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(viridis)
library(patchwork)
library(maftools)

load("DE_ICP.RData")
icp<-icp[which(abs(icp$logFC.normal) > 0.5 &icp$q.normal<0.01),] 
icp<-icp[which(abs(icp$logFC) > 0.5 & icp$q < 0.01),] 

icp.cold.maf = read.maf(maf ="LUAD_icp_Cold.maf")
oncoplot(maf = icp.cold.maf, top = 20,sortByAnnotation = TRUE)

icp.hot.maf = read.maf(maf ="LUAD_icp_Hot.maf")
oncoplot(maf = icp.hot.maf, top = 20,sortByAnnotation = TRUE)

##########Fig5B-C
risk.iPMS<-lapply(c("Progress"), function(x){
  ssv<-icp.exp[,c(1:4,grep(x,colnames(icp.exp)),14:46)]
  colnames(ssv)[c(5:6)]<-c("time","status")           
  loc<-union(which(is.na(ssv[,5])),
             which(is.na(ssv[,6])))
  if(length(loc)>0){ssv<-ssv[-loc,]}
  ssv$status=as.numeric(gsub(":.*","",ssv$status))
  ssv$time=as.numeric(as.character(ssv$time))

  res.cox <- coxph(Surv(time, status) ~ cg10122865+cg15138339+cg24887265+
                     cg00249511+cg18103859,data=ssv)
  riskScore<-res.cox$coefficients[1]*ssv[,names(res.cox$coefficients[1])]+
    res.cox$coefficients[2]*ssv[,names(res.cox$coefficients[2])]+
    res.cox$coefficients[3]*ssv[,names(res.cox$coefficients[3])]+
    res.cox$coefficients[4]*ssv[,names(res.cox$coefficients[4])]+
    res.cox$coefficients[5]*ssv[,names(res.cox$coefficients[5])]
  
  group<-function(x,cutoff){
    group<-rep(NA,length(x))
    group[which(x>=cutoff)]<-1
    group[which(x<cutoff)]<-0
    group
  }
  rg<-group(riskScore,cutoff = median(riskScore))
  ssv$rg=rg
  
  fit <- survfit(Surv(time, status) ~ rg, data=ssv)   
  res.cox <- coxph(Surv(time, status) ~ rg,data=ssv)
  p<-ggsurvplot(fit,data=ssv,
             pval = T, conf.int = T, 
             palette=c("#2980b9","#c0392b"), 
             risk.table.col = "strata", # Change risk table color by groups
             linetype = "strata", # Change line type by groups
             surv.median.line = "hv", # Specify median survival
             ggtheme = theme_few(),
             legend.labs=c("Low-risk","High-risk"),
             xlab = "Time (months)"
  )
 
  ggsave(paste0("iPMS_risk/",x,"_iPMS_risk_survival.pdf"),height = 5,width = 4.5,dpi=600)
  
  hr<-exp(res.cox$coefficients)
  c(hr,riskScore)
})

risk.iPCP<-lapply(c(x="Progress"), function(x){
  ssv<-icp.exp[,c(1:4,grep(x,colnames(icp.exp)),14:40)]
  colnames(ssv)[c(5:6)]<-c("time","status")           
  loc<-union(which(is.na(ssv[,5])),
             which(is.na(ssv[,6])))
  if(length(loc)>0){ssv<-ssv[-loc,]}
  ssv$status=as.numeric(gsub(":.*","",ssv$status))
  ssv$time=as.numeric(as.character(ssv$time))

  res.cox <- coxph(Surv(time, status) ~ AURKB+IL2RA+PD.1+RRM2+TIGIT+
                     CX3CL1+CXCL10+LILRA5+PLK1+SLAMF7+
                     TREM2+BTN3A3+LAG.3+LAIR.1,data=ssv)
  
  riskScore<-0
  for (i in names(res.cox$coefficients)) {
    riskScore<- riskScore+res.cox$coefficients[i]*ssv[,i]
  }
  group<-function(x,cutoff){
    group<-rep(NA,length(x))
    group[which(x>=cutoff)]<-1
    group[which(x<cutoff)]<-0
    group
  }  
  rg<-group(riskScore,cutoff = median(riskScore))
  ssv$rg<-rg
  res.cox <- coxph(Surv(time, status) ~ rg,data=ssv)
  fit <- survfit(Surv(time, status) ~ rg, data=ssv)
  p<-ggsurvplot(fit,data=ssv,
            pval = T, conf.int = F, 
            palette=c("#2980b9","#c0392b"), 
            risk.table.col = "strata", # Change risk table color by groups
            linetype = "strata", # Change line type by groups
            surv.median.line = "hv", # Specify median survival
            ggtheme = theme_few(),
            legend.labs=c("Low-risk","High-risk"),
            xlab = "Time (months)"
      )  
  ggsave(paste0("iPCP_risk/",x,"_iPCP_risk_survival.pdf"),height = 5,width = 4.5,dpi=600)
})








