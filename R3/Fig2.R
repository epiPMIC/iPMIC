library(ggplot2)
library(ggthemes)
library(scatterpie)
library(ggpubr)

####Fig2A
gd<-read.delim("OR.EWAS.txt",stringsAsFactors = F)
gd$class=factor(gd$class,levels=rev(c("cancer","lung","immune","other")))
ggdotchart(gd, x = "Trait", y = "lgOR",
           color = "class",                                # Color by groups
           palette = rev(c("#00AFBB", "#E7B800", "#FC4E07","black")), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           rotate = TRUE,                                # Rotate vertically
           group = "class",                                # Order by groups
           dot.size = 8,             # Adjust label parameters
           ggtheme = theme_pubclean()                        # ggplot2 theme
)
ggsave("OR.pdf",height=12,width = 8)

pie.gd<-reshape2::melt(gd[,c("Trait","Cold","Constitutive","Hot")])
pie.gd$value=pie.gd$value*100
pie.gd<-rbind(pie.gd,data.frame(Trait="legend",
                                variable=c("Cold","Constitutive","Hot"),
                                value=33.333))
ggplot(pie.gd,aes(x = '', y = value, fill = variable)) + 
  geom_bar(stat = 'identity', width = 1,color="white",size=1) + 
  theme_clean() + 
  coord_polar(theta = 'y', start = 0, direction = 1)+
  facet_wrap(.~Trait,ncol=3)+
  scale_fill_manual(values=c(Cold="#6AB5DD",
                             Hot="#EA614A",
                             Constitutive="#f6b93b"))
ggsave("OR.pie.pdf",height=12,width = 8)

####Fig2B
idmc<-read.delim(file = "iDMC_beta_0.3.txt",stringsAsFactors = F)
idmc$rho.q=p.adjust(idmc$X2,method = "BH")

ggplot(idmc,aes(y=X1,x=dis/1000,size=-log10(rho.q),color=Class))+
  geom_point(alpha=0.5)+
  geom_hline(yintercept = 0,linetype="dashed",color="gray20")+
  theme_few()+
  scale_size_area(max_size = 8)+
  ylab("Spearman's correlation coefficient")+
  xlab("Distant to the nearest PCG (kb)")+
  facet_wrap(Class~.)+
  ylim(-0.6,0.6)+
  scale_color_manual(values = c(Hot="#F7B731",Cold="#8854D0",Constitutive="#FA8231"))
ggsave("Rho.gene.pdf",height = 4,width = 9)

####Fig2C
gdd<-read.delim("FuncEnrich.txt",stringsAsFactors = F)
class<-read.delim("FuncEnrich_Class..txt",stringsAsFactors = F)
gdd$Description=factor(gdd$Description,levels = class$Description)
ggplot(gdd,aes(x=Description,y=NodeA,size=-LogP))+
  geom_point(shape=15,alpha=0.5)+
  geom_point(shape=0,color="black")+
  theme_pander()+
  scale_size_area(max_size = 10)+
  scale_color_gradient(low = "#ffffff",high = "#326280")+
  xlab(NULL)+ylab(NULL)+
  theme(axis.text.x = element_text(vjust = 0.5, 
                                   hjust = 1, angle = 45))
ggsave("function_enrich.pdf",height = 4,width = 8.5)
