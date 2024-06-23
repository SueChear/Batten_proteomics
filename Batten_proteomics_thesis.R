
#https://ab604.github.io/docs/bspr_workshop_2018/index.html#requirements####
#1. general volcano for W2####

library(magrittr)
library(ggplot2)


library(reshape2)
library(dplyr)
library(ggrepel)
setwd("P:/PhD project/60. Batten proteomics/analysis 2")
df<-read.csv("W2_forvolcanoplot.csv", header=T, sep=",")

View(df)

#replace NA with false
#df[is.na(df)]<-"FALSE"

View(df)
sig_pos <- df%>%dplyr::filter(log2FC>0)%>%dplyr::filter(SIGNIFICANT=="TRUE")
sig_neg <- df%>%dplyr::filter(log2FC<0)%>%dplyr::filter(SIGNIFICANT=="TRUE")
neutral<-df%>%dplyr::filter(SIGNIFICANT=="FALSE")


View(sig_pos)
#make a dataframe for just sig protein
df_sig<-rbind(sig_pos,sig_neg)
View(df_sig)

write.csv(df_sig,"sig_W2_RAW_fromR.csv", row.names=F)







df$genelabels <- ifelse( df$Gene == "TPP1"
                        |df$Gene == "PSAP"
                        | df$Gene == "NEFL" | df$Gene == "NAGA"| df$Gene == "DPP7"
                        | df$Gene == "NAGLU"
                        |df$Gene == "PPT1",
                        TRUE,FALSE)


p <-ggplot(df, aes(x=log2FC, y=-log10(pvalue)))+
  geom_point(data=sig_pos,aes(x=log2FC,y=-log10(pvalue)),color='red', size=1.5, alpha=0.4)+
  geom_point(data=sig_neg,aes(x=log2FC,y=-log10(pvalue)),color='spring green', size=1.5, alpha=0.4)+
  geom_point(data=neutral,aes(x=log2FC,y=-log10(pvalue)),color='grey',size=1.5, alpha=0.1)+
  
  labs(y="-log10 (p-value)")  +labs(x="Log2 Fold change (CLN3/CLN3-Cor)")+
  theme_classic()+
theme(plot.title = element_text(hjust = 0.5), 
      legend.position="right", 
      legend.title = element_blank())+
  geom_vline(xintercept=c(-0.014,0.014),lty=4,col="black",lwd=0.8)+
  geom_hline(yintercept = 0.65,lty=4,col="black",lwd=0.8)+
  theme(axis.text.x=element_text( face="bold",size=20,color='black'),
        axis.text.y=element_text(size=20, face="bold",color='black',vjust=-0.5),
        axis.title.y=element_text(size=30, face="bold"),
        axis.title.x=element_text(size=30, face="bold"),
        axis.line = element_line(colour = 'black', size = 2),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length = unit(.5, "cm"))+
  geom_text_repel(aes(x=log2FC, y=-log10(pvalue)),max.overlaps = Inf,label = ifelse(df$genelabels == TRUE, df$Gene,""), 
                  box.padding = unit(0.1, "lines"),hjust=1, color="black",fontface = 'bold', size=8)+
coord_cartesian(xlim=c(-0.5,0.5), clip="off")+scale_x_continuous(breaks=seq(-0.5,0.5,by=0.1))
p



#save as eps and png
ggsave("NW2_volcano.eps", device=cairo_ps, width = 35, height = 23, units = "cm",fallback_resolution = 600)
ggsave("NW2_volcano.png",width = 15, height = 10, dpi = 300, units = "in", device='png')

#2. general volcano for W6####

library(magrittr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)
setwd("P:/PhD project/60. Batten proteomics/analysis 2")
df<-read.csv("W6_forvolcanoplot.csv", header=T, sep=",")

View(df)

#replace NA with false
#df[is.na(df)]<-"FALSE"

View(df)
sig_pos <- df%>%dplyr::filter(log2FC>0)%>%dplyr::filter(SIGNIFICANT=="TRUE")
sig_neg <- df%>%dplyr::filter(log2FC<0)%>%dplyr::filter(SIGNIFICANT=="TRUE")
neutral<-df%>%dplyr::filter(SIGNIFICANT=="FALSE")

#make a dataframe for just sig protein
df_sig<-rbind(sig_pos,sig_neg)
View(df_sig)

write.csv(df_sig,"sig_W6_RAW_fromR.csv", row.names=F)


df$genelabels <- ifelse( df$Gene == "TPP1"
                             |df$Gene == "PSAP"
                       | df$Gene == "NEFL" | df$Gene == "NAGA"| df$Gene == "DPP7"
                        | df$Gene == "NAGLU"
                        |df$Gene == "PPT1",
                        TRUE,FALSE)


p <-ggplot(df, aes(x=log2FC, y=-log10(pvalue)))+
  geom_point(data=sig_pos,aes(x=log2FC,y=-log10(pvalue)),color='red', size=1.5, alpha=0.4)+
  geom_point(data=sig_neg,aes(x=log2FC,y=-log10(pvalue)),color='spring green', size=1.5, alpha=0.4)+
  geom_point(data=neutral,aes(x=log2FC,y=-log10(pvalue)),color='grey',size=1.5, alpha=0.1)+
  
  labs(y="-log10 (p-value)")  +labs(x="Log2 Fold change (CLN3/CLN3-Cor)")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())+scale_x_continuous(breaks=seq(-0.5,0.5,by=0.1))+
  geom_vline(xintercept=c(-0.014,0.014),lty=4,col="black",lwd=0.8)+
  geom_hline(yintercept = 0.65,lty=4,col="black",lwd=0.8)+
  theme(axis.text.x=element_text( face="bold",size=20,color='black'),
        axis.text.y=element_text(size=20, face="bold",color='black',vjust=-0.5),
        axis.title.y=element_text(size=30, face="bold"),
        axis.title.x=element_text(size=30, face="bold"),
        axis.line = element_line(colour = 'black', size = 2),
        axis.ticks = element_line(colour = "black", size = 2),
        axis.ticks.length = unit(.5, "cm"))+
  geom_text_repel(aes(x=log2FC, y=-log10(pvalue)),max.overlaps = Inf,label = ifelse(df$genelabels == TRUE, df$Gene,""), 
                  box.padding = unit(0.1, "lines"),hjust=1, color="black",fontface = 'bold',size=8) + 
  coord_cartesian(xlim=c(-0.3,0.5), clip="off")

p

#save as eps and png
ggsave("NW6_volcano.eps", device=cairo_ps, width = 35, height = 23, units = "cm",fallback_resolution = 600)
ggsave("NW6_volcano.png",width = 15, height = 10, dpi = 300, units = "in", device='png')


#3.k means clustering and pca####
library(factoextra)
library(tidyverse)
library(cluster)

setwd("P:/PhD project/60. Batten proteomics/analysis 2")
df<-read.csv("cleandatatransposed.csv", header=T, sep=",")
View(df)
df1<-na.omit(df)
df2<-df1[,-1]
df3<-scale(df2)

#find optimal kmean: try all three methods find the most common k mean
#method 1; elbow..k=4 is optimal
set.seed(123)
fviz_nbclust(df3, kmeans, method = "wss")

#method 2: average silhouette..k=2 is optimal
set.seed(123)
fviz_nbclust(df3, kmeans, method = "silhouette")

#method 3:gap statistic...k=4 is optimal
set.seed(123)
gap_stat <- clusGap(df3, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
print(gap_stat, method = "firstmax")
fviz_gap_stat(gap_stat)

#optimal k mean is 4
#visualize clustering with kmean=4
k2 <- kmeans(df3, centers = 4, nstart = 500)
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
fviz_cluster(k2, data = df2)

#to visualize distance matrix, red large dissimilarities vs teal fairly similar
distance <- get_dist(df2) #Euclidean
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

#rename df first column
x<-c("CLN3_W2","CLN3_W2","CLN3_W2","Control_W2","Control_W2","Control_W2","CLN3_W6","CLN3_W6","CLN3_W6","Control_W6","Control_W6","Control_W6")
df[["Gene"]]<-x
colnames(df)[colnames(df) == "Gene"] <- "Group"
View(df)

#do pca
set.seed(123)
final <- kmeans(df3, 4, nstart = 25)
res.pca<-prcomp(df3)
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
ind.coord$cluster <- factor(final$cluster)
ind.coord$Group <- df$Group

eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent


library(ggpubr)
library(ggplot2)


a<-ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "Group", size = 4,  legend = "right", ggtheme = theme_bw(),repel = T,
  font.label=c(20,"bold"),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) +
  theme(axis.text.x=element_text( face="bold",size=15,color='black'),
        axis.text.y=element_text(size=14, face="bold",color='black',vjust=-0.5),
        axis.title.y=element_text(size=30, face="bold"),
        axis.title.x=element_text(size=30, face="bold"))+
  stat_mean(aes(color = cluster), size = 2)+theme(plot.margin = unit(c(1,1,1,1), "lines"))
a

#save as eps and png
ggsave("pca.eps", device=cairo_ps, width = 20, height = 20, units = "cm",fallback_resolution = 600)
ggsave("pca.png",width = 20, height = 20, dpi = 300, units = "in", device='png')




#4-heatmap for significant protein only####
setwd("P:/PhD project/60. Batten proteomics/analysis 2")
df<-read.csv("anova_sig.csv", header=T, sep=",")
df$Gene<-gsub("X_","",as.character(df$Gene))
View(df)
df[, -c(1)] <- scale(df[, -c(1)])

rownames(df) = make.names(df$Gene, unique=TRUE)
df<-df[,-1]
View(df)
df_matrix <- data.matrix(df)%>%na.omit
View(df_matrix)

max(df_matrix)
min(df_matrix)
library(ggplot2)
library(ComplexHeatmap)
#https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
library(dendextend)
library(circlize)




b<-Heatmap(df_matrix,
           show_row_names = F,clustering_distance_rows = "spearman",
           row_dend_side="left",
           row_dend_width=unit(1.5, "cm"),
           column_names_gp = gpar(fontsize = 14),
           heatmap_legend_param = list(title = "Z score"),
           column_names_rot = 45,
           column_title = "",
           column_title_gp = gpar(fontsize = 20, fontface = "bold"),cluster_row_slices = FALSE,use_raster = FALSE,
           col=colorRamp2(c(-5,0,5),c("green", "black", "magenta")))
b

ggsave("heatmap_sig.eps", device=cairo_ps, width = 35, height = 50, units = "cm",fallback_resolution = 600)
ggsave("heatmap_sig.png",width = 15, height = 50, dpi = 300, units = "in", device='png')

 #save as pdf width 12 height 20

ggsave("heatmap_W2_sigprotein.eps", device=cairo_ps, width = 35, height = 26, units = "cm",fallback_resolution = 600)
ggsave("heatmap_W2_sigprotein.png",width = 15, height = 10, dpi = 300, units = "in", device='png')


#5. kegg_Week 2####
setwd("P:/PhD project/60. Batten proteomics/analysis 2")
df<-read.csv("kegg_w2.csv", header=T, sep=",")
View(df)

library(magrittr)
df2<-df%>%dplyr::select(c("Term","logp","Group"))
View(df2)

library(ggplot2)
library(reshape2)
library(dplyr)
library(viridis)

df2$Term <- factor(df2$Term, levels = df2$Term)
p<-ggplot(df2, aes(x=Term, y=logp,fill=Group)) + geom_bar(stat = "identity", width=0.6)+
  scale_fill_manual(values=c("#000080", "#FF0000"))+
  scale_y_continuous(breaks=seq(-20,60,by=10),labels=abs(seq(-20,60,by=10)))+
  labs(y=bquote(bold("-"~log[10]~"(P value)")))+theme_classic()+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(color="black",face="bold", size =14))+theme(legend.position="none")+
  theme(axis.title.x=element_text( face="bold",size=20,color='black'),
        axis.text.x=element_text( face="bold",size=16,color='black'))+
  coord_flip()

p

ggsave("kegg_W2.eps", device=cairo_ps, width = 35, height = 26, units = "cm",fallback_resolution = 600)
ggsave("kegg_W2.png",width = 15, height = 10, dpi = 300, units = "in", device='png')

#6. kegg_Week 6####
setwd("P:/PhD project/60. Batten proteomics/analysis 2")
df<-read.csv("kegg_w6.csv", header=T, sep=",")
View(df)

library(magrittr)
df2<-df%>%dplyr::select(c("Term","logp","Group"))
View(df2)

library(ggplot2)
library(reshape2)
library(dplyr)
library(viridis)

df2$Term <- factor(df2$Term, levels = df2$Term)
p<-ggplot(df2, aes(x=Term, y=logp,fill=Group)) + geom_bar(stat = "identity", width=0.6)+
  scale_fill_manual(values=c("#000080", "#FF0000"))+
  scale_y_continuous(breaks=seq(-20,60,by=10),labels=abs(seq(-20,60,by=10)))+
  labs(y=bquote(bold("-"~log[10]~"(P value)")))+theme_classic()+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(color="black",face="bold", size =14))+theme(legend.position="none")+
  theme(axis.title.x=element_text( face="bold",size=20,color='black'),
        axis.text.x=element_text( face="bold",size=16,color='black'))+
  coord_flip()

p

ggsave("kegg_W6.eps", device=cairo_ps, width = 35, height = 26, units = "cm",fallback_resolution = 600)
ggsave("kegg_W6.png",width = 15, height = 10, dpi = 300, units = "in", device='png')





#7. bar graph for lysosome W2####
setwd("P:/PhD project/60. Batten proteomics/analysis 2")
df<-read.csv("sig_w2.csv", header=T, sep=",")
View(df)

df1<-df%>%dplyr::select("Gene","m2")
View(df1)

#select lysosome genes
df2<-df1%>%filter(Gene%in% c("SCARB2", "HEXB", "HEXA", "CLTC", "LIPA", "SGSH", "NAGLU", "GM2A", "LAMP1", "CTSL", "SMPD1", "LAMP2", 
                             "ACP2", "CTSD", "CTSC", "CTSB", "ARSA", "MANBA", "ATP6AP1", "GAA", "AP3B1", "PLA2G15", "NPC1", "NPC2", 
                             "DNASE2", "ATP6V0D1", "ASAH1", "IDUA", "GBA", "GNS", "CLN5", "GGA2", "NEU1", "PSAP", "ATP6V0A1", "CTSA", 
                             "SORT1", "M6PR", "AP3D1", "NAGA", "IGF2R", "GALC", "GALNS", "GLB1", "MAN2B1", "TPP1", "PPT2", "GLA", "LGMN"))
View(df2)
write.csv(df2,"sig_W2_lysosome_fromR.csv", row.names=F)


#lysosome W6

df_L<-read.csv("sig_w6.csv", header=T, sep=",")
View(df_L)

dfL1<-df_L%>%dplyr::select("Gene","M6")
View(dfL1)

#select lysosome genes
dfL2<-dfL1%>%filter(Gene%in% c("SCARB2", "IDUA", "HEXB", "HEXA", "GBA", "LIPA", "GNS", "SGSH", "CLN5", "GM2A", "LAMP1", "NAGLU", 
                               "CTSL", "SMPD1", "LAMP2", "PSAP", "ACP2", "CTSD", "ARSB", "CTSC", "CTSB", "CTSA", "ARSA", "MANBA", 
                               "SORT1", "GAA", "M6PR", "NAGA", "IGF2R", "GALNS", "PLA2G15", "NPC2", "GLB1", "PPT1", "MAN2B1", "TPP1", 
                               "DNASE2", "PPT2", "GLA", "LGMN"))
View(dfL2)


write.csv(dfL2,"sig_W6_lysosome_fromR.csv", row.names=F)


dt<-read.csv("sig_W2W6_lysosome.csv", header=T, sep=",")

group_order <- factor(dt$Gene, level = c('ATP6V0A1', 'NEU1', 'NPC1', 'GGA2', 'ATP6AP1', 'AP3B1', 'ATP6V0D1', 'GALC', 'AP3D1', 'ASAH1',
                                         'CLTC', 'PPT1',
                                         'ACP2', 'PSAP', 'CTSD', 'SCARB2', 'GBA', 'M6PR', 'PLA2G15', 'SMPD1', 'ARSB', 'GAA', 'PPT2', 
                                         'HEXA', 'CTSB', 'NAGA',
                                         'CTSL', 'GM2A', 'CLN5', 'SORT1', 'IGF2R', 'LAMP2', 'LAMP1', 'GNS', 'HEXB', 'MANBA', 'ARSA', 
                                         'NPC2', 'TPP1', 'SGSH',
                                         'LIPA', 'GLA', 'GLB1', 'LGMN', 'GALNS', 'NAGLU', 'CTSA', 'MAN2B1', 'DNASE2', 'CTSC', 'IDUA'
                                         
))


View(dt)


ggplot(dt, aes(x=group_order, y=Values, color=Time))+geom_point(cex=5,alpha=0.5)+theme_classic()+
  
  ylab("")+xlab("")+
  theme(axis.text.x=element_text( face="bold",size=26,color='black'),
        axis.text.y=element_text(size=22, face="bold",color='black'),
        axis.title.y=element_text(size=50, face="bold"),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(.5, "cm"))+theme(legend.title = element_blank())+theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"))+coord_flip(ylim=c(1,1.25), clip="off")

ggsave("lysosome.jpg",width = 10, height = 16, dpi = 600, units = "in")

         
ggsave("lysosome.eps", device=cairo_ps, width = 35, height =25, units = "cm",fallback_resolution = 600)
ggsave("lysosome.png",width = 20, height = 15, dpi = 300, units = "in", device='png')


#8.bar graph for axon guidance W2####
setwd("P:/PhD project/60. Batten proteomics/analysis 2")
df<-read.csv("sig_w2.csv", header=T, sep=",")
View(df)

df1<-df%>%dplyr::select("Gene","m2")
View(df1)

#select axon genes
df2<-df1%>%filter(Gene%in% c("ROBO2", "GSK3B", "NRP1", "ROBO1", "CDC42", "PAK1", "PPP3CB", "EFNB3", "RGS3", "DPYSL5",
                             "DPYSL2", "CFL1", "PLXNA2", "RAC3", "FYN", "PLXNC1", "RAC1", "EPHB2", "SRGAP2", "HRAS", "PAK2", 
                             "SRGAP1", "EPHB1", "PLXNA3", "NTNG1", "ARHGEF12", "SEMA4D", "DCC", "L1CAM", "RHOA", "CDK5", "EPHA3"))
View(df2)
write.csv(df2,"sig_W2_axon_fromR.csv", row.names=F)


#axon W6

df_L<-read.csv("sig_w6.csv", header=T, sep=",")
View(df_L)

dfL1<-df_L%>%dplyr::select("Gene","M6")
View(dfL1)

#select axon genes
dfL2<-dfL1%>%filter(Gene%in% c("ROBO2", "SEMA3C", "GNAI1", "ROBO1", "CDC42", "PAK1", "ABLIM1", "PPP3CB", "NRAS", "EFNB3", "RGS3", 
                               "DPYSL5", "DPYSL2", "ABL1", "RAC3", "FYN", "PLXNC1", "SRGAP3", "RAC1", "PAK3", "SRGAP2", "HRAS", 
                               "SRGAP1", "EPHB1", "NTNG1", "ARHGEF12", "SEMA4D", "DCC", "SEMA4C", "RHOA", "RASA1", "KRAS"))
View(dfL2)


write.csv(dfL2,"sig_W6_axon_fromR.csv", row.names=F)


dt<-read.csv("sig_W2W6_axon.csv", header=T, sep=",")

group_order <- factor(dt$Gene, level = c('CDK5','CFL1','EPHA3','EPHB2','HRAS', 'L1CAM', 'PAK1', 'PAK2', 'PLXNA2','PLXNA3','NRP1','DCC',  
                                         'EFNB3',
                                         'SEMA3C', 'RGS3','PLXNC1', 'SRGAP1','ROBO2', 'ROBO1', 'RAC3', 'GNAI1','SEMA4D','SEMA4C', 'ABL1',
                                         'DPYSL5','RAC1','SRGAP2','KRAS', 'CDC42', 'GSK3B','SRGAP3','PAK3','FYN','RASA1', 'NRAS','ARHGEF12',
                                         'DPYSL2','RHOA','PPP3CB','ABLIM1','EPHB1',
                                         'NTNG1'
                                         

))


View(dt)


ggplot(dt, aes(x=group_order, y=Values, color=Time))+geom_point(cex=5,alpha=0.5)+theme_classic()+
  
  ylab("")+xlab("")+
  theme(axis.text.x=element_text( face="bold",size=26,color='black'),
        axis.text.y=element_text(size=22, face="bold",color='black'),
        axis.title.y=element_text(size=50, face="bold"),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(.5, "cm"))+theme(legend.title = element_blank())+theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"))+coord_flip(ylim=c(0.85,1), clip="off")

ggsave("axon.jpg",width = 11, height = 18, dpi = 600, units = "in")




ggsave("axon.eps", device=cairo_ps, width = 35, height =25, units = "cm",fallback_resolution = 600)
ggsave("axon.png",width = 20, height = 15, dpi = 300, units = "in", device='png')


#9.bar graph for endocytosis W2####
setwd("P:/PhD project/60. Batten proteomics/analysis 2")
df<-read.csv("sig_w2.csv", header=T, sep=",")
View(df)

df1<-df%>%dplyr::select("Gene","m2")
View(df1)

#select endo genes
df2<-df1%>%filter(Gene%in% c("ARF1", "SRC", "VPS4A", "ARPC1A", "NEDD4L", "AGAP1", "ARRB1", "CBL", "CDC42", "PSD2", 
                             "KIF5C", "KIF5A", "MVB12B", "PIP5K1C", "HRAS", "RAB11FIP5", "IQSEC1", "IST1", "EPS15L1", 
                             "STAM", "RBSN", "RHOA", "EHD1", "EHD3", "RABEP1", "BIN1", "AMPH", "CHMP6", "SMAP1", "VPS25", 
                             "HSPA1B", "ARF5", "HSPA1A"))
View(df2)
write.csv(df2,"sig_W2_endo_fromR.csv", row.names=F)


#endocytosis W6

df_L<-read.csv("sig_w6.csv", header=T, sep=",")
View(df_L)

dfL1<-df_L%>%dplyr::select("Gene","M6")
View(dfL1)

#select endo genes
dfL2<-dfL1%>%filter(Gene%in% c("RAB5C", "SRC", "AGAP2", "NEDD4L", "AGAP1", "ARRB1", "AGAP3", "CDC42", "CYTH2", "PSD2", 
                               "KIF5C", "PSD3", "MVB12B", "PIP5K1C", "RAB11FIP2", "HRAS", "AP2M1", "GIT1", "SH3GL2", 
                               "SH3GL1", "SH3GLB1", "SH3GLB2", "IQSEC1", "VTA1", "IQSEC3", "RHOA", "ARFGAP1", "EPN2", "EHD1", 
                               "AMPH", "VPS45", "CHMP6", "SMAP1", "ARF6"))
View(dfL2)


write.csv(dfL2,"sig_W6_endo_fromR.csv", row.names=F)



#select upregulated endo genes
df_U<-read.csv("sig_w6.csv", header=T, sep=",")



df_U<-df_U%>%dplyr::select("Gene","M6")



dfU2<-df_U%>%filter(Gene%in% c("VPS29", "TFRC", "SH3KBP1", "ARPC1B", "VPS4B", "ASAP3", "EGFR", "SNX3", "EEA1", "SNX1", "SNX2", 
                              "CHMP1B", "SNX5", "RAB11FIP5", "SNX6", "ARFGEF1", "GIT2", "GBF1", "KIAA0196", "HLA-A", "IGF2R", 
                              "PML", "DNM2", "RUFY1", "RAB31", "TRAF6", "CHMP2A", "CHMP3", "VPS25", "SPG20"))
View(dfU2)


write.csv(dfU2,"sig_W6_UP_endo_fromR.csv", row.names=F)




dt<-read.csv("sig_W2W6_downonly_endo.csv", header=T, sep=",")

group_order <- factor(dt$Gene, level = c('ARF1','ARF5','ARPC1A', 'BIN1', 'CBL', 'EHD3', 'EPS15L1','IST1', 'KIF5A','RAB11FIP5', 'RABEP1',
                                         'RBSN', 'STAM','VPS25','VPS4A', 'AGAP2','SH3GL2','PSD3','RAB11FIP2','ARRB1','IQSEC3','AGAP3',
                                         'SRC', 'ARF6','SMAP1','CYTH2','AGAP1','AMPH', 'CDC42', 'CHMP6','HRAS','PSD2','IQSEC1','MVB12B',
                                         'EHD1','GIT1','RHOA','RAB5C','VTA1','VPS45','SH3GLB2','PIP5K1C','SH3GL1','KIF5C','NEDD4L',
                                         'EPN2','AP2M1','ARFGAP1','SH3GLB1'
                                         
))


View(dt)


ggplot(dt, aes(x=group_order, y=Values, color=Time))+geom_point(cex=5,alpha=0.5)+theme_classic()+
  
  ylab("")+xlab("")+
  theme(axis.text.x=element_text( face="bold",size=26,color='black'),
        axis.text.y=element_text(size=22, face="bold",color='black'),
        axis.title.y=element_text(size=50, face="bold"),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(.5, "cm"))+theme(legend.title = element_blank())+theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"))+coord_flip(ylim=c(0.85,1), clip="off")

ggsave("endo.jpg",width = 10, height = 19, dpi = 600, units = "in")


ggsave("endo.eps", device=cairo_ps, width = 35, height =25, units = "cm",fallback_resolution = 600)
ggsave("endo.png",width = 20, height = 15, dpi = 300, units = "in", device='png')



#10. bar graph for ribosome W2####
setwd("P:/PhD project/60. Batten proteomics/analysis 2")
df<-read.csv("sig_w2.csv", header=T, sep=",")
View(df)

df1<-df%>%dplyr::select("Gene","m2")
View(df1)

#select ribosome genes
df2<-df1%>%filter(Gene%in% c("RPL4", "RPL5", "RPL30", "RPL3", "RPL31", "RPL34", "RPLP0", "MRPS10", "RPL8", "RPL9", 
                             "RPL6", "MRPL4", "RPS4X", "RPS15", "RPS14", "RPL7A", "RPS17", "RPS16", "MRPL1", "RPS19", 
                             "RPL18A", "RPS18", "RPL36", "RPL35", "RPS11", "RPS10", "RPS13", "RPS9", "RPL21", "RPS8", 
                             "RPL23", "RPL22", "RPS6", "RPS3A", "MRPS6", "MRPS9", "RPL37A", "RPL24", "RPL27", "RPL26", 
                             "RPL29", "RPL28", "RPS4Y2", "RPL10", "RPL12", "RPL11", "MRPL17", "RPS27L", "MRPL13", "MRPL11", 
                             "RPS15A", "RPS3", "RPL14", "RPL13", "RPS2", "RPL17", "RPL35A", "MRPL28", "RPL23A", "MRPL23", 
                             "RPS26", "RPS25", "RPS28", "RPL27A", "RPS20", "RPS24", "RPS23"))
View(df2)
write.csv(df2,"sig_W2_ribosome_fromR.csv", row.names=F)


#ribosome W6

df_L<-read.csv("sig_w6.csv", header=T, sep=",")
View(df_L)

dfL1<-df_L%>%dplyr::select("Gene","M6")
View(dfL1)

#select ribosome genes
dfL2<-dfL1%>%filter(Gene%in% c("RPL4", "RPL5", "RPL30", "RPL3", "RPL32", "RPL31", "RPL34", "RPLP0", "RPL8", 
                               "RPL9", "RPL6", "RPL7", "RPS4X", "RPS15", "RPS14", "RPL7A", "RPS17", "RPS16", 
                               "RPS19", "RPL18A", "RPS18", "RPL36", "RPL35", "RPL38", "RPS11", "RPS13", "RPS9", 
                               "RPL21", "RPS7", "RPS8", "RPL23", "RPS5", "RPL22", "RPS6", "RPL13A", "RPS3A", 
                               "RPSA", "RPL37A", "RPL24", "RPL27", "RPL26", "RPL29", "RPL28", "RPS4Y2", "RPL10", 
                               "RPL11", "RPS27L", "RPS15A", "RPS3", "RPL14", "RPL13", "RPS2", "RPL15", "RPS27A", 
                               "RPL18", "RPL17", "RPL19", "RPL35A", "RPL23A", "MRPL23", "RPS26", "RPS25", "RPL27A", 
                               "FAU", "RPS24", "RPS23"))
View(dfL2)


write.csv(dfL2,"sig_W6_ribosome_fromR.csv", row.names=F)


dt<-read.csv("sig_W2W6_ribosome.csv", header=T, sep=",")

group_order <- factor(dt$Gene, level = c('RPL12','RPS10','RPS20','RPS28','MRPS10','MRPS6','MRPS9','MRPL28','MRPL13', 'MRPL4',
                                         'MRPL1','MRPL17','MRPL11','RPL30','RPS17','RPSA','RPS27A','RPS18','RPLP0','RPS14',
                                         'RPS13','RPS19','RPS15','RPS7','RPL22','RPL15','RPS3','RPL5','RPL18','RPS25','RPL23A',
                                         'RPL23','MRPL23','RPS2','RPL38','RPL9','RPL36','RPS3A','RPS27L','RPS15A','RPL11','RPL21',
                                         'RPS24', 'RPL32','RPL27A','RPS5','RPS4X','RPS26','RPL24','RPL35A','RPL17','FAU','RPL7',
                                         'RPL10','RPS9','RPL26','RPS11', 'RPL27','RPL6','RPL31','RPL3','RPL4','RPL35','RPL34',
                                         'RPL7A','RPL29','RPS16','RPL13A','RPS8','RPL13','RPS23','RPL19','RPS6','RPL14','RPL8',
                                         'RPL18A','RPL37A','RPL28','RPS4Y2'
                                         
))


View(dt)


ggplot(dt, aes(x=group_order, y=Values, color=Time))+geom_point(cex=5,alpha=0.5)+theme_classic()+
  
  ylab("")+xlab("")+
  theme(axis.text.x=element_text( face="bold",size=26,color='black'),
        axis.text.y=element_text(size=22, face="bold",color='black'),
        axis.title.y=element_text(size=18, face="bold"),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.ticks.length = unit(.5, "cm"))+theme(legend.title = element_blank())+theme(legend.position = "none")+
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"))+coord_flip(ylim=c(1,1.15), clip="off")

ggsave("ribosome.jpg",width = 10, height =22, dpi = 600, units = "in")
ggsave("ribosome.eps", device=cairo_ps, width = 35, height =25, units = "cm",fallback_resolution = 600)
ggsave("ribosome.png",width = 20, height = 15, dpi = 300, units = "in", device='png')
