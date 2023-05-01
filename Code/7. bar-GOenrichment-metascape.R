rm(list = ls())
gc()
setwd("Path to data")
###set file
go_list=read.csv("tox2-new/Enrichment_GO/_FINAL_GO.csv",h=T, sep = ",",comment="",check=F,stringsAsFactors = T)
go_list<-go_list[,c(7,9,11)]
###获取三个分类
BP<-go_list[which(go_list$Category=="GO Biological Processes"),]
CC<-go_list[which(go_list$Category=="GO Cellular Components"),]
MF<-go_list[which(go_list$Category=="GO Molecular Functions"),]
MF <- MF[1:10,]
BP <- BP[1:10,]
CC <- CC[1:10,]
go_list <-rbind(BP,CC) 
go_list<-rbind(go_list,MF)
###绘制GO富集柱状图
## numbers as data on x axis
go_list$number <- factor(rev(1:nrow(go_list)))
labels=as.character(go_list$Description) 
names(labels) = rev(1:nrow(go_list))
library(ggplot2)
go_list$padj <- ifelse(go_list$LogP < -1.301, '*', '')
p <- ggplot(data = go_list, aes(x = number, y = -LogP, fill = Category)) +
  scale_fill_manual(values=c("GO Cellular Components"= "#F8766D", "GO Biological Processes"= "#00BFC4", "GO Molecular Functions"= "#00BA38"))+
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  geom_text(aes(label = padj), vjust=0.75, hjust = -0.25, color = 1) +
  #scale_fill_manual(values = CPCOLS) + 
  theme_bw() + 
  scale_x_discrete(labels=labels) +
  #  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "Tox2")+
  xlab(" ")+
  scale_y_continuous(name = "-log10(P)",limits = c(0,16),breaks=seq(0,16,4))+
  theme(panel.grid.major=element_line(colour=NA),           
        panel.grid.minor=element_blank())+                  #移除网格线
  theme(panel.background = element_rect(colour = "black",size = 1.5))+    #边框加粗
  theme(plot.title = element_text(size = 17),
        text = element_text(size = 16),
        axis.title = element_text(size = 15),
        axis.text.x=element_text(size = 15))
p
ggsave(filename = "Metascape-tox2-go_enrich_of_diffgene.pdf", plot = p, height=6, width=10,dpi = 300)#5.5/9..2/8..
