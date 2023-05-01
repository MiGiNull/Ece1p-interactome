rm(list = ls())
gc()
setwd("Path to data")
###set file
data=read.csv("gene_with_interaction_strength.csv",h=T, sep = ",",comment="",check=F,stringsAsFactors = T,row.names = 1)
data[data!=0]<-1

###fisher检验（超几何）
##两个组合
m_two<- matrix(nr=8,nc=8)
for (i in 1:7){
  for (j in (i+1):8) {
    q<-length(which(data[,i]==1 & data[,j]==1))
    if (length(which(data[,j]==1))>length(which(data[,i]==1))){
      m<-length(which(data[,j]==1))
      k<-length(which(data[,i]==1))
    }else{
      m<-length(which(data[,i]==1))
      k<-length(which(data[,j]==1))
    }
    n<-nrow(data)-m
    pvalue<-phyper(q-1, m, n, k, lower.tail=F) 
    m_two[i,j]<-p.adjust(pvalue,method='fdr')
  }
}
m_two<-as.data.frame(m_two)
colnames(m_two)<-colnames(data)
rownames(m_two)<-colnames(data)
m_two[m_two>0.05]<-1

library(corrplot)
# #P-下三角
m_two[is.na(m_two)]<-0
p<-corrplot(as.matrix(-log10(m_two)),method = "color",is.corr = F,p.mat =as.matrix(m_two),insig = "label_sig",
            sig.level =0.05 ,rect.lwd=1,
            diag= F,
            col = colorRampPalette(c("white","steelblue"))(50),tl.col = "black",tl.cex=2,,type="upper",
            cl.cex = 1.7,cl.ratio=0.3)
dev.off()


#tox3&fl
p.adjust(phyper(29, 83, 522-83, 77, lower.tail=F),method='fdr')
 









##三个组合[都不显著]
a<-t(combn(8,3))
m_three<- matrix(nr=8,nc=nrow(a))
for (i in 1:nrow(a)) {
  col=a[i,]
  q<-length(which(data[,col[1]]==1 & data[,col[2]]==1 & data[,col[3]]==1))
  x<-length(which(data[,col[1]]==1 | data[,col[2]]==1 | data[,col[3]]==1))
  for (j in col) {
    m<-(x-length(which(data[,j]==1))+q)
    k<-length(which(data[,j]==1))
    n<-nrow(data)-m
    if(m<k){   #备择假设为less：在次品个数是多数的情况下，overlap越低，从超几何分布来看，发生的概率越小。
      pvalue<-phyper(q, m, n, k, lower.tail=T)
    }else{
      pvalue<-phyper(q-1, m, n, k, lower.tail=F) #备择假设为more：在次品个数是少数的情况下，overlap越高，从超几何分布来看，发生的概率越小。
    }
    m_three[j,i]<-pvalue
  }
}
m_three<-t(m_three)
m_three[is.na(m_three)]<-0
mean<-apply(m_three,1, function(x){sum(x)/3})
final_three<-as.matrix(m_three[which(apply(m_three,1,mean)<0.05),])
##四个组合[不显著]
a<-t(combn(8,4))
m_four<- matrix(nr=8,nc=nrow(a))
for (i in 1:nrow(a)) {
  col=a[i,]
  q<-length(which(data[,col[1]]==1 & data[,col[2]]==1 & data[,col[3]]==1 & data[,col[4]]==1))
  print(q)
  x<-length(which(data[,col[1]]==1 | data[,col[2]]==1 | data[,col[3]]==1 | data[,col[4]]==1))
  for (j in col) {
    m<-(x-length(which(data[,j]==1))+q)
    k<-length(which(data[,j]==1))
    n<-nrow(data)-m
    if(m<k){   #备择假设为less：在次品个数是多数的情况下，overlap越低，从超几何分布来看，发生的概率越小。
      pvalue<-phyper(q, m, n, k, lower.tail=T)
    }else{
      pvalue<-phyper(q-1, m, n, k, lower.tail=F) #备择假设为more：在次品个数是少数的情况下，overlap越高，从超几何分布来看，发生的概率越小。
    }
    m_four[j,i]<-pvalue
  }
}
m_four<-t(m_four)
m_four[is.na(m_four)]<-0
mean<-apply(m_four,1,function(x){sum(x)/4})
final_four<-as.matrix(m_four[which(apply(m_four,1,mean)<0.05),])
##五个组合[不显著]
a<-t(combn(8,5))
m_five<- matrix(nr=8,nc=nrow(a))
for (i in 1:nrow(a)) {
  col=a[i,]
  q<-length(which(data[,col[1]]==1 & data[,col[2]]==1 & data[,col[3]]==1 & data[,col[4]]==1 & data[,col[5]]==1))
  x<-length(which(data[,col[1]]==1 | data[,col[2]]==1 | data[,col[3]]==1 | data[,col[4]]==1 | data[,col[5]]==1))
  for (j in col) {
    m<-(x-length(which(data[,j]==1))+q)
    k<-length(which(data[,j]==1))
    n<-nrow(data)-m
    if(m<k){   #备择假设为less：在次品个数是多数的情况下，overlap越低，从超几何分布来看，发生的概率越小。
      pvalue<-phyper(q, m, n, k, lower.tail=T)
    }else{
      pvalue<-phyper(q-1, m, n, k, lower.tail=F) #备择假设为more：在次品个数是少数的情况下，overlap越高，从超几何分布来看，发生的概率越小。
    }
    m_five[j,i]<-pvalue
  }
}
m_five<-t(m_five)
m_five[is.na(m_five)]<-0
mean<-apply(m_five,1,function(x){sum(x)/5})
final_five<-as.matrix(m_five[which(apply(m_five,1,mean)<0.05),])
##六个组合[不显著]
a<-t(combn(8,6))
m_six<- matrix(nr=8,nc=nrow(a))
for (i in 1:nrow(a)) {
  col=a[i,]
  q<-length(which(data[,col[1]]==1 & data[,col[2]]==1 & data[,col[3]]==1 & data[,col[4]]==1 & data[,col[5]]==1 & data[,col[6]]==1))
  x<-length(which(data[,col[1]]==1 | data[,col[2]]==1 | data[,col[3]]==1 | data[,col[4]]==1 | data[,col[5]]==1 | data[,col[6]]==1))
  for (j in col) {
    m<-(x-length(which(data[,j]==1))+q)
    k<-length(which(data[,j]==1))
    n<-nrow(data)-m
    if(m<k){   #备择假设为less：在次品个数是多数的情况下，overlap越低，从超几何分布来看，发生的概率越小。
      pvalue<-phyper(q, m, n, k, lower.tail=T)
    }else{
      pvalue<-phyper(q-1, m, n, k, lower.tail=F) #备择假设为more：在次品个数是少数的情况下，overlap越高，从超几何分布来看，发生的概率越小。
    }
    m_six[j,i]<-pvalue
  }
}
m_six<-t(m_six)
m_six[is.na(m_six)]<-0
mean<-apply(m_six,1,function(x){sum(x)/6})
final_six<-as.matrix(m_six[which(apply(m_six,1,mean)<0.05),])
##七个组合[不显著]
a<-t(combn(8,7))
m_seven<- matrix(nr=8,nc=nrow(a))
for (i in 1:nrow(a)) {
  col=a[i,]
  q<-length(which(data[,col[1]]==1 & data[,col[2]]==1 & data[,col[3]]==1 & data[,col[4]]==1 & data[,col[5]]==1 & data[,col[6]]==1 & data[,col[7]]==1))
  x<-length(which(data[,col[1]]==1 | data[,col[2]]==1 | data[,col[3]]==1 | data[,col[4]]==1 | data[,col[5]]==1 | data[,col[6]]==1 | data[,col[7]]==1))
  for (j in col) {
    m<-(x-length(which(data[,j]==1))+q)
    k<-length(which(data[,j]==1))
    n<-nrow(data)-m
    if(m<k){   #备择假设为less：在次品个数是多数的情况下，overlap越低，从超几何分布来看，发生的概率越小。
      pvalue<-phyper(q, m, n, k, lower.tail=T)
    }else{
      pvalue<-phyper(q-1, m, n, k, lower.tail=F) #备择假设为more：在次品个数是少数的情况下，overlap越高，从超几何分布来看，发生的概率越小。
    }
    m_seven[j,i]<-pvalue
  }
}
m_seven<-t(m_seven)
m_seven[is.na(m_seven)]<-0
mean<-apply(m_seven,1,function(x){sum(x)/7})
final_seven<-as.matrix(m_six[which(apply(m_seven,1,mean)<0.05),])
##八个组合[不显著]
a<-t(combn(8,8))
m_eight<- matrix(nr=8,nc=nrow(a))
for (i in 1:nrow(a)) {
  col=a[i,]
  q<-length(which(data[,col[1]]==1 & data[,col[2]]==1 & data[,col[3]]==1 & data[,col[4]]==1 & data[,col[5]]==1 & data[,col[6]]==1 & data[,col[7]]==1 & data[,col[8]]==1))
  x<-length(which(data[,col[1]]==1 | data[,col[2]]==1 | data[,col[3]]==1 | data[,col[4]]==1 | data[,col[5]]==1 | data[,col[6]]==1 | data[,col[7]]==1 | data[,col[8]]==1))
  for (j in col) {
    m<-(x-length(which(data[,j]==1))+q)
    k<-length(which(data[,j]==1))
    n<-nrow(data)-m
    l<-phyper(q-1, m, n, k, lower.tail=F) 
    m_eight[j,i]<-l
  }
}
