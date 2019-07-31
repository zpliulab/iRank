rm(list=ls())
library(R.matlab)#D:/shx_bioinformatics/TCGA_Project/Copy_of_reg_s
#setwd("D:/shx_bioinformatics/TCGA_Project/Copy_of_reg_s")
#setwd("C:/Users/shx/Desktop/used_for_plot_ROC_and_static")
#setwd("D:/shx_bioinformatics/TCGA_Project/used_for_plot_ROC_and_static")

#data <- readMat("reg_s_33_new.mat") #print(data)
setwd("D:/shx_bioinformatics/TCGA_Project/RWH_result")
setwd("D:/git_data/procedure/HCC_procedure/Other_comparing_algorithm/RWH")
data <- readMat("RWH_new.mat") #print(data)
Auc <-as.matrix(unlist(data[1]))#D <-as.matrix(unlist(data[2]))
pdf(paste("RWH_","re0.pdf",sep = ""))
boxplot(Auc,col=c("red"))
dev.off() 
#N_id <-data.frame(as.matrix(unlist(data[3]), nrow=33,ncol=1000))
D <-as.data.frame(data[2])# N_id <-as.data.frame(data[3])# N_pr <-as.data.frame(data[4])
N_ra <-as.data.frame(data[5])

ROC <-function(NORMA,DISEA){
  library(e1071)
  library(pROC) 
  label<-matrix(c(rep(0,nrow(NORMA)),rep(1,nrow(DISEA))))
  mdata<-rbind(NORMA,DISEA)
  svm <-NULL
  for(i in 1:nrow(label))
  {
    mtrain  <-label[-i,]
    predata <-mdata[i,]
    traindata <-mdata[-i,]
    model <-svm(traindata,mtrain,type="C-classification",kernel="radial")
    pred  <-predict(model,t(predata),decision.values=TRUE)
    f     <-attr(pred,"decision.values")
    prb   <-1/(1+exp(-f))
    svm   <-rbind(svm,prb)
  }
  re <-cbind(label,svm)
  return (re)
} 
Se_Sp <-function(Input1)
{
  library(caret)  
  threshold=0.5
  predicted_values<-ifelse(as.numeric(Input1[,2])<threshold,1,0)#max1[,2]
  actual_values<-Input1[,1]#max1[,1]
  conf_matrix<-table(predicted_values,actual_values)
  conf_matrix
  Sse <-sensitivity(conf_matrix) # Sse
  Ssp <-specificity(conf_matrix)# Ssp
  Spre <-precision(conf_matrix)
  F1_Score <-2*(Sse*Ssp)/(Ssp+Sse)#F1_Score
  rer <-confusionMatrix(conf_matrix)#mode(rer)
  Acc <-rer[["overall"]][["Accuracy"]]
  re <-cbind(cbind(Sse,Ssp),cbind(F1_Score,Acc))
  return(re)
}

sp <-cbind()
se <-cbind()
se_sp <-rbind()
Auc_r <-rbind()
#dim(Auc)
for(i in 1:length(Auc)){
#d1 <-N_ra[,i]
 disea<-as.matrix(D[,3])#疾病是1 疾病的rank
norma1<-as.matrix(N_ra[,i])#正常是0 disease 的rank
max1 <-ROC(norma1,disea)
y1 <-roc(as.factor(max1[,1]),as.vector(max1[,2]),plot=F,print.auc=F,print.thres=F,col="blue",lty=1,lwd=2)
Auc_r <-rbind(Auc_r,y1[["auc"]])
data <-data.frame(x=y1$specificities, y=y1$sensitivities)#1-
se <-cbind(se,data[,1])#
sp <-cbind(sp,data[,2])
#se <-cbind(se,as.matrix(unlist(y1[["sensitivities"]])))
#sp <-cbind(sp,as.matrix(unlist(y1[["specificities"]])))
re <-Se_Sp(max1)
se_sp <-rbind(se_sp,re)
}
length(Auc_r)
max(Auc_r)
mean(Auc_r)
sd(Auc_r)
min(Auc_r)
pdf(paste("RWH_","re01.pdf",sep = ""))
boxplot(Auc_r,col=c("red"))
dev.off() 
me <-cbind(mean(se_sp[,1]),mean(se_sp[,2]),mean(se_sp[,3]),mean(se_sp[,4]),mean(Auc_r))
#me
mode(me)
st <-cbind(sd(se_sp[,1]),sd(se_sp[,2]),sd(se_sp[,3]),sd(se_sp[,4]),sd(Auc_r))
#st
library(gridExtra)
library(grid)
pdf(paste("RWH_","re11.pdf",sep = ""))
da <-rbind(cbind(cbind("mean",round(me,5))),cbind("sd",cbind(round(st,5))))
#da
colnames(da) <- c("Index","Se","Sp","F1","Acc","Auc")
tt3 <- ttheme_minimal(
  core=list(bg_params = list(fill = blues9[1:4], col=NA),
            fg_params=list(fontface=3)),
  colhead=list(fg_params=list(col="navyblue", fontface=4L)),
  rowhead=list(fg_params=list(col="orange", fontface=3L)))
grid.table(da,theme=tt3)
dev.off() 
me <-cbind(mean(se_sp[,1]),mean(se_sp[,2]),mean(se_sp[,3]),mean(se_sp[,4]))
#me
mode(me)
st <-cbind(sd(se_sp[,1]),sd(se_sp[,2]),sd(se_sp[,3]),sd(se_sp[,4]))

#
pdf(paste("RWH_","re1.pdf",sep = ""))
da <-rbind(cbind(cbind("mean",round(me,5)),round(mean(Auc),5)),cbind("sd",cbind(round(st,5),round(sd(Auc),5))))
#da
colnames(da) <- c("Index","Se","Sp","F1","Acc","Auc")
tt3 <- ttheme_minimal(
  core=list(bg_params = list(fill = blues9[1:4], col=NA),
            fg_params=list(fontface=3)),
  colhead=list(fg_params=list(col="navyblue", fontface=4L)),
  rowhead=list(fg_params=list(col="orange", fontface=3L)))
grid.table(da,theme=tt3)
dev.off() 
##  最大值 AUC_r 最小AUC_r Median AUC_r
pdf(paste("RWH_R_","re22.pdf",sep = ""))
disea<-as.matrix(D[,3])#疾病是1 疾病的rank
norma11<-as.matrix(N_ra[,which(Auc_r==max(Auc_r))])#正常是0 disease 的rank
max11 <-ROC(norma11,disea)
norma22<-as.matrix(N_ra[,which(Auc_r==min(Auc_r))])
min11 <-ROC(norma22,disea)
kkk <-which(Auc_r==median(Auc_r))
kkk
#kkk <-which(Auc==mean(Auc))
if(length(kkk)>0){
  norma33 <-as.matrix(N_ra[,which(Auc_r==median(Auc_r))[1]])#1
}else{
  Auuc <-order(Auc_r)[(length(Auc_r)/2)+1]#排序之后的中间的位置
  norma33 <-as.matrix(N_ra[,Auuc])  
}
med <-ROC(norma33,disea)
y11 <-roc(as.factor(max11[,1]),as.vector(max11[,2]),plot=T,print.auc=F,print.thres=F,col="blue",lty=1,lwd=2)
y22 <-roc(as.factor(min11[,1]),as.vector(min11[,2]),plot=T,print.auc=F,print.thres=F,col="green",lty=1,lwd=2,add=T)
y33 <-roc(as.factor(med[,1]),as.vector(med[,2]),plot=T,print.auc=T,print.thres=F,col="red",lty=1,lwd=2,add=T)
dev.off() 
##  最大值 AUC 最小AUC Median AUC
pdf(paste("Reg_S","re2.pdf",sep = ""))
disea<-as.matrix(D[,3])#疾病是1 疾病的rank
norma11<-as.matrix(N_ra[,which(Auc==max(Auc))])#正常是0 disease 的rank
max11 <-ROC(norma11,disea)
norma22<-as.matrix(N_ra[,which(Auc==min(Auc))])
min11 <-ROC(norma22,disea)
kkk <-which(Auc==median(Auc))
kkk
#kkk <-which(Auc==mean(Auc))
if(length(kkk)>0){
norma33 <-as.matrix(N_ra[,which(Auc==median(Auc))[2]])#1
}else{
  Auuc <-order(Auc)[(length(Auc)/2)+1]#排序之后的中间的位置
  norma33 <-as.matrix(N_ra[,Auuc])  
}
med <-ROC(norma33,disea)
y11 <-roc(as.factor(max11[,1]),as.vector(max11[,2]),plot=T,print.auc=F,print.thres=F,col="blue",lty=1,lwd=2)
y22 <-roc(as.factor(min11[,1]),as.vector(min11[,2]),plot=T,print.auc=F,print.thres=F,col="green",lty=1,lwd=2,add=T)
y33 <-roc(as.factor(med[,1]),as.vector(med[,2]),plot=T,print.auc=T,print.thres=F,col="red",lty=1,lwd=2,add=T)
dev.off() 

library(ggplot2)
dim(se)# y轴
dim(sp)# x轴
dim(se_sp)
x1 <-rbind()
yy1 <-rbind()
yy2 <-rbind()
yy3 <-rbind()
seq1<-seq(from=1,to=0,by=-0.015)
seq1
for (x in seq1)
{ 
  x1 <-rbind(x1,x)
  mini <-rbind()
  maxx <-rbind()
  mea  <-rbind()
  for(jj in 1:dim(sp)[2])
  {
    l1 <-which(abs(sp[,jj]-x)<0.04)
    l1
    if(length(l1)>0)
    {
      yy <-se[l1,jj]
      mini<-rbind(mini,min(yy))
      maxx<-rbind(maxx,max(yy))
    }else{
      mini <-rbind(mini,0)
      maxx<-rbind(maxx,0)
    }
  }
  l2 <-cbind(min(mini),max(maxx))
  yy1 <-rbind(yy1,min(mini))
  yy2 <-rbind(yy2,max(maxx))
  yy3 <-rbind(yy3,(min(mini)+max(maxx))/2)
}
library(ggplot2)
x <-(1-x1)# 1-sp
dat <-as.data.frame(cbind(x,cbind(yy1,yy2,yy3)))
names(dat) = c("x","Model 1","Model 2","Model 3")
pdf(paste("RWH_","f1.pdf",sep = ""))
ggplot()+
geom_ribbon(data=dat,aes(x,ymin=`Model 1`,ymax=`Model 2`), fill="pink")+
geom_line(aes(x=x,y=yy3,color=""))+
scale_color_manual(values=c("purple"))+
scale_size_manual(values=c(1.5))+
labs(x="1 - Specificity", y="Sensitivity", colour="")+
theme_classic()
dev.off() 



#library(ROCR)
#pred <- prediction(as.vector(max1[,2]),as.factor(max1[,1]))
#perf <- performance(pred,"spec","sens")
#mode(perf)
#plot(perf)#这个是反的

library(MLmetrics)
if(F){
  library(caret)  
  threshold=0.5
  predicted_values<-ifelse(as.numeric(max1[,2])<threshold,1,0)
  actual_values<-max1[,1]
  conf_matrix<-table(predicted_values,actual_values)
  conf_matrix
  Sse <-sensitivity(conf_matrix)
  Sse
  Ssp <-specificity(conf_matrix)
  Ssp
  Spre <-precision(conf_matrix)
  F1_Score <-2*(Sse*Ssp)/(Ssp+Sse)
  F1_Score
  rer <-confusionMatrix(conf_matrix)
  mode(rer)
  Acc <-rer[["overall"]][["Accuracy"]]
 

#pdf("f1.pdf")
#par(mfcol=c(5,3))
#install.packages("ggplot2")
library(ggplot2)
sp <-rbind()
se <-rbind()
y1 <-roc(as.factor(max1[,1]),as.vector(max1[,2]),plot=F,print.auc=F,print.thres=F,col="blue",lty=1,lwd=2)
unlist(y1[["sensitivities"]])
unlist(y1[["specificities"]])
}
#y2=roc(as.factor(min1[,1]),as.vector(min1[,2]),plot=T,print.auc=F,print.thres=F,col="green",lty=1,lwd=2,add = TRUE)
#polygon(x,y,col=gray(0.8))
#roc(as.factor(med1[,1]),as.vector(med1[,2]),plot=T,print.auc=T,print.thres=F,col="red",lty=1,lwd=2,add = TRUE)
# legend("bottomright", legend = c("Max", "Min", "Median"), col = c("blue", "green", "red"),lwd = 1,fill = 1:3, ncol = 1)
# lty 控制线的类型，col控制线的粗细
#legend("bottomright", legend=c("empirical"),col=c("blue"), lwd=c(2))
#dev.off()
#rocobj <- roc(as.factor(label1), as.vector(svm1))
#auc(rocobj)
#rets <-c("threshold","sensitivity","specificity","accuracy","ppv","tn","fp")
#show1 <-ci.coords(rocobj,x="best",input = "threshold",best.policy='random',ret=rets)


if(F)
{
  label1<-matrix(c(rep(0,nrow(norma1)),rep(1,nrow(disea))))
  mdata1<-rbind(norma1,disea)
  svm1<-NULL
  # 使用留1法绘制ROC curve
  # 每次都排除当前的i的值
  for(i in 1:nrow(label1))
  {
    mtrain  <-label1[-i,]
    predata <-mdata1[i,]
    traindata <-mdata1[-i,]
    model  <-svm(traindata,mtrain,type="C-classification",kernel="radial")
    pred  <-predict(model,t(predata),decision.values=TRUE)
    f     <-attr(pred,"decision.values")
    prb1   <-1/(1+exp(-f))
    svm1   <-rbind(svm1,prb1)
    #prb1<-cbind(svm,m1)
  }
}