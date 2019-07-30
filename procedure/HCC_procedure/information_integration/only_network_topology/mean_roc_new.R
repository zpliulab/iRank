rm(list=ls())
library(R.matlab)#D
#setwd("D:/shx_bioinformatics/TCGA_Project/used_for_plot_ROC_and_static")
#data <- readMat("reg_s_33_new.mat")
setwd("D:/shx_bioinformatics/TCGA_Project/Copy_of_Reg_original")
#data <- readMat("RDC_3_33_new.mat")
#setwd("I:\\shx_bioinformatics\\TCGA_Project\\Copy_of_reg_c")
data <- readMat("ORI_disease_33_new.mat")
##直接把pr value 输入看看可否出来roc
D <-as.data.frame(data[2])
N_pr <-(data[["Normal.pr"]])
D_pr <-(D[,"D.2"])
index_all <-rbind()
library(cutpointr)                                  
#library(pROC)
#tp_list <-rbind()
#fp_list <-rbind()
#Thre_list <-rbind()
for(i in 1:dim(N_pr)[2])#
{
  #i <-1
  data_t <- as.matrix(D_pr)
  #length(data_t)
  data_n <- as.matrix(N_pr[,i])
  label<-matrix(c(rep(1,nrow(data_t)),rep(0,nrow(data_n))))
  ##疾病标1正常标0
  #label
  use <-c(data_t,data_n)
  # k1 <-roc(as.factor(label),as.vector(use),plot=T,print.auc=T,print.thres=T,col="blue",lty=1,lwd=2)
  #ll1 <-t(as.matrix(k1[["specificities"]]))
  #ll1
  #length(ll1)
  #dim(ll1)
  #SE_list <-rbind(SE_list,t(as.matrix(k1[["specificities"]])))
  #SP_list <-rbind(SP_list,t(as.matrix(k1[["sensitivities"]])))
  #Thre_list <-rbind(Thre_list,t(as.matrix(k1[["thresholds"]])))
  #k1[["specificities"]]
  opt_cu <- cutpointr(x =as.vector(use) , class = as.vector(label))
  opt_cu
  zz <-opt_cu[["roc_curve"]]
  z1 <-zz[[1]]
  #mode(z1)
  #tp_list <-rbind(tp_list,t(as.matrix(z1[["tpr"]])))
  #fp_list <-rbind(fp_list,t(as.matrix(z1[["fpr"]])))
  acc <-round(opt_cu[["acc"]],4)
  se <-round(opt_cu[["specificity"]],4)
  sp <-round(opt_cu[["sensitivity"]],4)
  F1 <-round(2*(se*sp)/(sp+se),4)
  auc <-round(opt_cu[["AUC"]],4)
  used_index <-cbind(se,sp,acc,F1,auc)
  index_all<-rbind(index_all,used_index)
}
dim(index_all)
sd(index_all)
mm <-cbind()
sdd <-cbind()
for(i in 1:dim(index_all)[2])
{
   mm <-cbind(mm,round(mean(index_all[,i]),4))
  sdd <-cbind(sdd,round(sd(index_all[,i]),4))
}
info <-c("Se","Sp","acc","F1","AUC")
info
mm
sdd

if(F)
{
#library(rapport)
va <-rbind()
xx <-rbind()
mm <-rbind()
mn <-rbind()
#lk <-unique(Thre_list)
sd <-seq(0,1,by=1/60)
for(i in 1:length(sd))
{
  #画格子 
  #落在0~0.001之间的所有点的se的均值
  #i <-0
  # i <-1
  l1 <-which(fp_list<sd[i+1])
  l2 <-which(fp_list>=sd[i])
  in1 <-intersect(l1,l2)
  ##若是非空
  if(length(in1)>0)
  {
  mm1 <-max(tp_list[in1])
  mn1 <-min(tp_list[in1])
  #va <-rbind(va,mean(SE_list[l1]))
  mm <-rbind(mm,mm1)
  mn <-rbind(mn,mn1)
  va <-rbind(va,(mn1+(mm1-mn1)/2))
  xx <-rbind(xx,sd[i])
  }
}
length(xx)
length(va)
plot(c(0,xx),c(0,va),type = 'l',col = "red",lty=4,ylim = c(0.0,1),xlab = "False positive rate",
     ylab = "True Positive rate")
library(ggplot2)
x11 <-c(0,xx)
yy1 <-c(0,va)
yy2 <-c(0,mn)
yy3 <-c(0,mm)
data <-data.frame(x11,yy1,yy2,yy3)
ggplot(data,aes(x=x11,y=yy1))+
  geom_line(colour = "blue",size=1)+
  geom_ribbon(aes(ymin=yy2,ymax=yy3),alpha=0.2,fill="red")+
  #labs(x="1 - Specificity", y="Sensitivity", colour="")+,axis.title.x=element_blank()
  xlab("1-Specificity")+theme(axis.text.x = element_text(colour = 'black', angle = 0, size = 13, hjust = 0.5, vjust = 0.5)) + 
  ylab("Sensitivity") + theme(axis.text.y = element_text(colour = 'black', size = 12), axis.title.y = element_text(size = 12, hjust = 0.5, vjust = 0.2)) + 
  theme(strip.text.y = element_text(size = 11, hjust = 0.0, vjust = 0.2, face = 'bold'))
##这个图可以了
}