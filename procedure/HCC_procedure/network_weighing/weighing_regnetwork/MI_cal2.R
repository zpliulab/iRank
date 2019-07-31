###############################
##搭建MI calculation 程序框架
###############################
rm(list = ls())
setwd("D:/git_data/procedure/HCC_procedure/network_weighing/weighing_regnetwork")
##读取RNA-sq data and Gene regulation network
##################
##RNA-seq unique data 
R_da <-as.matrix(read.table("RNA_seq_unique.txt",header=T,sep="\t"))

mode(R_da)
dim(R_da)
R_na <-(as.matrix(R_da[,1]))
dR <-R_da[,2:dim(R_da)[2]]
Tse <-seq(1,73,by=2)#Tumor data
Nse <-seq(2,74,by=2)#Normal data
mode(R_na)
dim(R_na)
##读取Gene regulatory network data 
RT <-as.matrix(read.table("Regnetwork_in_RNA_10_2.txt",header=T,sep="\t"))
dim(RT)
#SR <-RT[,1]
#TR <-RT[,3]
#LL <-length(union(unique(SR),unique(TR)))
##17020
##加载infotheo package 
##用于MI 的计算
##############################
####采用离散化整体，之后逐个取出来
#############################
library(infotheo)
dR1 <-discretize(dR)#整体离散化数据
dR1[1:10,]
DR_data1 <-dR1[,Tse]
NR_data1 <-dR1[,Nse]
DMI1 <-rbind()
for(i in 1:dim(RT)[1])
{
 o1 <-which(RT[i,1]==R_na)##source node
 #o1
 D11<-t(as.matrix(DR_data1[o1,]))
 C11<-t(as.matrix(NR_data1[o1,]))
 
 o2 <-which(RT[i,3]==R_na)## target node
 #o2
 D21<-t(as.matrix(DR_data1[o2,]))
 C21<-t(as.matrix(NR_data1[o2,]))
 
 DM1 <-abs(mutinformation(D11,D21,method="emp") - mutinformation(C11,C21,method="emp"))
 DMI1<-rbind(DMI1,DM1)#D1: row :samples, column: feature,故需要转置一下
}

length(DMI1)
length(which(DMI1>0))
length(unique(DMI1))

write.table(cbind(RT,DMI1),file="Reg_in_RNA_MI_C1.txt",sep='\t',row.names=F,quote=F)

