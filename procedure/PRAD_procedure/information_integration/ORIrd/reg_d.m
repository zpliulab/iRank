clear;clc;close all;
%% gene level 的信息
% 加载RegNetwork 信息
G_data=importdata('Reg_intersect.txt');
used_G=G_data.data(:,:);
Sg=unique(used_G(:,1));
Tg=unique(used_G(:,2));
Go_node=unique(union(Sg,Tg));
%% Reg level Transition Matrix Construction
Reg_node=unique(Go_node);
NR=length(Reg_node);
R=sparse(NR,NR);
RW=sparse(NR,NR);
RCS=sparse(NR,NR);
% 加入gene regulation 信息
for(i=1:length(used_G))
    l1=find(used_G(i,1)==Reg_node);
    l2=find(used_G(i,2)==Reg_node);
    if(~isempty(l1))
        if(~isempty(l2))
            R(l2,l1)=1;
            RW(l2,l1)=used_G(i,3);
            RCS(l2,l1)=1;
        end
    end
end

% 转移矩阵正规化
sum_R=sum(R,1);
for(i=1:length(sum_R))
    if(sum_R(i)>0)
    R(:,i)= R(:,i)/sum_R(i);
    end
end
% 转移矩阵加权
RTW=R.*RW;
%%  加载DNA methylation 信息
D_data=importdata('DNA_intersect.txt');
used_D=D_data.data(:,:);
% % 创建DNA--->RegNetwork
lamda=0.1;
%  存在Reg_node
D_num=(intersect(unique(used_D(:,2)-0.3),(Reg_node-0.3)));
DR=sparse(NR,length(D_num));
DWR=sparse(NR,length(D_num));
DRCS=sparse(NR,length(D_num));
for(i=1:length(used_D))
    l1=find((used_D(i,1)-0.2)==D_num);
    l2=find(used_D(i,2)==Reg_node);
    if(~isempty(l1))
        if(~isempty(l2))
           DR(l2,l1)=1;
           DWR(l2,l1)=used_D(i,3);
           DRCS(l2,l1)=1;
        end
    end
end
% 转移矩阵正规化
sum_D=sum(DR,1);
for(i=1:length(sum_D))
    if(sum_D(i)>0)
        DR(:,i)=DR(:,i)/sum_D(i);
    end
end
% 转移矩阵增加权重
DRTW=DR.*DWR*(lamda/2);
%% 构建总体的转移矩阵
Node_all=[Reg_node;D_num];
number=length(Node_all);
Trans=sparse(number,number);
Trans(1:NR,1:NR)=RTW;
Trans((1:NR),(NR+1):(NR+length(D_num)))=DRTW;

CS=sparse(number,number);
CS(1:NR,1:NR)=RCS;
CS((1:NR),(NR+1):(NR+length(D_num)))=DRCS;
%%  输入程序
r=0.85;
threshold=1e-10;%
N = length(Trans);
PR =( 1/N*ones(N,1));%
restart =PR;
iter = 1;
delta_PR = Inf; 
while (delta_PR > threshold || iter>200)    %iterate until error is less than 1e-10 1e-7
    tic;
    prev_PR = PR;               %save previous PageRank vector (t-1)
    CST=CS.*(1*Trans);%*1/N Trns'
    CST(find(isnan(CST)==1))=0;
    PR = r*CST* PR + (1-r)*restart;     %calculate new Pa
    delta_PR= norm(PR-prev_PR);%calculate new error
    t(iter)=toc;
    iter = iter + 1;
end
length(unique(PR))
[Rank,index]=sort(unique(PR'),'descend');%
rank=1:length(unique(PR));
%% gene 在网络 1中单独的排序 RegNetwork 中的排序
[LRPR,LRrank]=sort(unique(PR(1:NR)),'descend');
length(unique(LRPR))
%% gene 在网络2中单独的排序  DNA methylation中的排序
[LDPR,LDrank]=sort(unique(PR((NR+1):(NR+length(D_num)))),'descend');
length(unique(LDPR))

HD =importdata('29.txt');
D=[]; ii=1;
for(i=1:length(HD))
    l1=find((HD(i)+0.3)==Reg_node);
    l2=find((HD(i))==D_num);
    if(~isempty(l1))
            if(~isempty(l2))
               D(ii,1)=HD(i);
               D(ii,2)=PR(l1);
               D(ii,3)=find(PR(l1)==LRPR);
               ii=ii+1;
             end
    end
end
Auc=[];
Normal_rank=[];
Normal_pr=[];
Normal_id=[];
Circle=1000;
normal_genes_all=intersect((D_num),(Reg_node-0.3));
Candidate_normal=setdiff(normal_genes_all,D(:,1));

for(count=1:Circle)
Rand_loc_p=randperm(length(Candidate_normal),length(D));
Rand_gene_ID=Candidate_normal(Rand_loc_p);
nor=[];ii=1;
for(i=1:length(Rand_loc_p))
    l1=find((Rand_gene_ID(i)+0.3)==Reg_node);
    if(~isempty(l1))
        nor(ii,1)=Rand_gene_ID(i);
        nor(ii,2)=PR(l1);
        nor(ii,3)=find(PR(l1)==LRPR);
        ii=ii+1;
    end
end
% 存储一下rank
for(j=1:length(nor))
Normal_rank(j,count)=nor(j,3);
Normal_pr(j,count)=nor(j,2);
Normal_id(j,count)=nor(j,1);
end
nor_rank_f=nor(:,3);
D_rank_f=D(:,3);
data_p=[];
for(i=1:(2*(length(D_rank_f))))
    if(i<=length(D_rank_f))
   data_p(i,1) =D_rank_f(i);%disease 
    else
       data_p(i,1) = nor_rank_f(i-length(D_rank_f)); % normal 
    end
end
y_p=(1:(2*(length(D_rank_f))))'<(length(D_rank_f));%>
fe =[D(:,2)',nor(:,2)'];
[Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(y_p,fe,'true');
Auc(count)=AUCsvm;
end
boxplot(Auc)
max(Auc)
mean(Auc)
min(Auc)
std(Auc)
save('reg_d_29.mat','D','Auc','Normal_id','Normal_pr','Normal_rank','-v6')
median(Auc)
