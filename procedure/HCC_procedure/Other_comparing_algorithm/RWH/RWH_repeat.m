clear;
close all;clc;
% RegNetwork的信息
Reg  =importdata('Regnetwork_only_gene.txt');
R_s=Reg.textdata(2:end,2);
R_s_num=[];
for(j=1:length(R_s))
    R_s_num(j,1)=str2num(R_s{j});
end
R_t =Reg.data;
R_node=unique(union(R_s_num,R_t));
[m,n]=find(isnan(R_node));
 R_node(m,:)=[];
% 构建邻接矩阵
R_matrix=sparse(length(R_node),length(R_node));
for(i=1:length(R_t))
    l1=find(R_s_num(i)==R_node);
    l2=find(R_t(i)==R_node);
    if(~isempty(l1))
        if(~isempty(l2))
            R_matrix(l2,l1)=1;
        end
    end
end
%矩阵正规化
s_R=sum(R_matrix,1);%按列求和
for(jj =1:length(s_R))
    if(s_R(jj)>0)
        R_matrix(:,jj)=R_matrix(:,jj)/s_R(jj);
    else
        R_marix(:,jj)=1/length(R_node);
    end
end
%
%读取ppi 网络的信息
pp1=xlsread('my_pair_10_22.csv');
pp1s=pp1(:,1);
pp1t=pp1(:,3);
pp2=xlsread('string_t_905_10_2.csv');
pp2s=pp2(:,1);
pp2t=pp2(:,3);
P_node=unique((union(union(pp1s,pp1t),union(pp2s,pp2t))));
[m,n]=find(isnan(P_node));
P_node(m,:)=[];
P_matrix=sparse(length(P_node),length(P_node));  
for(i=1:length(pp1s))
    l1=find(pp1s(i) == P_node);
    l2=find(pp1t(i) == P_node);
    if(~isempty(l1))
        if(~isempty(l2))
           P_matrix(l1,l2)=1;
           P_matrix(l2,l1)=1;
        end
    end
end
for(i=1:length(pp2s))
    l1=find(pp2s(i) ==P_node);
    l2=find(pp2t(i) ==P_node);
    if(~isempty(l1))
        if(~isempty(l2))
            P_matrix(l1,l2)=1;
            P_matrix(l2,l1)=1;
        end
    end
end
% 矩阵正规化
sum_P=sum(P_matrix,1);
for(i=1:length(sum_P))
    if(sum_P(i)>0)
        P_matrix(:,i)=P_matrix(:,i)/sum_P(i);
    else
        P_matrix(:,i)=1/length(P_node);
    end
end
% 构建Reg--Pro 的转移关系
RP=sparse(length(P_node),length(R_node));
ww=intersect(P_node,R_node);
for(i=1:length(ww))
    l1=find(ww(i)==R_node);
    l2=find(ww(i)==P_node);
    RP(l2,l1)=1;
end
% Reg--Pro 矩阵正规化
lamda=0.5;
sum_RP=sum(RP,1);
for(i=1:length(sum_RP))
    if(sum_RP(i)>0)
        RP(:,i)=lamda*RP(:,i)/sum_RP(i);
    end
end
% Regnetwork 
for i=1:length(ww)
    l1=find(ww(i)==R_node);
    if(~isempty(l1))
        R_matrix(:,l1)=(1-lamda)*R_matrix(:,l1);
    end
end
%% 构建转移矩阵
Trans=sparse((length(P_node)+length(R_node)),(length(P_node)+length(R_node)));
% 逐一赋值
Trans((1:length(P_node)),(1:length(P_node)))=P_matrix;
Trans((length(P_node)+1):((length(P_node)+length(R_node))),((length(P_node)+1):(length(P_node)+length(R_node))))=R_matrix;
Trans(1:length(P_node),(length(P_node)+1):(length(P_node)+length(R_node)))=RP;
%
r=0.85;%0.7,0.6,0.5,0.4,0.8
threshold=1e-10;
N = length(Trans);
PR = 1/N*ones(N,1);
iter = 1;
delta_PR = Inf; 
while (delta_PR > threshold)    %iterate until error is less than 1e-10 1e-7
    tic;
    prev_PR = PR;               %save previous PageRank vector (t-1)
    PR = (1-r)*Trans*PR + r*PR;     %calculate new PageRank (t)
    delta_PR = norm(PR-prev_PR);%calculate new error
    t(iter)=toc;
    iter = iter + 1;
end
% 排序
% 提取排序后的结果
% protein 的结果
result_P=PR(1:length(P_node));
rank_p=sort(unique(result_P),'descend');
% regnetwork的结果
result_R=PR((length(P_node)+1):end);
rank_r=sort(unique(result_R),'descend');
%rank aggregation
% 两者相比，取更小的
all_node=unique(union(P_node,R_node));
rank_all=[];
pr_all=[];
for(i=1:length(all_node))
    l1=find(all_node(i)==P_node);
    l2=find(all_node(i)==R_node);
    if(~isempty(l1))
        if(isempty(l2))
            rank_all(i,1)=find(result_P(l1)==rank_p);
            pr_all(i,1)=result_P(l1);
        else
            k1=find(result_P(l1)==rank_p);
            k2=find(result_R(l2)==rank_r);
            rank_all(i,1)=max([k1,k2]);
            if rank_all(i,1)==k1
                pr_all(i,1)=result_P(l1);
            else
                pr_all(i,1)=result_R(l2);
            end
        end
    else
         rank_all(i,1)=find(result_R(l2)==rank_r);
         pr_all(i,1)=result_R(l2);
    end     
end

% gain information of disease genes
disease = importdata('hcc_disease_genes.txt');
HD = disease.data;
D=[]; ii=1;
for(i=1:length(HD))
    l1=find(HD(i)==all_node);
    if(~isempty(l1))
        D(ii,1)=HD(i);
        D(ii,2)=pr_all(l1);
        D(ii,3)=rank_all(l1);
        ii=ii+1;
    end
 end
Auc=[];
Normal_rank=[];
Normal_pr=[];
Normal_id=[];
Circle=1000;
normal_genes_all=all_node;%intersect((Reg_node-0.4),(S_num-0.3));
Candidate_normal=setdiff(normal_genes_all,D(:,1));

for(count=1:Circle)
Rand_loc_p=randperm(length(Candidate_normal),length(D));
Rand_gene_ID=Candidate_normal(Rand_loc_p);
nor=[];ii=1;
for(i=1:length(Rand_loc_p))
    l1=find((Rand_gene_ID(i))==all_node);
    if(~isempty(l1))
        nor(ii,1)=Rand_gene_ID(i);
        nor(ii,2)=pr_all(l1);
        nor(ii,3)=rank_all(l1);
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
%之前贴错标签了
%正常是0，异常是1
mdlSVM = fitcsvm(data_p,y_p,'Standardize',true);
mdlSVM = fitPosterior(mdlSVM);
[~,score_svm] = resubPredict(mdlSVM);
[Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(y_p,score_svm(:,mdlSVM.ClassNames),'true');
Auc(count)=AUCsvm;
end
boxplot(Auc)
max(Auc)
mean(Auc)
min(Auc)
length(find(Auc>0.7))
length(find(Auc<0.6))
std(Auc)
median(Auc)
save('RWH_new.mat','D','Auc','Normal_id','Normal_pr','Normal_rank','-v6')





