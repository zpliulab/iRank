# -*- coding: cp936 -*-
#import os
#os.chdir("D:\Cygwin\Cygwin_load\home\sea_sunshine\PHENORank\phenorank\data_phenorank")
import cPickle
import pkg_resources
dir_data="data_phenorank"
#con = open(pkg_resources.resource_filename("phenorank", dir_data + "/gc_h.pickle"), "rb")
#con = open(pkg_resources.resource_filename("phenorank", dir_data + "/gc_m.pickle"), "rb")
#con = open(pkg_resources.resource_filename("phenorank", dir_data + "/pheno_condition_ic_matrix.pickle"), "rb")
#con = open(pkg_resources.resource_filename("phenorank", dir_data + "/pheno_cooccur.pickle"), "rb")

#con1 =open("gc_h.pickle","rb")
gc_h = cPickle.load(con)#Compressed Sparse Row matrix
gc_h
# 把数据注意读取进去
# unique 一下
# 读取 biogrid_hiii14_hprd_intact.tsv

'''f1=open('biogrid_hiii14_hprd_intact.tsv','r')
dat_1=[]
dat1=[]
dat2=[]
for dd in f1:
    dat_1=dd.split(',')
    dat1.append(dat_1[0])
  #  dat2.append(dat_1[1])
print(len(dat_1))
#print(len(dat_2))
print(len(list(set(dat_1))))'''
'''
#print(len(list(set(dat_2))))
# 使用pandas
import os
os.chdir("D:\\Cygwin\\Cygwin_load\\home\\sea_sunshine\\PHENORank\\phenorank\\data_phenorank")
import pandas as pd
#import pkg_resources
#biogrid_hiii14_hprd_intact
#condition_ic
#cp_h_omim
#cp_h_omim_ids
#cp_m
#cp_m_ids
#gc_h_omim
#gc_h_omim_ids
#gc_m
#gc_m_ids
#phenorank_conditions
#phenorank_genes
#phenorank_phenotypes
#phenotype_ancestors
#phenotype_ic
data1 =pd.read_csv('phenorank_genes.tsv',sep='\t')
#data1[1,1]
type(data1)
#data1.loc[[1],[1]]
print(data1.iloc[0:5,1])
Souce  =len(list(set(data1.iloc[:,0])))
Target =len(list(set(data1.iloc[:,1])))
l1 =data1.iloc[:,0]
l2 =data1.iloc[:,1]
allgene =len(list(set(l1.append(l2))))

#head(1)
#con = open(pkg_resources.resource_filename("phenorank", dir_data + "/condition_ic.tsv"), "r")
#omim_ic = np.array(list(pd.read_csv(con, header=None)[0]))'''

