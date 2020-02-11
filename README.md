# iRank
## a network-based method for ranking cancer genes by integrative ranks from multi-omics data
## 
## data:
### Omcis data:
#### Hepatocellular carcinoma (HCC) data from TCGA database
#### including DNA methylation, somatic mutation, copy number variant, mRNA-seq, miRNA-seq from 37 samples with both tumor and normal samples
#### Prostate adenocarcinoma (PRAD) data from TCGA database
#### including DNA methylaion, copy number variant, mRNA-seq, miRNA-seq from 35 samples with both tumor and normal samples
### Network
#### Gene regulation network (GRN), from RegNetwork
#### Protein-protein interaction network (PPIN), from several databases 
## Procedure:
### (a) prepare data: 
####    Mutual information calculation 
### (b) implement iRank:
####    Integrating different omics data, ORI, ORIr, ORIrd,ORIrds, ORIrdms..... 
####     Saving results as .mat ('-v6')
### (c) results show:
####    Input .mat, implement .R, get the ROC curves and the boxplots of different integrations
###
## Giude:
#### iRANK_guide gives an example on how to implement our iRANK (ORIrd), a .txt guide in this file can guide you implement our procedure easily.
