#!/usr/bin/env python
# coding: utf-8

# ### 将基因与肽段的互作表输出为矩阵

# In[5]:


import  pandas  as pd
import numpy as np


# #### 读取每个肽段和每个基因互作的强度

# In[13]:


#读取基因互作强度表单
#将有多个基因的行拆开
def gene_count(gene_list,interaction):
    genes={}
    for i in range(len(gene_list)):
        spilt_genes=gene_list[i].split(",")
        if len(spilt_genes)>1:
             for gene in spilt_genes:
                    if gene.strip() not in genes.keys():
                        genes[gene.strip()]=round(interaction[i],3)
        else:
            if spilt_genes[0].strip() not in genes.keys():
                genes[spilt_genes[0].strip()]=round(interaction[i],3)
    return genes

#读取每个sheet（peptide）的互作基因和其互作强度
def read_gene_interation(filename,sheet):
    df=pd.read_excel(filename,sheet_name=sheet,usecols=['BLAST genes','Strength of interaction'])
    gene_list=list(df['BLAST genes'].values)
    interaction=list(df['Strength of interaction'].values)
    gene_interaction=gene_count(gene_list,interaction)
    return gene_interaction



filename='Data'
tox1=read_gene_interation(filename,0)
tox2=read_gene_interation(filename,1)
tox3=read_gene_interation(filename,2)
tox4=read_gene_interation(filename,3)
tox5=read_gene_interation(filename,4)
tox6=read_gene_interation(filename,5)
tox7=read_gene_interation(filename,6)
tox8=read_gene_interation(filename,7)

#合并字典为dataframe
final={'Ece1-I':tox1,'Ece1-II':tox2,'Ece1-III':tox3,'Ece1-IV':tox4,'Ece1-V':tox5,'Ece1-VI':tox6,'Ece1-VII':tox7,'Ece1-VIII':tox8}
final=pd.DataFrame(final)
final=final.fillna(0)
print(final)


# #### 输出矩阵，在R里面做overlap的分析（超级几何检验）

# In[15]:


final.to_csv(r'Filename')## 如不输出位置：index=false


# 超几何分布在基因overlap中的使用：<https://blog.csdn.net/linkequa/article/details/86491665>
