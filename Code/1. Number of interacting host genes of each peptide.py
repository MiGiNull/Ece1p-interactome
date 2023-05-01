#!/usr/bin/env python
# coding: utf-8

# ### 1.Counting interacting genes of each peptide

# In[1]:


import  pandas  as pd
import numpy as np


# ##### Read by specifying the form name

# In[34]:


df=pd.read_excel('/Users/cwn/Desktop/待处理/LNN/20220407 Final list of genes for all peptides.xlsx',sheet_name=0)#You can specify the form to be read by sheet_name, or you can specify the form by value, '0' is the first form
data=df.head()#Read the first 5 rows of data by default
print("获取到所有的值:\n{0}".format(data))#Formatted output


# #### pandas manipulation of Excel rows

# In[3]:


print("输出列标题",df.columns.values)
genes=df['BLAST genes'].values#Read the value of the genes column for all rows
#print("输出值\n",list(genes))

#Read the interacting genes for each sheet (peptide)
def read_sheet(filename,sheet):
    df=pd.read_excel(filename,sheet_name=sheet)
    gene_list=list(df['BLAST genes'].values)
    return gene_list


# #### Count the number of occurrences of each gene

# In[4]:


gene_count = pd.value_counts(genes)
print(gene_count)

#Splitting rows with multiple genes
def gene_count(gene_list):
    genes=[]
    for i in gene_list:
        i=i.split(",")
        for gene in i:
            if gene.strip() not in genes:
                genes.append(gene.strip())
    return genes


# #### Final: counting for each peptide

# In[35]:


filename='Data to be analysed'
tox1=gene_count(read_sheet(filename,0))
tox2=gene_count(read_sheet(filename,1))
tox3=gene_count(read_sheet(filename,2))
tox4=gene_count(read_sheet(filename,3))
tox5=gene_count(read_sheet(filename,4))
tox6=gene_count(read_sheet(filename,5))
tox7=gene_count(read_sheet(filename,6))
tox8=gene_count(read_sheet(filename,7))
#counting interacting genes for each peptide
data=[len(tox8),len(tox7),len(tox6),len(tox5),len(tox4),len(tox3),len(tox2),len(tox1)]
data
#data=[len(FL),len(tox3)]


# In[36]:


tox=''
for i in read_sheet(filename,7):
    tox= tox+"'"+i+"',"
print(tox)


# #### Plotting bar graphs

# In[51]:


import matplotlib.pyplot as plt
#labels = ['FL','Tox8', 'Tox7', 'Tox6', 'Tox5', 'Tox4','Tox3','Tox2','Tox1']

plt.figure(figsize=(4,1), dpi=300)  ##figure size
labels = ['FL','Tox3']
plt.barh(range(len(data)), data, tick_label=labels,ec='k',  lw=0.7,color="cornflowerblue",height=0.7)##width
plt.margins(x=0.03)

#Border
ax = plt.gca()
bwith=1
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)

#Set the size of the scale value and the font of the scale value
plt.tick_params(labelsize=13)
labels = ax.get_xticklabels() + ax.get_yticklabels()
#[label.set_fontname('Times New Roman') for label in labels]
 
#set names for axis
font = {
'weight' : 'normal',
'size'   : 13,
}

#set labels for x and y axis  
plt.xlabel("Number of interacting host genes",font)  

#Save and plotting
plt.savefig('File name',bbox_inches = 'tight')
plt.show()


# In[37]:


print(FL)
string=''
for i in FL:
    string+=i.strip()+","
print(string)


# In[38]:


print(tox8)
string=''
for i in tox8:
    string+=i.strip()+","
print(string)


# In[39]:


print(tox7)
string=''
for i in tox7:
    string+=i.strip()+","
print(string)


# In[40]:


print(tox6)
string=''
for i in tox6:
    string+=i.strip()+","
print(string)


# In[41]:


print(tox5)
string=''
for i in tox5:
    string+=i.strip()+","
print(string)


# In[42]:


print(tox4)
string=''
for i in tox4:
    string+=i.strip()+","
print(string)


# In[43]:


print(tox3)
string=''
for i in tox3:
    string+=i.strip()+","
print(string)


# In[44]:


string=''
for i in tox2:
    string+=i.strip()+","
print(string)


# In[45]:


string=''
for i in tox1:
    string+=i.strip()+","
print(string)


# In[30]:


print("CDKN2C|CRYBA1|DMTN|SFN|PDE4D|PPP1R8|PPP2CA|PPP2R1A|PRKAA2|QARS1|RB1|TRIM27|SMARCB1|CTDSP2|CTDSPL|PNKP|KDM1A|VPS28|PRKAG2|INPP5K|OTUB1|SMYD3|HEXIM2|ZNF675|SPRED2".replace("|"," "))


# In[1]:


version


# In[ ]:




