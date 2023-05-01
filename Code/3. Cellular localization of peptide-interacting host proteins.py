#!/usr/bin/env python
# coding: utf-8

# ### 1.统计各个肽段的基因亚细胞定位

# In[2]:


import  pandas  as pd
import numpy as np


# #### 读取 表单

# In[6]:


#读取metascape_result.xlsx的"注释"表单
#df=pd.read_excel('/Users/cwn/workdir/project/peptide-host proteins interction/result/all/metascape_result.xlsx',sheet_name=0)#可以通过sheet_name来指定读取的表单,也可以通过数值指定表单，『0』则为第一个表单
#df=pd.read_excel('/Users/cwn/workdir/project/peptide-host proteins interction/result/FL-tox3/metascape_result.xlsx',sheet_name=0)
df=pd.read_excel('Path to data',sheet_name=0)

data=df.head()#默认读取前5行的数据
print("获取到所有的值:\n{0}".format(data))#格式化输出


# #### pandas操作Excel的行列，读取注释表单中的细胞亚定位注释：'Subcellular Location (Protein Atlas)'

# In[7]:


print("输出列标题",df.columns.values)
location=list(df['Subcellular Location (Protein Atlas)'].values)#读所有行的genes列的值
#print("输出值\n",gene_list)
#gene_list=list(df["Tox1"].values)
#print([i for i,x in enumerate(gene_list) if x==1])
#读取每个peptide的互作基因,输出互作基因的亚细胞定位
def read_gene(peptide,annotation_location):
    gene_list=list(df[peptide].values)
    position=[i for i,x in enumerate(gene_list) if x==1]##列表中某值的位置，这里是指该肽段的互作基因的位置，以找到肽段的互作基因是哪些
    location=[annotation_location[i] for i in range(len(annotation_location)) if i in position]
    return location
#print(read_gene("Tox1",location))


# #### 统计每个细胞定位的比例

# In[8]:


##计算nan值，为未知(unknow)定位
def count_unknow(peptide,annotation_location):
    a=read_gene(peptide,annotation_location)
    a_no_nan = [a_ for a_ in a if a_ == a_]
    nan=len(a)-len(a_no_nan)
    return nan

##对定位目标名称进行简化并统计
def count_location(peptide,annotation_location):
    a=read_gene(peptide,annotation_location)
    a_no_nan = [a_ for a_ in a if a_ == a_]
    clear_location1=[]
    for i in a_no_nan:
        if ";" in i:
            x=i.split(";")
            for n in x:
                clear_location1.append(n.strip())
        else:
            clear_location1.append(i)
        
    clear_location2=[]
    for i in clear_location1:
        if "(" in i:
            x=i.split("(")
            clear_location2.append(x[0].strip())
        else:
            clear_location2.append(i)

    clear_location3=[]
    for i in clear_location2:
        if ":" in i:
            x=i.split(":")
            clear_location3.append(x[1].strip())
        else:
            clear_location3.append(i)
    ##名称分组
    nucleus=["Nucleoplasm","Nuclear bodies","Nucleoli","Nuclear speckles","Nuclear membrane","Nucleoli fibrillar center"]
    membranes=["Plasma membrane","Vesicles"]
    others=["Actin filaments","Centrosome","Microtubules","Focal adhesion sites","Intermediate filaments","Peroxisomes","Microtubule ends","Cell Junctions","Midbody ring",
           "Endosomes","Lysosomes","Centriolar satellite","Cytokinetic bridge"]
    cytoplasm=["Cytosol","Cytoplasmic bodies"]
    golgi=["Golgi apparatus"]
    ER=["Endoplasmic reticulum"]
    for i in range(len(clear_location3)):
        if clear_location3[i] in nucleus:
            clear_location3[i]="Nucleus"
        elif clear_location3[i] in  membranes:
            clear_location3[i]="Membranes"
        elif clear_location3[i] in others:
            clear_location3[i]="Others"
        elif clear_location3[i] in cytoplasm:
            clear_location3[i]="Cytoplasm"
        elif clear_location3[i] in golgi:
            clear_location3[i]="Golgi"
        elif clear_location3[i] in ER:
            clear_location3[i]="ER"   
    ##计数
    location_count=dict(pd.value_counts(clear_location3))
    location_count["Unknow"]=count_unknow(peptide,annotation_location)
   # location_count=pd.DataFrame(location_count)
    return location_count
#print(pd.value_counts(a))
#print("\n",nan,"\n")
#print(clear_location2)
#print(clear_location3)
#print(len(clear_location3))


# #### 最终：统计每个肽段

# In[12]:


tox1=count_location("Ece_I",location)
tox2=count_location("Ece_II",location)
tox3=count_location("Ece_III",location)
tox4=count_location("Ece_IV",location)
tox5=count_location("Ece_V",location)
tox6=count_location("Ece_VI",location)
tox7=count_location("Ece_VII",location)
tox8=count_location("Ece_VIII",location)
#fl=count_location("FL",location)
##将每个肽段的细胞定位数值合并为df
final={'Ece-I':tox1,'Ece-II':tox2,'Ece-III':tox3,'Ece-IV':tox4,'Ece-V':tox5,'Ece-VI':tox6,'Ece-VII':tox7,'Ece-VIII':tox8}
#final={'Tox3':tox3,'FL':fl}
final=pd.DataFrame(final)

##计算比例
final=final.T
final=final.div(final.sum(axis=1), axis=0)
final=final.fillna(0)
print(final)
print(tox1)
print(tox2)
print(tox3)
print(tox4)
print(tox5)
print(tox6)
print(tox7)
print(tox8)


# #### 绘制堆叠图

# In[16]:


import matplotlib.pyplot as plt

##获取每一组的值
Nucleus=np.array(final['Nucleus'].values.tolist())
Cytoplasm=np.array(final['Cytoplasm'].values.tolist())
Membranes=np.array(final['Membranes'].values.tolist())
Others=final['Others'].values.tolist()
Golgi=final['Golgi'].values.tolist()
Mitochondria=final['Mitochondria'].values.tolist()
ER=np.array(final['ER'].values.tolist())
Unknown=final['Unknow'].values.tolist()

plt.figure(figsize=(7,6), dpi=300)  ##通过设置图片大小来调整条形图间距

labels = ['Ece-I','Ece-II','Ece-III', 'Ece-IV', 'Ece-V', 'Ece-VI', 'Ece-VII','Ece-VIII']
plt.bar(labels, Nucleus, color='#ABB584', label='Nucleus')
plt.bar(labels, Cytoplasm, bottom=Nucleus, color='#C2959A', label='Cytoplasm')      #横向的堆叠图：bar——barh；bottom——left
plt.bar(labels, Membranes, bottom=Nucleus+Cytoplasm, color='#7E9E9B', label='Membranes')
plt.bar(labels, Others, bottom=Nucleus+Cytoplasm+Membranes, color='#E6D9C1', label='Others')
plt.bar(labels, Golgi, bottom=Nucleus+Cytoplasm+Membranes+Others, color='#B7C0A4', label='Golgi')
plt.bar(labels, Mitochondria, bottom=Nucleus+Cytoplasm+Membranes+Others+Golgi, color='#FFD2D2', label='Mitochondria')
plt.bar(labels, ER, bottom=Nucleus+Cytoplasm+Membranes+Others+Golgi+Mitochondria, color='#969393', label='ER')
plt.bar(labels, Unknown, bottom=Nucleus+Cytoplasm+Membranes+Others+Golgi+Mitochondria+ER, color='#E0B192', label='Unknown')
#plt.title("Cellular localization ")                                                        #图片标题
plt.legend(loc=[1.06,0.3])                                                             #图例的显示位置设置

#设置边框
ax = plt.gca()#获取边框
bwith=1
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)

#设置坐标刻度值的大小以及刻度值的字体
plt.tick_params(labelsize=13)
labels = ax.get_xticklabels() + ax.get_yticklabels()
#[label.set_fontname('Times New Roman') for label in labels]
 
#设置横纵坐标的名称以及对应字体格式
font = {
'weight' : 'normal',
'size'   : 13,
}

#设置X轴Y轴名称  
plt.xlabel("Peptide",font)  
plt.savefig('Filename',bbox_inches = 'tight')          #保存图片命令一定要放在plt.show()前面
plt.show()


# In[ ]:




