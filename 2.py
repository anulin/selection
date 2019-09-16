import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import matplotlib as mpl
x={}

text=[]
values=[]
i=0
k=0
with open("KEGG_database.txt") as file:#открываем файл с хорошими и плохими названиями сетей и записываем их в словарь
    for line in file:
            if(len(line.split())>1):
                x[line.split()[1]]=line.split()[2]
with open("метрики.txt") as file: #открываем файл с метриками сетей, запоминаем данные
    for line in file:
        if k==1:
            if float(line.split()[1])<=0.01:
                text.append(x[line.split()[0]] + '\t' +'\t'.join(line.split()[1:]) + '\t1e')
            else:
                text.append(x[line.split()[0]] + '\t' + '\t'.join(line.split()[1:]) + '\t0e')
            values.append(float(line.split()[2]) )
        else:
            k=1

values2=st.zscore(values)# считаем z-score
with open('метрики3.txt','w') as file:# записываем хорошие названия сетей,  их метрики и данные а также z-score
    for i in range(len(text)):
        file.write(text[i][:-1]+'\t'+str(values2[i])+'\n')

