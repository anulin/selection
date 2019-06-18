import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import matplotlib as mpl
x={}

text=[]
values=[]
i=0
k=0
with open("KEGG_database.txt") as file:
    for line in file:
            if(len(line.split())>1):
                x[line.split()[1]]=line.split()[2]
with open("метрики.txt") as file:
    for line in file:
        if k==1:
            if float(line.split()[1])<=0.01:
                text.append(x[line.split()[0]] + ' ' +' '.join(line.split()[1:]) + ' 1e')
            else:
                text.append(x[line.split()[0]] + ' ' + ' '.join(line.split()[1:]) + ' 0e')
            values.append(float(line.split()[2]) )
        else:
            k=1

values2=st.zscore(values)
with open('метрики3.txt','w') as file:
    for i in range(len(text)):
        file.write(text[i][:-1]+' '+str(values2[i])+'\n')

