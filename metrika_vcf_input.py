import glob
import networkx as nx
import os
import math
import matplotlib.pyplot as plt
import numpy as np
import random as ra
import matplotlib as mpl

x={}
y={}
z={}
"""z - использовалось в старой программе, х-метрики снипов, y - названия генов у белков"""
"считает метрики сети как сумму метрик генов в сети"

metrics={}#словарь списков метрик используемый для нахождения медианной метрики гена
k=0
#код от старой программы
with open("results.txt") as file:
    for line in file:
            z[line.split()[1]]=line.split()[0]
#конец

with open("scoresvcf.txt") as file:#открываем файл с метриками и снипами и записываем в словарь х
    for line in file:
        x[line.split()[0]]=line.split()[1]
G=nx.Graph() #объявляем массив генов и их метрик
with open("genespos.txt") as file:# открываем файл аннотации нужных сипов и генов и записываем их метрики, используя словарь х
    for l in file:
        if( (l.split()[2] in G.nodes or l.split()[2] in z.keys() and z[l.split()[2]] in G.nodes)): #если название уже есть в словаре добавляем метрику по ключу в список
            if (l.split()[2] in z.keys()):# часть старого кода. Проверяет если название хорошее(hgnc)
                metrics[z[l.split()[2]]].append(float(x[l.split()[1]]))
            else:
                metrics[l.split()[2]].append(float(x[l.split()[1]]))
                z[l.split()[2]] = l.split()[2]
        elif(  float(x[l.split()[1]])>0): #иначе добавляем название(ключ) в словарь
            if(l.split()[2] in z.keys()):
                G.add_node(z[l.split()[2]], chi=float(x[l.split()[1]]))
                metrics[z[l.split()[2]]]=[float(x[l.split()[1]])]
            else:
                G.add_node(l.split()[2], chi=float(x[l.split()[1]]))
                metrics[l.split()[2]]=[float(x[l.split()[1]])]
                z[l.split()[2]]=l.split()[2]

"""with open("mart_export proteins.txt") as file:
    for line in file:
        if (line.split()[0] in G.nodes or line.split()[0] in z.keys() and z[line.split()[0]] in G.nodes):
            if(line.split()[0]in  z.keys()):
                y[line.split()[2]] = z[line.split()[0]]
            else:
                y[line.split()[2]] = line.split()[0]
print(G.nodes['CD36']['chi'],z['ENSG00000135218'])"""
"""print(G.nodes['PSCA']['chi'])"""
pvalues={}
for node in G.nodes:#находим медианы записываем их в массив G генов и метрик

    metrics[node].sort()
    G.nodes[node]['chi']=np.median(metrics[node])
    if(math.isnan(G.nodes[node]['chi'])):
        print(metrics[node])
path = 'KEGG_gml'
numer=0
for filename in glob.glob(os.path.join(path,'*.xml')):#открываем все xml файлы (файлы генных сетей) в папке KEGG_gml
    metric=0


    names={}
    i=10
    relations=0
    edges=0
    with open(filename) as file:
        for line in file:
            if(not relations):
                if (line.split()[0]=="<entry"):#после "<entry" содержится индекс (переменная id), а в четвертой строке содержится название гена
                    i=2
                    id = line.split('"')[1];#индекс не используется здесь. Он использовался в программе которая вычисляла метрики ребер
                if (line.split()[0]=="<relation"):
                    relations=1
                    i=1
                if i==4:#проверка что строка четвертая
                    if(line[24:].split(',')[0].split('"')[0] in G.nodes): #если такой ген есть в списке метрик и генов (G) то запоминаем его имя, прибавляем к метрике сети его метрику
                        names[id] = line[24:].split(',')[0].split('"')[0]
                        metric+=G.nodes[names[id]]['chi']
                i += 1
            else:

                if line.split()[0]=='</pathway>':# дальше нужной информации нет
                    break
                if line.split()[0]=="<relation":
                    break
        counter=0
        k=0
        if(metric==0):
            pvalues[filename]='1.0 0 '+str(len(names))
        else: #если метрика сети не 0 то генерируем случайные сети
            while(counter<10000 and k < 100000):

                RandMetric=0
                for i in range(len(names)):#берем столько случайных генов, сколько удалось опознать в сети

                    RandMetric+=G.nodes[ra.choice(list(G.nodes))]['chi']
                if(metric<=RandMetric):#если метрика больше то +1 к счетчику таких метрик. Далее counter делится на число генераций k с целью получения p-value
                    counter+=1
                k+=1

            print(filename)
            pvalues[filename]=str(counter/k)+' ' +str((metric/len(names)))+' '+ str(len(names))#записываем количество генов метрику и p-value

with open("метрики.txt",'w') as file:#записываем результат в файл.
    file.write("Name P-value Metric vertices genes"+'\n')
    for i in pvalues:# i это название сети и ключ в словаре pvalues. по нему получаем строку с остальными значениями
        file.write(i[9:-4] +' '+str(pvalues[i]) +'\n')