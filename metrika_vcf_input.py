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
print(np.percentile((1,2,5,8,9),25))
metrics={}#словарь списков метрик используемый для нахождения медианной метрики гена
k=0
#код от старой программы
# with open("results.txt") as file:
#     for line in file:
#             z[line.split()[1]]=line.split()[0]
used=set()
with open("scoresvcfExomeCeuMaf.txt") as file:#открываем файл с метриками и снипами и записываем в словарь х
    for line in file:
        x[line.split()[0]]=line.split()[1]
G=nx.Graph() #объявляем массив генов и их метрик
with open("annotated_Exome.txt") as file:# открываем файл аннотации нужных сипов и генов и записываем их метрики, используя словарь х
    for l in file:
        if l.split()[1] in used or len(l.split())<3:#избегаем повторения снипов
            continue

        if( l.split()[1] in x and (l.split()[2] in G.nodes or l.split()[2] in z.keys() and z[l.split()[2]] in G.nodes)): #если название уже есть в словаре добавляем метрику по ключу в список
            if (l.split()[2] in z.keys()):# часть старого кода. Проверяет если название хорошее(hgnc)
                metrics[z[l.split()[2]]].append(float(x[l.split()[1]]))
            else:
                metrics[l.split()[2]].append(float(x[l.split()[1]]))
                z[l.split()[2]] = l.split()[2]
        elif  l.split()[1] in x : #иначе добавляем название(ключ) в словарь
            G.add_node(l.split()[2], chi=float(x[l.split()[1]]))
            metrics[l.split()[2]]=[float(x[l.split()[1]])]
            z[l.split()[2]]=l.split()[2]
        used.add(l.split()[1])

pvalues={}

for node in metrics:#находим медианы записываем их в массив G генов и метрик
    if (len(metrics[node]) < 1):
        G.remove_node(node)
    else:
        G.nodes[node]['chi'] = np.percentile(metrics[node],50)
    # elif (max(metrics[node]) == np.median(metrics[node])):
    #     G.nodes[node]['chi'] = 1
    # else:
    #     metrics[node].sort()
    #     G.nodes[node]['chi']=(max(metrics[node])-np.median(metrics[node]))/(max(metrics[node])-metrics[node][len(metrics[node])-2])
    #     if(math.isnan(G.nodes[node]['chi'])):
    #         print(metrics[node])
PathwayGenes={}
with open("KEGG_database.txt") as file:
    for line in file:
        if line.split()[1] not in PathwayGenes:
            PathwayGenes[line.split()[1]]=[line.split()[0]]
        else:
            PathwayGenes[line.split()[1]].append(line.split()[0])

path = 'KEGG_gml'
numer=0
ListOfVals=[G.nodes[node]['chi']for node in G]
# ListOfVals=[float(i)for i in x.values()]
for filename in glob.glob(os.path.join(path,'*.xml')):#открываем все xml файлы (файлы генных сетей) в папке KEGG_gml
    metric = 0
    names = set()
    i = 10
    relations = 0
    edges = 0
    with open(filename) as file:
        for gene in PathwayGenes[filename[9:-4]]:
            if gene in G.nodes:
                metric += G.nodes[gene]['chi']
                names.add(gene)
        counter=0
        k=0
        if (metric == 0):
            pvalues[filename] = '1.0 0 ' + str(len(names))
        else:  # если метрика сети не 0 то генерируем случайные сети
            while (counter < 290 and k < 80000 and not (counter >= 10 and counter / k >= 0.1)):
                RandMetric = 0
                for i in range(len(names)):  # берем столько случайных генов, сколько удалось опознать в сети
                    RandMetric += ra.choice(ListOfVals)
                if (metric <= RandMetric):  # если метрика больше то +1 к счетчику таких метрик. Далее counter делится на число генераций k с целью получения p-value
                    counter += 1
                k += 1

            print(filename,metric/len(names),counter/k, len(names))
            pvalues[filename]=str(counter/k)+'\t' +str((metric)/len(names))+'\t'+ str(len(names))#записываем количество генов метрику и p-value
with open("метрики.txt",'w') as file:#записываем результат в файл.
    file.write("Name P-value Metric vertices genes"+'\n')
    for i in pvalues:# i это название сети и ключ в словаре pvalues. по нему получаем строку с остальными значениями
        file.write(i[9:-4] +'\t'+str(pvalues[i]) +'\n')