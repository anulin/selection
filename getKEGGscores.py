import glob
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import random as ra
import matplotlib as mpl

def union(a,b):
    for i in b:
        a+=b[i]

x={}
y={}
z={}


metrics={}#словарь списков метрик  для  гена/ dictionary of snp metrics lists for genes
k=0

used=set()
c=2
print('''getKEGGscores scoresfile --a annotationfile --n networksfile'
--a (optional) --n (optional)
--h - help''')
exit()
if '--h' in sys.argv:
    print('''getKEGGscores scoresfile --a annotationfile --n networksfile
    scoresfile is produced by scores.dll 
    --a (optional) - custom file with snp to gene annotation
    --n (optional) - custom file with gene to pathway annotation
    --p (optional) - percentile. Default is 15%
    --h - help''')
    exit()
if '--a' in sys.argv and len(sys.argv)>=sys.argv.index('--a')+2:
    annotfile = sys.argv[sys.argv.index('--a')+1]
    c+=2
else:
    annotfile = "annotated_ExomeCeu.txt"
if '--n' in sys.argv and len(sys.argv)>=sys.argv.index('--n')+2:
    networkfile = sys.argv[sys.argv.index('--n')+1]
    c += 2
else:
    networkfile="KEGG to Genes Upd_ShortNams.txt"

if len(sys.argv)!=c or sys.argv[1]=='--a' or sys.argv[1]=='--n':
    raise Exception('Wrong arguments. Use --h to see how to run')#frqfile,  annot file, networkfile, flags
frqFile=sys.argv[1]


with open(frqFile) as file:#открываем файл с метриками и снипами и записываем в словарь х 
    for line in file:
        x[line.split()[0]]=line.split()[1]
# with open("scoresvcfExomeCeuMaf.txt") as file:
#     for line in file:
#         x[line.split()[0]+'_'+line.split()[2]]=line.split()[7]

#annotated_exomeCeuS annotated_ExomeS
genes=set()
with open(annotfile) as file:# открываем файл аннотации нужных сипов и генов и записываем их метрики, используя словарь х
    for l in file:
        if l.split()[1] in used or len(l.split())<3:#избегаем повторения снипов
            print(l)
            continue
        if( l.split()[1] in x and (l.split()[2] in genes or l.split()[2] in z.keys() and z[l.split()[2]] in genes)): #если название уже есть в словаре добавляем метрику по ключу в список
            if (l.split()[2] in z.keys()):# часть старого кода. Проверяет если название хорошее(hgnc)
                metrics[z[l.split()[2]]].append(float(x[l.split()[1]]))
            else:
                metrics[l.split()[2]].append(float(x[l.split()[1]]))
                z[l.split()[2]] = l.split()[2]
        elif  l.split()[1] in x : #иначе добавляем название(ключ) в словарь
            genes.add(l.split()[2])
            metrics[l.split()[2]]=[float(x[l.split()[1]])]
            z[l.split()[2]]=l.split()[2]
        used.add(l.split()[1])
Goodnames={}
pvalues={}
PathwayGenes={}
ListOfVals=[]
nnn=set()#genes present in networks
with open(networkfile) as file:#KEGG_database  #KEGG to Genes Upd
    file.readline()
    for line in file:
        if line.split()[1] not in PathwayGenes:
            Goodnames[line.split()[1]]=[line.split()[2]]
            PathwayGenes[line.split()[1]]=[line.split()[0]]
        else:
            PathwayGenes[line.split()[1]].append(line.split()[0])
        if line.split()[0] in metrics: #if gene is in network
            nnn.add( line.split()[0])
for i in nnn:
    ListOfVals+=metrics[i]
print(len(x),sum(len(i) for i in metrics))
print(len(ListOfVals))
exist='F'
for pathway in PathwayGenes:#открываем все xml файлы (файлы генных сетей) в папке KEGG_gml; walk through all the kegg names
    metric = []
    # metric=0
    names = set()
    for gene in PathwayGenes[pathway]:
        if gene in metrics:
            metric += metrics[gene]  # cумма листов - лист

            names.add(gene)
    counter = 0
    k = 0
    # print(filename,len(names), len(PathwayGenes[filename[9:-4]]))
    if (len(metric) == 0):
        pvalues[pathway] = '1.0 0\t0'
        continue
    metric =np.percentile(metric,15)
    snpCount = sum(len(metrics[gene]) for gene in names)

    if (metric == 0):
        pvalues[pathway] = '1.0 0 ' + str(snpCount)
    else:  # если метрика сети не 0 то генерируем случайные сети
        while (counter < 290 and k < 250000 and not (counter >= 9 and counter / k >= 0.05)):
            RandMetric = ra.sample(ListOfVals, snpCount)
            RandMetric = np.percentile(RandMetric, 15)
            if (metric <= RandMetric):  # если метрика больше то +1 к счетчику таких метрик. Далее counter делится на число генераций k с целью получения p-value
                counter += 1
            k += 1

        #print(filename, metric, counter / k, snpCount, len(names))

        pvalues[pathway] = str(counter / k) + '\t' + str((metric)) + '\t' + str(
            snpCount) + '\t' + str(len(names))  # записываем количество генов метрику и p-value
with open("pathway_scores.txt",'w') as file:#записываем результат в файл.
    file.write("Name\tID\tP-value\tMetric\tsnps\tgenes"+'\n')
    for i in pvalues:# i это название сети и ключ в словаре pvalues. по нему получаем строку с остальными значениями
        file.write(Goodnames[i]+'\t'+i +'\t'+str(pvalues[i]) +'\n')
