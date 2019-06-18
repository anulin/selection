import glob
import networkx as nx
import os
import matplotlib.pyplot as plt
import numpy as np
import random as ra
import matplotlib as mpl

x={}
y={}
z={}
"""z - хорошие названия генов, х-метрики снипов, y - названия генов у белков"""
"считает метрики сети как сумму метрик ребер gjltktyye. yf xbckj ht,th"

metrics={}

with open("results.txt") as file:
    for line in file:
            z[line.split()[1]]=line.split()[0]
with open("scores_rep.txt") as file:
    for line in file:
        x[line.split()[0]]=line.split()[1]
G=nx.Graph()
with open("genes.txt") as file:
    for l in file:
        if(len(l.split())==3 and l.split()[1] in x.keys() and (l.split()[2] in G.nodes or l.split()[2] in z.keys() and z[l.split()[2]] in G.nodes)):
            if (l.split()[2] in z.keys()):
                metrics[z[l.split()[2]]].append(float(x[l.split()[1]]))
                G.nodes[z[l.split()[2]]]['chi']+=float(x[l.split()[1]])
            else:
                G.nodes[l.split()[2]]['chi']+=float(x[l.split()[1]])
                metrics[l.split()[2]].append(float(x[l.split()[1]]))
                z[l.split()[2]] = l.split()[2]
        elif(len(l.split())==3  and l.split()[1] in x.keys() and float(x[l.split()[1]])>0):
            if(l.split()[2] in z.keys()):
                G.add_node(z[l.split()[2]], chi=float(x[l.split()[1]]))
                metrics[z[l.split()[2]]]=[float(x[l.split()[1]])]
            else:
                G.add_node(l.split()[2], chi=float(x[l.split()[1]]))
                metrics[l.split()[2]]=[float(x[l.split()[1]])]
                z[l.split()[2]]=l.split()[2]

with open("mart_export proteins.txt") as file:
    for line in file:
        if (line.split()[0] in G.nodes or line.split()[0] in z.keys() and z[line.split()[0]] in G.nodes):
            if(line.split()[0]in  z.keys()):
                y[line.split()[2]] = z[line.split()[0]]
            else:
                y[line.split()[2]] = line.split()[0]
print(G.nodes['CD36']['chi'],z['ENSG00000135218'])
print(G.nodes['PSCA']['chi'])
pvalues={}
for node in G.nodes:

    metrics[node].sort()
    G.nodes[node]['chi']=np.mean(metrics[node])
path = 'KEGG_gml'
numer=0
for filename in glob.glob(os.path.join(path,'*.xml')):
    metric=0

    using=set()
    names={}
    i=10
    relations=0
    edges=0
    with open(filename) as file:
        for line in file:
            if(not relations):
                if (line.split()[0]=="<entry"):
                    i=2
                    id = line.split('"')[1];
                if (line.split()[0]=="<relation"):
                    relations=1
                    i=1
                    a = names[line.split('"')[1]]
                    b = names[line.split('"')[3]]
                    if (a in G.nodes and b in G.nodes):
                        metric+= G.nodes[a]['chi'] * G.nodes[b]['chi']
                        using.add(a)
                        using.add(b)
                        edges+=1
                if i==4:
                    names[id]=line[24:].split(',')[0].split('"')[0]
                i += 1
            else:

                if line.split()[0]=='</pathway>':
                    break
                if line.split()[0]=="<relation":
                    a=names[line.split('"')[1]]
                    b=names[line.split('"')[3]]
                    if(a in G.nodes and b in G.nodes):
                        metric+= G.nodes[a]['chi'] * G.nodes[b]['chi']
                        if filename == 'hsa00230.xml' and a in using and b in using:
                            print(a,b)
                        using.add(a)
                        using.add(b)
                        edges+=1
        counter=0
        k=0
        if(metric==0):
            pvalues[filename]=(1,0)
        else:
            while(counter<80):
                RandVert = []
                for i in range(len(using)):
                    RandVert.append(ra.choice(list(G.nodes)))
                RandMetric=0
                for i in range(edges):
                    RandMetric+= G.nodes[ra.choice(RandVert)]['chi'] * G.nodes[ra.choice(RandVert)]['chi']
                if(metric<=RandMetric):
                    counter+=1
                k+=1
                if counter==0 and k>1000 or counter>0 and counter/k<10**-3:
                    print(edges)
            pvalues[filename]=str(counter/k)+' ' +str((metric/edges,len(using), edges, len(names)))

with open("метрики.txt",'w') as file:
    file.write("Name P-value Metric vertices edges genes"+'\n')
    for i in pvalues:
        file.write(i[9:-4] +' '+str(pvalues[i]) +'\n')