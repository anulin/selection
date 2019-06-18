import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import random as ra

x={}
y={}
z={}
"""z - хорошие названия генов, х-метрики снипов, y - названия генов у белков"""
denominator={}
k=0
stat=[]
with open("results.txt") as file:
    for line in file:
            z[line.split()[1]]=line.split()[0]
d=0
count=0
with open("scores_rep.txt") as file:
    for line in file:
        x[line.split()[0]]=line.split()[1]
        if float(line.split()[1])>d :
            d=float(line.split()[1])
        if float(line.split()[1]) > 70:
            print(line)

#print(d,count)
G=nx.Graph()
e=0
with open("genes.txt") as file:
    for l in file:
        if(e<6 and l.split()[0]=='g'):
            print(l)
            e+=1
        if(len(l.split())==3 and l.split()[1] in x.keys() and (l.split()[2] in G.nodes or l.split()[2] in z.keys() and z[l.split()[2]] in G.nodes)):
            if (l.split()[2] in z.keys()):
                denominator[z[l.split()[2]]].append(float(x[l.split()[1]]))
                G.nodes[z[l.split()[2]]]['chi']+=float(x[l.split()[1]])
            else:
                G.nodes[l.split()[2]]['chi']+=float(x[l.split()[1]])
                denominator[l.split()[2]].append(float(x[l.split()[1]]))
                z[l.split()[2]] = l.split()[2]
        elif(len(l.split())==3  and l.split()[1] in x.keys() and float(x[l.split()[1]])>0):
            if(l.split()[2] in z.keys()):
                G.add_node(z[l.split()[2]], chi=float(x[l.split()[1]]))
                denominator[z[l.split()[2]]]=[float(x[l.split()[1]])]
            else:
                G.add_node(l.split()[2], chi=float(x[l.split()[1]]))
                denominator[l.split()[2]]=[float(x[l.split()[1]])]
                z[l.split()[2]]=l.split()[2]
with open("mart_export proteins.txt") as file:
    for line in file:
        if (line.split()[0] in G.nodes or line.split()[0] in z.keys() and z[line.split()[0]] in G.nodes):
            if(line.split()[0]in  z.keys()):
                y[line.split()[2]] = z[line.split()[0]]
            else:
                y[line.split()[2]] = line.split()[0]

for node in G.nodes:
    denominator[node].sort()
    G.nodes[node]['chi']=np.median(denominator[node])
print(G.nodes['CD36']['chi'],z['ENSG00000135218'])
print(G.nodes['PSCA']['chi'])
with open("9606.protein.links.v11.0.txt") as file:
    for line in file:
        if(line.split()[0][5:] in y.keys() and line.split()[1][5:] in y.keys()):
            stat.append(G.nodes[y[line.split()[0][5:]]]['chi'] * G.nodes[y[line.split()[1][5:]]]['chi'])
    stat.sort()
    a=round((len(stat))*0.9999995)-1
    topedg=stat[a:]
    scr=[]
    count=0
    nodes=[b[1] for b in G.nodes('chi')]
    nodes.sort()
    for product in topedg:
        for node in nodes:
            i = (len(nodes) + 1) // 2
            j = i
            if node != 0:
                while i <= len(nodes) and node * nodes[i - 1] != product:
                    if (j != 1):
                        j = (j) // 2
                    if (node * nodes[i - 1] > product and node * nodes[i - 2] < product):
                        break
                    if (node * nodes[i - 1] > product):
                        i -= j
                    else:
                        i += j
                count += len(nodes) - i + 1

        scr.append(count / len(nodes) ** 2)
        count = 0

    print(scr)
    a=stat[a]
    print(a, len(stat), len(G.nodes))

with open("9606.protein.links.v11.0.txt") as file:
    for line in file:
        if (line.split()[0][5:] in y.keys() and line.split()[1][5:] in y.keys() and G.nodes[y[line.split()[0][5:]]][
            'chi'] * G.nodes[y[line.split()[1][5:]]]['chi'] >= a):
            G.add_edge(y[line.split()[0][5:]], y[line.split()[1][5:]], weight=int(line.split()[2]),
                       edge_color=G.nodes[y[line.split()[0][5:]]]['chi'] * G.nodes[y[line.split()[1][5:]]]['chi'])
A = list(G.nodes)
for i in A:
    if (G.degree[i] < 1):
        G.remove_node(i)

colors = list(G[i][j]['edge_color'] for (i, j) in G.edges())
nx.draw(G, with_labels=True, edge_cmap=plt.cm.Blues, edge_color=colors)

plt.show()