N1={}
alleles={}
N0={}
scores={}
counter=0
#получает скоры и чинит их
dict={'A' :'T', 'T':'A','C':'G','G':'C'}
with open("data.txt") as file:
    for line in file:
        l=line.split()
        N0[l[0]]=[int(l[1]),int(l[2]),int(l[3]),0,0]
with open("RepairedSkipogh2.txt") as file:
    for line in file:
        l=line.split()
        n0=N0[l[0]][0]*2
        max=0
        maxi=1
        n1=int(l[1])
        for i in range(4):

            p0=N0[l[0]][i+1]/n0
            p1=int(l[i+2])/n1
            if(l[7]==dict[l[6]] and abs(p1-p0)>0.17 and (abs(p0-int(l[i + 3]) / n1)<abs(p1-p0) or p0==0)):
                if(p0==0):
                    p0 = N0[l[0]][i + 2] / n0
                else:
                    p1 = int(l[i + 3]) / n1
                if p0==p1:
                    break
                if (p1 < 0.65 and p1 > 0.35):
                    counter += 1
                max = (p0 - p1) * (p0 - p1) / (
                            p0 * (1 - p0) / n0 + p0 * p0 - 2 * p0 * p1 + p1 * (1 - p1) / n1 + p1 * p1 - (p1 - p0) ** 2)
                maxi = i + 1
                break
            if(p1!=p0 and p1>1/n1 and p0>1/n0):
                a=(p0 - p1) * (p0 - p1) / ( p0 * (1 - p0)/n0 + p0 * p0 - 2 * p0 * p1 +  p1 * (1 - p1)/n1 +  p1  * p1-(p1-p0)**2)
                if(max<a):
                    max=a
                    maxi=i+1

        scores[l[0]]=[max,maxi]+l[5:]
    print(counter)
with open("scores_rep.txt",'w') as file:
    for key in scores.keys():
        file.write(key+' '+str(scores[key])+"\n")