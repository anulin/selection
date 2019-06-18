count=0
dict={'A' :'T', 'T':'A','C':'G','G':'C'}
text=[]

with open('dataSkipogh.txt') as file:
    for lin in file:
        abc=lin.split()
        if(abc[2]=='0' and abc[3]=='0'):
            if(dict[abc[6]]==abc[8]):
                line=abc[:2]+[abc[4] ,abc[5], '0 0']+abc[6:9]
            elif dict[abc[7]]==abc[8]:
                line = abc[:2] + [abc[5] , abc[4],'0 0']+ abc[6:9]
            text.append(' '.join(line))
        #elif(abc[2]=='0'):
         #   line = abc[:2]+[abc[3],abc[2],abc[5],abc[4]]+abc[6:9]
          #  text.append(' '.join(line))
        elif (not (abc[7] == dict[abc[6]] and abs(int(abc[2]) - int(abc[3])) < 687)):
            text.append(lin[:-1])
with open('RepairedSkipogh2.txt','w') as file:
    for lin in text:
        file.write(lin+'\n')
print(len(text))