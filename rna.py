# 20 angstrom max
import math
pdbIn = open ("4kzz.pdb", "r")
pdbOut = open ("4kzzC3.pdb", "w+")

myline = pdbIn.readline()
chains = dict()
seq= []
while myline:
    if 'ATOM' in myline[0:6]:
        if "C3'" in myline[12:16]:
            if not (myline[21:22] in chains) :            
                seq = []
                seq.append(myline[19])
                chains[myline[21:22]] = seq
            if myline[21:22] in chains : 
                seq.append(myline[19])
                chains[myline[21:22]] = seq   
            pdbOut.write(myline)
    myline = pdbIn.readline()

for key in chains:
    chains[key] = ''.join(chains[key])

pdbIn.close()
pdbOut.close()

pdbIn1 = open ("4kzzC3.pdb", "r")

mylinef = pdbIn1.readline()

i = 0
duos = dict()

fileList = []
while mylinef :
    fileList.append(mylinef)
    mylinef = pdbIn1.readline()


for key in chains:
    i = 0
    while i < len(fileList):
        
        j = i + 1  
        # print (j)      
        if key in fileList[i][21:22]:
            xi = float(fileList[i][30:38])
            yi = float(fileList[i][38:46])
            zi = float(fileList[i][46:54])
            while j < len(fileList):
                
                if key in fileList[j][21:22]:                    
                    xj = float(fileList[j][30:38])
                    yj = float(fileList[j][38:46])
                    zj = float(fileList[j][46:54])
                    dist = math.sqrt(((xj-xi)**2)+((yj-yi)**2)+((zj-zi)**2))
                    
                    duo = fileList[i][19] + fileList[j][19]
                    revDuo = fileList[j][19] + fileList[i][19]
                    if not duo in duos: 
                        duos[duo]=[]                       
                    if not revDuo in duos:
                        duos[revDuo] = []
                            
                    if dist < 20:
                        if duo in duos :                            
                                duos[duo].append(dist)
                        if revDuo in duos :                            
                                duos[revDuo].append(dist)
                j+=1
        i += 1
pdbIn1.close()
# print (duos['AA'])                  

ranges = dict()
for i in range (20):
    
    for key in duos:
        count = 0
        if not key in ranges:
            ranges[key] = []
        for val in duos[key]:
            if val > i and val < i+1: 
                count += 1    
        ranges[key].append(count)

# print (ranges)

tot = []
for key in ranges :
    tot.append(sum(ranges[key]))
# print (tot)

fobs = []
for i in range (20):
    j = 0
    rangeI = []
    for key in ranges :
        f = ranges[key][i]/tot[j]
        j += 1
        rangeI.append(f)
    fobs.append(rangeI)
# print(fobs)

fref = []

for i in range (20):
    count = 0
    for key in ranges:
        count += ranges[key][i]
    fref.append(count/sum(tot))

# print (fref)

e = []

for i in range (len(fobs)):
    eline = []
    for j in range (len(fobs[i])):
        # print (fobs[i][j])
        if fobs[i][j]==0.0 or fref[j] == 0.0:
            calc = 10
        else:
            calc = -1 * math.log(fobs[i][j]/fref[j])
        eline.append(calc)
    e.append(eline)
# print(e)
    
pairSeq = []

for key in ranges :
    pairSeq.append(key)

energies = {}
for list in e:
    for i in range (len(list)):
        if not pairSeq[i] in energies :
            energies[pairSeq[i]]=[]
        energies[pairSeq[i]].append(list[i])

print(energies)



from matplotlib import pyplot as plt
from matplotlib.pyplot import figure

figure(figsize=(20,15), dpi=80)
i=1

for pair, scores in energies.items():
    plt.subplot(4,5,i)
    plt.title(f"{pair}")
    plt.axhline(y=0)
    plt.plot(range(20), energies[pair]) 
    i+=1
    
plt.show()