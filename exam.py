pdbIn = open ("3iwm.pdb", "r")
# pdbOut = open ("3iwm.pdb", "w+")

myline = pdbIn.readline()

# The list that will contain our file
pdbList = []

# Putting in the list only the lines that refers to ATOMs
while myline : 
    if 'ATOM' in myline[0:6]:
        pdbList.append(myline)
    myline = pdbIn.readline()

# The dictionnary that will contain each chains line
chains = {}

# Fill the dictionnary. Each key is a chain, 
# and each value is the list of corresponding lines
for line in pdbList :
    if not (line[21:22] in chains) : 
        chains[line[21:22]] = []
        chains[line[21:22]].append(line)
    else : 
        chains[line[21:22]].append(line)

# Create a file for each chain (key of the dict)
# and fill it with the corresponding lines
for key in chains:
    pdbOut = open ('3iwm'+key+'.pdb', 'w+')
    for line in chains[key]:
        pdbOut.write(line)
    pdbOut.close()




