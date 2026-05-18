import sys

file = open(sys.argv[1],'r')
validlist = {}
for line in file:
    orf_info = line.strip().split("\t")
    IDs = orf_info[4].split(",")
    for i in IDs:
        validlist[i] = "Valid"
   

file2 = open(sys.argv[2],'r')
with open(sys.argv[3],'w') as f:
    for line in file2:
        orf_info = line.split("\t")
        ID = orf_info[0].split(':')[1]
        if ID in validlist.keys():
            f.write(line)
