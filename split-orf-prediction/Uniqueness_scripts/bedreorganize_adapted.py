import sys
uniquefile = open(sys.argv[1],'r')
uniquelist=[]
for line in uniquefile:
    elements = line.strip()
    elements = elements.split("\t")
    ids = elements[0].split(":")
    temp=[ids[0],elements[1],elements[2],ids[1],ids[2],ids[3]]
    uniquelist.append(temp)
with open(sys.argv[2],'w') as f:
    for i in uniquelist:
        #give the unique regions in terms of transcript not of the orf, needed for intersection with riboseq data
        f.write(i[0] + "\t" + str(int(i[1]) + int(i[4])) + "\t" + str(int(i[2]) + int(i[4])) + "\t" + i[3] +  ":" + i[4] + ":" + i[5] + "\n")