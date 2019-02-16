# addExternalIDs2 script written by Marcel Schulz 2019
# currently supports only Python 2.7
# This script can be used to add Ensembl gene IDs to refseq sequences

# usage addExternalIDs2.py IDtable.txt RefseqProteins.fa

#!/usr/bin/python

import sys
import re

from itertools import groupby


#read in fasta file
if len(sys.argv) < 2:
        print "usage addExternalIDs2.py IDtable.txt RefseqProteins.fa"
else :

        file=open(sys.argv[1],'r')
        ProteinMap = {}
        #go through protein fasta file to add IDs
        for line in file :
            #line=line.strip()
            elems = line.split("\t")
            if len(elems) == 2:
                elems[1]=elems[1].strip()
                #print  elems[1],elems[0]
                ProteinMap[elems[1]] = elems[0]

        
        file=open(sys.argv[2],'r')
        output = 0
        for line in file :
        	if line[0] == ">" :
				elems = re.split("[>.]",line)
				
				if elems[1] in ProteinMap:
					print "".join([">",ProteinMap[elems[1]],"|",elems[1]])
					output = 1
				else :
					output = 0
        	elif output == 1:
        		print line.strip()
                
        