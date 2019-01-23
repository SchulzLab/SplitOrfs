
#!/usr/bin/python
#This script parses the out file of DetectValidNMDOrfMatches.py to select the one best transcript-protein pair
#for further analysis. For each transcript it uses the some of ORF-protein lengths to determine, which protein
#would show the largest coverage. For this the file needs to be sorted according to column 3 (see runNMDOrf.sh):
# THis file is named ValidProteinORFPairs_sortCol3.txt in the NMDOrf workflow

import sys


#read in fasta file
if len(sys.argv) < 2:
        print "usage getLongestOrfMatches.py ValidProteinORFPairs_sortCol3.txt"
else :

        file=open(sys.argv[1],'r')
        lastElem="dummy" #variable that remembers the last transcript in use (all alignments are sorted by source sequence)
        bestProteinMatch = "" #rembers the line where the current 
        Matchlength = 0
        #reprint header for result file
        print "\t".join(["geneID", "targetTransID" ,"OrfTransID", "NumOrfs","OrfIDs","OrfPos","OrfLengths","OrfSeqIdents","MinSeqIdent","MaxSeqIdent","protAlignPos","protCoverage"])

        #skip header
        next(file)
        for line in file :

            elems = line.split("\t")
            if lastElem != elems[2] :
                #found a new source transcript
                #output best proteins that has longest matches from ORFs of the source
                if lastElem != "dummy":
                	print bestProteinMatch
                lastElem = elems[2]
            	bestProteinMatch = ""
            	Matchlength = 0
            if int(elems[11]) > Matchlength 	:
                bestProteinMatch = line.rstrip()
                Matchlength = int(elems[11])

        print bestProteinMatch  #print last source transcript in file
            
