#!/usr/bin/python
#create bed file for ORF matches of a transcript


#conversion line for Ensembl putput
#awk 'BEGIN{OFS="\t"}{print $2,$4,$5,$3}' EnsemblMouse95PFAMannotation.txt > PfamEnsemblMouse.bed

import sys
import re

#read in fasta file
if len(sys.argv) < 2:
        print "usage makeBed.py UniqueProteinORFPairs.txt"
else :

        file=open(sys.argv[1],'r')
        #skip header
        next(file)
        for line in file :
            elems = line.split("\t")
            matches=re.split(",",elems[10])
            orfIDs=re.split(",",elems[4])

            for i in range(0,len(matches)) :
                positions=re.split("-",matches[i])
                #joint=":".join([elems[2],orfIDs[i]])
                print "\t".join([elems[1],positions[0],positions[1],elems[2],orfIDs[i]])

#geneID targetTransID   OrfTransID  NumOrfs OrfIDs  OrfPos  OrfLengths  OrfSeqIdents    MinSeqIdent MaxSeqIdent protAlignPosprotCoverage
#ENSMUSG00000000561  ENSMUST00000010278  ENSMUST00000000572  4   ORF-115,ORF-102,ORF-497,ORF-321 88396-88542,80071-80241,132564-132839,120917-121237 147,171,276,321 100.0,100.0,100.0,93.243    93.243  100.0   240-286,108-159,361-451,280-353 264