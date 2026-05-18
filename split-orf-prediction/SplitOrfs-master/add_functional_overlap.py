#!/usr/bin/python
#This script parses the out file of getLongestMatches.py to add overlap with functional annotations

import sys
import re

#read in fasta file
if len(sys.argv) < 3:
        print("usage python addFunctionalOverlap_py3.py UniqueProteinORFPairs.txt intersectResults.txt")
else :
    ORFPairs=open(sys.argv[1],'r')
    overlap=open(sys.argv[2],'r')
    annotations = {}  #hashmap that records ORF IDs (keys) and annotation IDs (values)

    for line  in overlap :
        line = line.rstrip()
        elems = line.split("\t")
        annotations[elems[4]] = elems[8]   #remember the last functional annotation for an ORF that exists in file

    #print new header for result file
    print("\t".join(["geneID", "targetTransID" ,"OrfTransID", "NumOrfs","OrfIDs",\
                        "OrfPos","OrfLengths","OrfSeqIdents","MinSeqIdent","MaxSeqIdent",\
                        "protAlignPos","protCoverage","ORF-DomainAnnot","NumORFAnnot","AnnotPercent"]))

    #go through unique protein pairs and add functional overlap of ORFs per entry
    next(ORFPairs)
    for line in ORFPairs :
        line = line.rstrip()
        countAnnot = 0
        listAnnot = []
        elems = line.split("\t")
        orfs = elems[4].split(",")
        #go through orfs and check whether annotations for them exist
        for j in orfs:
            if j in annotations :
                countAnnot = countAnnot +1
                listAnnot.append(annotations[j])
        if not listAnnot :
            listAnnot.append("NA")
        annotPercentage = round(float(countAnnot)/float(elems[3]),3)
        print("\t".join([line,",".join(listAnnot),str(countAnnot),str(annotPercentage)]))

#Example line in UniqueProteinORFPairs.txt         
#geneID  targetTransID   OrfTransID  NumOrfs OrfIDs  OrfPos  OrfLengths  OrfSeqIdents    MinSeqIdent MaxSeqIdent protAlignPos    protCoverage
#ENSMUSG00000000561  ENSMUST00000010278  ENSMUST00000000572  4   ORF-115,ORF-102,ORF-497,ORF-321 88396-88542,80071-80241,132564-132839,120917-121237 147,171,276,321 100.0,100.0,100.0,93.243    93.243  100.0   240-286,108-159,361-451,280-353 264

#bedtools output
#ENSMUST00000010278      108     159     ENSMUST00000000572      ORF-102 ENSMUST00000010278      123     153     PF00400
#ENSMUST00000152088      1       62      ENSMUST00000001171      ORF-50  ENSMUST00000152088      4       75      PF03388
