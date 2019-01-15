#!/usr/bin/python
#This script parses a blastp output file and checks for NMD transcripts that have more than one ORF matching to the
#peptide sequence of one
#the script expects that the blastpoutput file was sorted on the database 
#sequence column (2nd column), e.g. via srt -k2 align.out

import sys
import re

from itertools import groupby

minOrfNum=2  #minimum number ORFs that need to align to a peptide sequence
identityCutoff = 90  #minimum percent identity of the protein alignment to be considered as a valid ORF-peptide match
colon=":"


def checkAlignments(Alignments,gene, target) :
    if Alignments :
                    for match in Alignments:
                        entries = Alignments[match]
                        if len(entries) > 1: #found more than one ORF 
                            orfIDs =[]
                            orfPos =[]
                            for i in range(0,len(entries)) :    
                                orfIDs.append(entries[i].split(":")[0])
                                orfPos.append("-".join([entries[i].split(":")[1],entries[i].split(":")[2]]))

                            print "\t".join([gene,target,match,str(len(entries)),",".join(orfIDs),",".join(orfPos)])

#read in fasta file
if len(sys.argv) < 2:
        print "usage DetectValidNMDOrfMatches.py BlastpAlign.out"
else :

        file=open(sys.argv[1],'r')

        lastElem="dummy" #variable that remembers the last transcript in use (all alignments are sorted by traget sequence)
        Alignments = {}  #hash map that stores all valid alignments for one target sequence
        #print header for result file
        print "\t".join(["geneID", "targetTransID" ,"OrfTransID", "NumOrfs","OrfIDs","OrfPos"])

        for line in file :
            elems = line.split("\t")
            orf=re.split("[|:]+",elems[0])
            target=elems[1].split("|")
            if lastElem != elems[1] :
                #found a new target transcript
                #go through all transcripts that have matches to the target
                checkAlignments(Alignments,target[0],target[1])
                lastElem = elems[1]
                Alignments={}

            if float(elems[2]) >= identityCutoff and (orf[0] == target[0]):
                if(orf[1] in Alignments) :  #transcript has at least one matching orf
                    Alignments[orf[1]].append(colon.join(orf[2:5]))
                else :
                    Alignments[orf[1]] = [colon.join(orf[2:5])]
            
        checkAlignments(Alignments,target[0],target[1])



#ENSMUSG00000031748|ENSMUST00000127900:ORF-189:67160:67360  ENSMUSG00000000001|ENSMUST00000000001   