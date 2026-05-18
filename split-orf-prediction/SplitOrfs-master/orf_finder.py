# Orfinder script written by Marcel Schulz 2019
# currently supports only Python 2.7
# Orffinder finds all possible ORFs of a minimal length in a given transcript and outputs
# to stream a fasta file with the protein sequences of these orfs
# each such orf has the header formatted like this: 
# >FastaId:ORF-X:Start:Stop
# ORF-X denotes the x-th orf found for sequence with name FastaId and the transcript relative start and stop positions

# usage OrfFinder_py3.py transcripts.fa

#!/usr/bin/python

import sys
import re
from Bio import SeqIO


minOrfLength=50  #minimum number of amino acids per ORF

#original codons functions by natasha.sernova obtained from Biostars:
#https://www.biostars.org/p/229060/
#code has been modified

def codons(seq,id,countOrfs):

        stops = ['TAA','TGA','TAG']    
        lst1 = [] #List for the start codons
        lst2 = [] #List for the stop codons
        start = 0 #The start position of the sequence.
        counter = 0 #Counter for 3 optional orfs.
        #initializes the lists for 3 optional orfs.
        for i in range (3):
            lst1.append([])
            lst2.append([])
        #Add to the lists the positions of the start and stop codons.
        while (seq and counter < 3):

            for i in range(start,len(seq),3):
                codon = seq[i:i+3] #The codon is 3 nucleotides.
                #print codon+ '\t'
                if(codon == 'ATG'): #If the codon is  a start codon.
                    lst1[start].append(i) #The position of the start codon.

                if(codon in stops): #if the codon is a stop codon.
                    lst2[start].append(i) #The position of the stop codon.


            start += 1 #promotes the starting position.
            counter += 1 #promotes the counter

        #for each reading frame go through the start site and extract and output possible proteins
        for frame in range (3):
                if len(lst1[frame])>0 and len(lst2[frame])>0:  #at least one start and stop codon per frame must exist
                
                    currentStart=0 #the current start position of a start codon in the frame
                    currentStop=0  #the current start position of a stop codon in the frame
                    while  currentStart < len(lst1[frame]) and currentStop < len(lst2[frame]): #codons are available
                        if lst1[frame][currentStart] < lst2[frame][currentStop] :  #found valid ORF
                                #translate orf sequence
                                prot=translate(seq[lst1[frame][currentStart]:lst2[frame][currentStop]])
                                if len(prot) >= minOrfLength :
                                        ##### added Frame to ID #####
                                        #print(''.join(['>',id,':ORF-',str(countOrfs),':',str(lst1[frame][currentStart]),':',str(lst2[frame][currentStop]+3)]))
                                        print(''.join(['>',id,':ORF-',str(countOrfs),':',str(lst1[frame][currentStart]),':',str(lst2[frame][currentStop]+3)]))
                                        #depends on the coordinate system to use, if one based then this with 3 would be correct, but we start counting the 
                                        #orfs from 0 and not from 1, so this should be +2 to stay in the 0-based coordinate system
                                        print(prot)
                                        countOrfs=countOrfs+1
                                #remove all other start codons that are nested between the current start and stop
                                while currentStart < len(lst1[frame]) and lst1[frame][currentStart] < lst2[frame][currentStop]:
                                        currentStart=currentStart+1
                                currentStop=currentStop+1
                                
                        elif lst1[frame][currentStart] > lst2[frame][currentStop]:
                                currentStop=currentStop+1 
        return(countOrfs)        

#translation code taken from https://www.geeksforgeeks.org/dna-protein-python-3/
def translate(seq): 
       
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
    } 
    protein ='' 
    if len(seq)%3 == 0 and not re.search('N',seq): #if seq contains N discard it
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            protein+= table[codon]
    #else:
    #    print 'not correct length'       
    return protein 


#read in fasta file
if len(sys.argv) < 2:
        print('usage OrfFinder_py3.py transcripts.fa')
else :

        transcripts=SeqIO.parse(sys.argv[1],'fasta')

        countOrfs=1
        for record in transcripts:
                countOrfs=codons(str(record.seq), record.id ,countOrfs)
                
        

