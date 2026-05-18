#!/usr/bin/python
# This script parses a blastp output file and checks for transcripts that have more than one ORF matching to the
# peptide sequence of one transcript of the same gene
# the script expects that the blastpoutput file was sorted on the database
# sequence column (2nd column), e.g. via srt -k2 align.out

import sys
import re

from itertools import groupby

minOrfNum = 2  # minimum number ORFs that need to align to a peptide sequence
# minimum percent identity of the protein alignment to be considered as a valid ORF-peptide match
identityCutoff = 95
minLength = 50  # minimum number of amino acids for an ORF to be considered
# rate of positions to be covered by the alignment of an ORF (to remove spurious local protein alignments)
minAlignmentRate = 0.5
colon = ":"
# require a minimum length of the aligned sequence
minAlignLength = 15


def checkAlignments(Alignments, gene, target):
    if Alignments:
        for match in Alignments:
            entries = Alignments[match]
            if len(entries) >= minOrfNum:  # found more than one ORF
                # compute max and min sequence identity values
                MinSeqIdent = float(100)
                MaxSeqIdent = float(0)

                orfIDs = []
                orfPos = []
                orfSeqIdents = []
                orfLengths = []
                pAlignPos = []
                totalAlignLength = 0
                for i in range(0, len(entries)):
                    # print("**",entries)
                    orfIDs.append(entries[i].split(":")[0])
                    orfPos.append(
                        "-".join([entries[i].split(":")[1], entries[i].split(":")[2]]))
                    currentSeqIdent = float(entries[i].split(":")[3])
                    if (currentSeqIdent < float(MinSeqIdent)):
                        MinSeqIdent = currentSeqIdent
                    if (currentSeqIdent > MaxSeqIdent):
                        MaxSeqIdent = float(currentSeqIdent)
                    orfSeqIdents.append(str(currentSeqIdent))
                    orfLengths.append(
                        str(int(entries[i].split(":")[2])-int(entries[i].split(":")[1])+1))
                    pAlignPos.append(entries[i].split(":")[5])
                    totalAlignLength = totalAlignLength + \
                        int(entries[i].split(":")[4])

                print("\t".join([gene, target, match, str(len(entries)), ",".join(orfIDs), ",".join(orfPos), ",".join(orfLengths),
                                 ",".join(orfSeqIdents), str(MinSeqIdent), str(MaxSeqIdent), ",".join(pAlignPos), str(totalAlignLength)]))


# read in fasta file
if len(sys.argv) < 2:
    print("usage DetectValidSplitOrfMatches_py3.py BlastpAlign.out")
else:

    file = open(sys.argv[1], 'r')

    # variable that remembers the last transcript in use (all alignments are sorted by target sequence)
    lastElem = "dummy|dummy"
    Alignments = {}  # hash map that stores all valid alignments for one target sequence
    # print header for result file
    print("\t".join(["geneID", "targetTransID", "OrfTransID", "NumOrfs", "OrfIDs",
                     "OrfPos", "OrfLengths", "OrfSeqIdents", "MinSeqIdent", "MaxSeqIdent", "protAlignPos", "protCoverage"]))

    for line in file:
        elems = line.split("\t")
        orf = re.split("[|:]+", elems[0])

        # print(lastElem + "***" + elems[1])
        if lastElem != elems[1]:
            # found a new target transcript
            # go through all transcripts that have matches to the target
            target = lastElem.split("|")
            checkAlignments(Alignments, target[0], target[1])
            lastElem = elems[1]
            Alignments = {}
            target = (elems[1]).split("|")

        orflength = orf[5].split("-")
        blastOrfLength = int(orflength[1])-int(orflength[0])
        # why +2 here? The format is one based so half-open, that means adding +1 should be sufficient!
        orfLenNuc = int(orf[4]) - int(orf[3])+1
        orfLenProt = orfLenNuc/float(3)
        alignLength = float(int(elems[9])-int(elems[8])+1)

        # and (orf[0] == target[0])
        if (float(elems[2]) >= identityCutoff) and (orfLenProt >= minLength) and ((alignLength/blastOrfLength) >= minAlignmentRate) and (orf[0] == target[0]) and (alignLength > minAlignLength):
            dummy = orf[2:5]
            dummy.append(str(elems[2]))
            dummy.append(str(int(alignLength)))
            dummy.append("-".join([elems[8], elems[9]]))
            if (orf[1] in Alignments):  # transcript has at least one matching orf
                entries = Alignments[orf[1]]
                add = True
                for i in range(0, len(entries)):
                    if (orf[2] == entries[i].split(":")[0]):
                        add = False
                        # Only the longest Alignment of an ORF against a Transcript is used
                        if (int(alignLength) > int(entries[i].split(":")[4])):
                            Alignments[orf[1]][i] = colon.join(dummy)
                if add:
                    Alignments[orf[1]].append(colon.join(dummy))
            else:
                Alignments[orf[1]] = [colon.join(dummy)]

        target = elems[1].split("|")
    # output the last transcript
    checkAlignments(Alignments, target[0], target[1])


# ENSMUSG00000031748|ENSMUST00000127900:ORF-189:67160:67360  ENSMUSG00000000001|ENSMUST00000000001
# 83.87   31      4       1       102     132     143     172     2e-06   51.2
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
