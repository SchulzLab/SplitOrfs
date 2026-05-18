#!/usr/bin/python
# create bed file for Ensmebl dump from Biomart

# When downloading PFAM (or similar) annotations for a species the file looks as below (except the header was already modified)
# This script converts that into proper bed format to be used with bedtools as part of the annotation

# GeneID  TranscriptID    PfamID  start   end
# ENSMUSG00000046716      ENSMUST00000226237      PF03402 14      224
# ENSMUSG00000090546      ENSMUST00000166381
# ENSMUSG00000092544      ENSMUST00000149782      PF01352 79      120


# conversion line for Ensembl putput
# awk 'BEGIN{OFS="\t"}{print $2,$4,$5,$3}' EnsemblMouse95PFAMannotation.txt > PfamEnsemblMouse.bed

import sys
import re

# read in fasta file
if len(sys.argv) < 2:
    print("usage convertEnsemblOutput2Bed.py EnsemblBiomartDownload.txt")
else:

    file = open(sys.argv[1], 'r')
    # skip header
    next(file)
    for line in file:
        line = line.rstrip()
        elems = line.split("\t")
        if len(elems) > 2:
            print("\t".join([elems[1], elems[3], elems[4], elems[2]]))
