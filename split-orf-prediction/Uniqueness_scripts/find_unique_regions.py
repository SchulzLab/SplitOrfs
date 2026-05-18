import sys
from Bio import SeqIO
import pybedtools

# Arguments:
#
# 1 = Mummer.mums
# 2 = non_unique.bed (will be created)
# 3 = reference.fasta
# 4 = reference.bed (will be created)
# 5 = unique.bed (will be created)

# read in the result file of MUMMER (sys.argv[1])
infile = open(sys.argv[1], "r")
lines = infile.readlines()

# create a bedfile (sys.argv[2]) using the result file of MUMMER to annotate the matched regions:
# Format:   Transript-ID    start   end
with open(sys.argv[2], "w") as out1:
    name = ""
    for i in lines:
        if i[0] == ">":
            name = i[2:-1]
        else:
            columns = i.split()
            start = int(columns[2])-1#1 has to be substracted as Mummers first position is 1 and bedfile positions start at 0
            # MUMMER gives the start position and length of the match, therefor the end position needs to be calculated
            end = start+int(columns[3])
            out1.write(name + "\t" + str(start) + "\t" + str(end) + "\n")

# create a bedfile (sys.argv[4]) of the split-ORF protein fasta (sys.argv[3])
# by simply annotating the whole sequence for each transcript
# Format:   Transcript ID   start (will always be 0)    end
reference = SeqIO.parse(sys.argv[3], "fasta")
with open(sys.argv[4], "w") as out2:
       for line in reference:
           out2.write(line.id + "\t" + "0" + "\t" + str(len(str(line.seq))) + "\n")

# Use the bedfile of the split-ORF proteins (sys.argv[4]) as reference and use
# the inverse intersect function of bedtools to create a bedfile (sys.argv[5])
# that now contains only the unique regions. (The regions not annotated in (sys.argv[4]))
a = pybedtools.BedTool(sys.argv[2])
b = pybedtools.BedTool(sys.argv[4])
b.subtract(a).saveas(sys.argv[5])
