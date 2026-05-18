import sys
from Bio import SeqIO

transcriptfile = sys.argv[1]
transcript_fastas = SeqIO.parse(open(transcriptfile),'fasta')
transcriptsequences = {}
for fasta in transcript_fastas:
        name, sequence = fasta.id, str(fasta.seq)
        transcriptsequences[name] = sequence

file = open(sys.argv[2],'r')
with open(sys.argv[3],'w') as out:
    for line in file:
        elements = line.strip()
        elements = elements.split("\t")
        temp = elements[0].split(":")
        id = temp[0]
        out.write(">" + id + ":" + temp[1] + ":" + temp[2] + ":" + temp[3] + "\n")
        out.write(transcriptsequences[id][int(temp[2]):int(temp[3])] + "\n")
