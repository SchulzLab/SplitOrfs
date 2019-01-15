#This script is the master script to execute the NMDOrf workflow.


output=$1
proteins=$2
transcripts=$3

echo "run the NMDOrf script on: " $output $proteins $transcripts

#create output folder if it does not exist
mkdir -p $output
#create Orf sequences
python OrfFinder.py ${transcripts} > ${output}/OrfProteins.fa
makeblastdb -in ${proteins} -out ${output}/ProteinDatabase -dbtype prot

#use BlastP to align the translated ORFs to the proteins (currently using 20 threads, change -num_threads otherwise
#-outfmt keyword standard results in file format:
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

blastp -query ${output}/OrfProteins.fa -db ${output}/ProteinDatabase -out ${output}/OrfsAlign.txt -outfmt "6 std" -evalue 0.001 -num_threads 20

#sort the blastp output by the second column to group by proteins
srt -k2 ${output}/OrfsAlign.txt > ${output}/OrfsAlign_sorted.txt
rm ${output}/OrfsAlign.txt

#run the detection script to parse
python DetectValidNMDOrfMatches.py ${output}/OrfsAlign_sorted.txt > ${output}/ValidProteinORFPairs.txt


