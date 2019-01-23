#This script is the master script to execute the NMDOrf workflow.


output=$1
proteins=$2
transcripts=$3
annotations=$4

if len(sys.argv) < 3:
        print "usage getLongestOrfMatches.py ValidProteinORFPairs_sortCol3.txt"
else :
echo "run the NMDOrf script on: " $output $proteins $transcripts $annotations

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
sort -k2 ${output}/OrfsAlign.txt > ${output}/OrfsAlign_sorted.txt
#rm ${output}/OrfsAlign.txt

#run the detection script to parse
python DetectValidNMDOrfMatches.py ${output}/OrfsAlign_sorted.txt > ${output}/ValidProteinORFPairs.txt

#sort file per Orf-transcript ID on column 3. Here it is important to omit the head while sorting
cat ${output}/ValidProteinORFPairs.txt | awk 'NR<2{print ;next}{print | "sort -k3"}'  > ValidProteinORFPairs_sortCol3.txt

python getLongestOrfMatches.py ValidProteinORFPairs_sortCol3.txt > 


