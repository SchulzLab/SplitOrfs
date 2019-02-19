#This script is the master script to execute the SplitOrfs workflow.


output=$1
proteins=$2
transcripts=$3
annotations=$4



if [[ ($# -ne 3 && $# -ne 4) ]]; then
        echo "usage runSplitOrfs.sh outputFolder proteins.fa transcripts.fa annotation.bed"
    	echo  "output Folder will be created if not existing"
else 
	echo "run the SplitOrfs script on: " $output $proteins $transcripts $annotations

	#create output folder if it does not exist
	mkdir -p $output
	#create Orf sequences
	python OrfFinder.py ${transcripts} > ${output}/OrfProteins.fa
	makeblastdb -in ${proteins} -out ${output}/ProteinDatabase -dbtype prot

	#use BlastP to align the translated ORFs to the proteins (currently using 20 threads, change -num_threads otherwise
	#-outfmt keyword standard results in file format:
	#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

	blastp -query ${output}/OrfProteins.fa -db ${output}/ProteinDatabase -out ${output}/OrfsAlign.txt -outfmt "6 std" -evalue 10 -num_threads 20

	#sort the blastp output by the second column to group by proteins
	sort -k2 ${output}/OrfsAlign.txt > ${output}/OrfsAlign_sorted.txt
	#rm ${output}/OrfsAlign.txt

	#run the detection script to parse
	python DetectValidSplitOrfMatches.py ${output}/OrfsAlign_sorted.txt > ${output}/ValidProteinORFPairs.txt

	#sort file per Orf-transcript ID on column 3. Here it is important to omit the head while sorting
	cat ${output}/ValidProteinORFPairs.txt | awk 'NR<2{print ;next}{print | "sort -k3"}'  > ${output}/ValidProteinORFPairs_sortCol3.txt
	python getLongestOrfMatches.py ${output}/ValidProteinORFPairs_sortCol3.txt > ${output}/UniqueProteinORFPairs.txt

	python makeBed.py ${output}/UniqueProteinORFPairs.txt > ${output}/UniqueProteinMatches.bed

fi

if [[ "$#" -eq 4 ]]; then
	#only executed if you give a annotation.bed file with the script
	echo "Run annotation part"
        #if you want to run this for other bedfiles than the ones provided in the annotations folder look into convertEnsemblOutput2Bed.py
		#python convertEnsemblOutput2Bed.py EnsemblDOwnloadannotation.txt > EnsemblMouse95PFAMannotation.bed


		# Each ORF match has to overlap  by 100% of one annotated entry to be considered by bedtools
		bedtools intersect -a ${output}/UniqueProteinMatches.bed -b ${annotations} -wa -F 1  -wb > ${output}/intersectResults.txt
		python addFunctionalOverlap.py ${output}/UniqueProteinORFPairs.txt ${output}/intersectResults.txt > ${output}/UniqueProteinORFPairs_annotated.txt

fi	
