#!/bin/bash

# ----- This script is the master script to run the extended split-ORF pipeline. The Inputs are described below,  ----- #
# ----- the Outputs are all stored in a new folder which is created at the given location. The pipeline predicts  ----- #
# ----- potential Split-ORF transcripts based on the supplied sequences of the transcripts of interest            ----- #
# ----- unique regions of the potential Split-ORFs which can be used for Ribo-seq validation are calculated       ----- #

# ----- Help message: ----- #
usage="
Usage: ./run_splitorfs_pipeline.sh [-h] json-file \\
  

Arguments in JSON file:
  file_path
    Path to directory where Input files are stored.
    All Input files should be in the same directory.
  proteins
      A multi-FASTA file containing the amino acid sequences of the proteins
      that are used as references (e.g. protein coding sequences from Ensembl 
      with certain transcript support level).

  transcripts
      A multi-FASTA file containing the DNA sequences of the reads/transcripts
      that shall be analyzed (input transcripts).

  annotation
      A TAB delimited file containing the PFAM annotations for the used genome build.
      These were downloaded from Ensembl and modified with a custom script.
      Format:
        GeneID TranscriptID PFAM ID PFAM start PFAM end
        and remodel with Input_scripts/convert_ensmebl_output_to_bed.py

  reference_transcripts
      A multi-FASTA file of all reference transcripts (should not contain the 
      input transcripts).

  exon_positions
      A BED file containing the chromosome-specific positions of all exons.
      These can be downloaded from Biomart using the structures functionality.
      Format:
        downlaod from Ensembl with Chromosome/scaffold name Gene stable ID
         Transcript stable ID Genomic coding start  Genomic coding end  Strand
        and remodel with Input_scripts/CDS_coordinates_to_bedfile.py

  align_method
      blast or diamnond: diamond is faster while blast recovers more hits

  cds_coordinate_bed
      genomic coordinates of the protein coding CDS reference regions
      Format:
        Gene stable ID	Transcript stable ID	Exon region start (bp)	
        Exon region end (bp)	Transcript start (bp)	Transcript end (bp)
        	Strand	Chromosome/scaffold name

  output_dir
      Directory where Output folder with time stamp will be placed
      Please use absolute paths!!!

Options:
  -h    Show this help message.
"


# ----- available options for the programm -----#
while getopts ':h' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done

# ----- Check that all arguments are provided and are in the correct file format.                         ----- #
# ----- Give error messages according to errors in the call and exit the programm if something goes wrong ----- #
RED='\033[0;31m' #Red colour for error messages
NC='\033[0m'	 #No colour for normal messages

# ----- check for right number of arguments ----- #
if [ "$#" -ne 1 ]; then 
  echo -e "${RED}
ERROR while executing the Pipeline!
Wrong number of arguments.${NC}"
  echo "$usage" >&2
  exit 1
fi

# ----- Assign the given arguments to variables ----- #
CONFIG="$1"

echo "CONFIG is: $CONFIG"
ls -l "$CONFIG"

file_path=$(jq -r '.file_path' "$CONFIG")
proteins=$(jq -r '.proteins' "$CONFIG")
transcripts=$(jq -r '.transcripts' "$CONFIG")
annotation=$(jq -r '.annotation' "$CONFIG")
reference_transcripts=$(jq -r '.reference_transcripts' "$CONFIG")
exon_positions=$(jq -r '.exon_positions' "$CONFIG")
align_method=$(jq -r '.align_method' "$CONFIG")
cds_coordinate_bed=$(jq -r '.cds_coordinate_bed' "$CONFIG")
output_dir=$(jq -r '.output_dir' "$CONFIG")


# Construct full paths
proteins="$file_path/$proteins"
transcripts="$file_path/$transcripts"
annotation="$file_path/$annotation"
reference_transcripts="$file_path/$reference_transcripts"
exon_positions="$file_path/$exon_positions"
cds_coordinate_bed="$file_path/$cds_coordinate_bed"

# ----- check if any argument is a directory ----- #
if  [ -d "$proteins" ] || [ -d "$transcripts" ] || [ -d "$annotation" ] || [ -d "$reference_transcripts" ] || [ -d "$exon_positions" ] || [ -d "$align_method" ] || [ -d "$cds_coordinate_bed" ]; then
  echo -e "${RED}
ERROR while executing the Pipeline!
One or more of the arguments are directories.${NC}"
  echo "$usage" >&2
  exit 1
  
# ----- check if every argument has the correct file extension ----- #
elif ! [[ $proteins =~ \.(fa|fasta)$ ]] || ! [[ $transcripts =~ \.(fa|fasta)$ ]] || ! [[ $annotation =~ \.(txt|bed|tsv)$ ]] || ! [[ $reference_transcripts =~ \.(fa|fasta)$ ]] || ! [[ $exon_positions =~ \.(txt|bed|tsv)$ ]]; then
  echo -e "${RED}
ERROR while executing the Pipeline!
One or more of the arguments are not in the specified file format.${NC}"
  echo "$usage" >&2
  exit 1
fi


# ----- get the script directory ----- #
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


# ----- create the directory "Output" if it does not already exist ----- #
if [ -d "$output_dir" ]
then
	echo "Directory ${output_dir} already exists."
else
	echo "Directory ${output_dir} does not exist, creating ..."
	mkdir "${output_dir}"
	echo "Directory ${output_dir} created."
fi

# ----- create a run specific output folder ----- #
timestamp=$(date "+%d.%m.%Y-%H.%M.%S")
if [ -d "${output_dir}/run_$timestamp" ]
then 
	echo "Directory ${output_dir}/run_$timestamp already exists." >&2
	exit 1
else
	mkdir "${output_dir}"/run_$timestamp
  mkdir "${output_dir}"/run_$timestamp/tests
fi

# ----- Create a log file that saves all text outputs of the pipeline ----- #
exec > >(tee -i "${output_dir}/run_$timestamp/Logfile.txt")
exec 2>&1

echo "Log Location should be: [ ${output_dir}/run_$timestamp ]"

output="${output_dir}/run_$timestamp"
echo "*********"${output}"**********"
echo "run the Pipeline with: " "${output}" "${proteins}" "${transcripts}" "${annotation}" "${reference_transcripts}" "${exon_positions}"

# ----- create Orf sequences using the OrfFinder script ----- #
echo "searching for possible ORFs"
python "${script_dir}"/SplitOrfs-master/orf_finder.py "${transcripts}" > "${output}"/OrfProteins.fa

# ----- Determine Unique Protein regions by exact matching with Mummer3 ----- #
# ----- Regions of ORFs not matching exactly ot protein coding reference are unique ----- #
echo "Align predicted ORFs (AA) to protein coding transcripts mummer -maxmatch -l 8"
mummer -maxmatch -l 8 "${proteins}" "${output}"/OrfProteins.fa > "${output}"/Proteins_maxmatch_l8.mums

echo "Select the non matching regions as unique regions and save as bedfile"
python "${script_dir}"/Uniqueness_scripts/find_unique_regions.py \
 "${output}"/Proteins_maxmatch_l8.mums \
 "${output}"/Protein_non_unique.bed \
 "${output}"/OrfProteins.fa \
 "${output}"/OrfProteins.bed \
 "${output}"/Unique_Protein_Regions.bed

python "${script_dir}"/Uniqueness_scripts/filter_small_regions.py\
 "${output}"/Unique_Protein_Regions.bed \
 8 \
 "${output}"/Unique_Protein_Regions_gt8.bed

echo "Select only the non-unique regions for the upcoming blast"
bedtools subtract\
 -a "${output}"/OrfProteins.bed \
 -b "${output}"/Unique_Protein_Regions_gt8.bed > "${output}"/ProteinRegionsforBlast.bed

bedtools getfasta \
 -fi "${output}"/OrfProteins.fa \
 -bed "${output}"/ProteinRegionsforBlast.bed \
 -fo "${output}"/ProteinsforBlast.fa

if [[ "$align_method" == "blast" ]]; then
    # ----- create a blast database for all proteins using BLAST+ ----- #
    echo "Create BlastDB using the protein coding peptide fasta"
    makeblastdb -in "${proteins}" -out "${output}"/ProteinDatabase -dbtype prot

    # ----- use BlastP to align the translated ORFs to the proteins ----- #
    echo "Blast the non-unique regions against the DB"
    blastp -query "${output}"/ProteinsforBlast.fa \
    -db "${output}"/ProteinDatabase \
    -out "${output}"/OrfsAlign.txt \
    -evalue 10 \
    -num_threads 8 \
    -outfmt "6 std"

    # ----- sort the blastp output by the second column to group by proteins and delete unsorted file ----- #
    sort -k2 "${output}"/OrfsAlign.txt > "${output}"/OrfsAlign_sorted.txt
    rm "${output}"/OrfsAlign.txt
elif [[ "$align_method" == "diamond" ]]; then
    #----- create a blast database for all proteins using BLAST+ ----- #
    echo "Create BlastDB using the protein coding peptide fasta"
    diamond makedb --in "${proteins}" -d "${output}"/ProteinDatabase

    # ----- use diamond BlastP to align the translated ORFs to the proteins ----- #
    echo "Blast the non-unique regions against the DB"
    diamond blastp -q "${output}"/ProteinsforBlast.fa \
    -d "${output}"/ProteinDatabase \
    -o "${output}"/OrfsAlign.tsv \
    --evalue 10 \
    --very-sensitive

    # ----- sort the blastp output by the second column to group by proteins and delete unsorted file ----- #
    sort -k2 "${output}"/OrfsAlign.tsv > "${output}"/OrfsAlign_sorted.txt
    rm "${output}"/OrfsAlign.tsv
fi

# ----- run the detection script to identify the valid ORFs ----- #
# ----- at least two ORFs need to map to the same protein coding reference ----- #
echo "Select the valid ORFs"
python "${script_dir}"/SplitOrfs-master/detect_valid_split_orf_matches.py \
"${output}"/OrfsAlign_sorted.txt > "${output}"/ValidProteinORFPairs.txt

# ----- sort file per Orf-transcript ID on column 3, omit the head while sorting ----- #
cat "${output}"/ValidProteinORFPairs.txt | awk 'NR<2{print ;next}{print | "sort -k3"}' \
 > "${output}"/ValidProteinORFPairs_sortCol3.txt

# ----- select the best ORF matches and convert to bedfile ----- #
python "${script_dir}"/SplitOrfs-master/get_longest_orf_matches.py \
 "${output}"/ValidProteinORFPairs_sortCol3.txt \
 > "${output}"/UniqueProteinORFPairs.txt

python "${script_dir}"/SplitOrfs-master/make_bed.py \
 "${output}"/UniqueProteinORFPairs.txt \
 > "${output}"/UniqueProteinMatches.bed

# ----- Add the domain overlap to the selected ORFs        ----- #
echo "Add functional overlap (PFAM)"
#intersect results are empty for the new run...
bedtools intersect -a "${output}"/UniqueProteinMatches.bed \
 -b "${annotation}" \
 -wa \
 -F 1  \
 -wb > "${output}"/intersectResults.txt

python "${script_dir}"/SplitOrfs-master/add_functional_overlap.py \
 "${output}"/UniqueProteinORFPairs.txt \
 "${output}"/intersectResults.txt \
 > "${output}"/UniqueProteinORFPairs_annotated.txt

# ----- create a report displaying the results of the original pipeline steps ----- #
echo "Create Split-ORF report"

##need to remove . in front of the reltive path to the transcripts file, otherwise this file cannot be found when the path is put together
if [[ "${transcripts}" == ./* ]]; then
  transcripts_R="${transcripts:1}"
else
  transcripts_R="${transcripts}"
fi


R -e "rmarkdown::render('${script_dir}/Split-ORF_Report.Rmd',
output_file='${output}/Split-ORF_Report.html',
params=list(args = c('${output}/ValidProteinORFPairs_sortCol3.txt',
'${output}/UniqueProteinORFPairs_annotated.txt',
'$transcripts_R',
knit_root_dir='${output}')))"


# ----- extract the Valid ORF sequences from the given transcripts ----- #
echo 'Select Split-ORF DNA Sequences'
# now look back at all ORFS, bc until now the unique regions were omitted
# when valid (unique best match) -> take the whole ORF in the new file: Valid_ORF_Proteins.bed
python "${script_dir}"/Uniqueness_scripts/select_valid_orf_sequences.py \
 "${output}"/UniqueProteinORFPairs.txt \
 "${output}"/OrfProteins.bed \
 "${output}"/Valid_ORF_Proteins.bed

# via transcript IDs: select the valid transcript DNA sequences by knowing valid IDs 
python "${script_dir}"/Uniqueness_scripts/select_valid_orf_dna_sequences.py \
 "${transcripts}" "${output}"/Valid_ORF_Proteins.bed \
 "${output}"/ValidORF_DNA_Sequences.fa 

# ----- Determine Unique DNA regions by exact matching of input transcripts to reference transcripts ----- #
# ----- with default paramters (-l = 20) for the minimum matchlength  ----- #
echo "Align ORF-transcripts(DNA) to protein coding transcripts with mummer -maxmatch"
mummer -maxmatch "${reference_transcripts}" "${output}"/ValidORF_DNA_Sequences.fa > "${output}"/DNA_maxmatch.mums

echo "Select the non matching regions as unique regions and save as bedfile"
python "${script_dir}"/Uniqueness_scripts/find_unique_regions.py \
 "${output}"/DNA_maxmatch.mums \
 "${output}"/DNA_non_unique.bed \
 "${output}"/ValidORF_DNA_Sequences.fa \
 "${output}"/ValidORF_DNA_Sequences.bed \
 "${output}"/Unique_DNA_Regions.bed

# ----- Filter out unique regions smaller than the Mummer3 length parameter  ----- #
python "${script_dir}"/Uniqueness_scripts/filter_small_regions.py \
 "${output}"/Unique_DNA_Regions.bed \
 20 \
 "${output}"/Unique_DNA_Regions_gt20.bed

# ----- Extract the valid ORF-Proteins annotated in UniqueProteinORFPairs.txt  ----- #
echo "Extract only the valid regions"
python "${script_dir}"/Uniqueness_scripts/select_valid_orf_sequences.py \
 "${output}"/UniqueProteinORFPairs.txt \
 "${output}"/Unique_Protein_Regions_gt8.bed \
 "${output}"/Unique_Protein_Regions_gt8_valid.bed

# ----- Filter the protein and DNA unique regions for weird SplitORF positions
python "${script_dir}"/Uniqueness_scripts/filter_unique_regions.py \
 "${output}"/UniqueProteinMatches.bed \
 "${output}"/Unique_Protein_Regions_gt8_valid.bed\
 protein\
 "${output}"

python "${script_dir}"/Uniqueness_scripts/filter_unique_regions.py \
 "${output}"/UniqueProteinMatches.bed \
 "${output}"/Unique_DNA_Regions_gt20.bed \
 DNA \
 "${output}"


echo "Finishing Steps"
# ----- get the ORF seqeunces of those predicted proteins with unique regions  ----- #
python "${script_dir}"/Uniqueness_scripts/get_protein_seqs_fasta.py \
 "${output}"/UniqueProteinORFPairs.txt \
 "${output}"/OrfProteins.fa \
 "${output}"/Proteins_for_masspec.fa

# ----- Reorganize Unique_DNA_Regions_gt20.bed for later intersection with riboseq Alignment ----- #
# Note: The position of the unique region is given with respect to the transcript and not the ORF!
python "${script_dir}"/Uniqueness_scripts/bedreorganize_adapted.py \
    "${output}"/Unique_DNA_Regions_gt20_filtered.bed \
    "${output}"/Unique_DNA_Regions_reorganized.bed

python "${script_dir}"/Uniqueness_scripts/bedreorganize_proteins.py \
    "${output}"/Unique_Protein_Regions_gt8_valid_filtered.bed \
    "${output}"/Unique_Protein_Regions_transcript_coords.bed


exon_positions_sorted=$(basename "${exon_positions}" .bed)_sorted.bed
awk '$7 == "1" || $7 == "+"' "${exon_positions}" | sort -k1,1 -k2,2 -k3,3n > "${output}"/plus_strand.bed
awk '$7 == "-1" || $7 == "-"' "${exon_positions}" | sort -k1,1 -k2,2 -k4,4nr > "${output}"/minus_strand.bed
cat "${output}"/plus_strand.bed "${output}"/minus_strand.bed > "${output}"/"${exon_positions_sorted}"
rm "${output}"/plus_strand.bed
rm "${output}"/minus_strand.bed


# ----- Get the genomic positions for the unique DNA and protein regions ----- #
echo "Change the Positions of the unique DNA and Protein as well as the ValidORF bed to their positions within the unspliced transcript"
exon_positions_transcript=$(basename "${exon_positions}" .bed)_transcript_positions.bed
echo "Change the positions to their genomic equivalent and transform into UCSC format"
python "${script_dir}"/Genomic_scripts_18_10_24/exon_to_transcript_positions.py \
 "${output}"/"${exon_positions_sorted}" \
  "${output}"/"${exon_positions_transcript}"


# ----- Test functionality of the exon coordinate conversion ----- #
mkdir "${output}"/tests
python "${script_dir}"/Genomic_scripts_18_10_24/exon_to_transcript_positions.py \
 "${script_dir}"/Genomic_scripts_18_10_24/test/exon_transcript_positions_unit_test.bed \
  "${output}"/tests/exon_transcript_positions_results.bed

python "${script_dir}"/Genomic_scripts_18_10_24/test/test_transcript_exon_positions.py \
 "${output}"/tests/exon_transcript_positions_results.bed \
 "${output}"/"${exon_positions_sorted}" \
 "${output}"/"${exon_positions_transcript}"
 

python "${script_dir}"/Genomic_scripts_18_10_24/genomic_dna_regions_polars.py \
 "${output}"/Unique_DNA_Regions_reorganized.bed \
 "${output}"/"${exon_positions_transcript}" \
 "${output}"/Unique_DNA_Regions_genomic.bed

python "${script_dir}"/Genomic_scripts_18_10_24/genomic_dna_regions_polars.py \
 "${output}"/Unique_Protein_Regions_transcript_coords.bed \
 "${output}"/"${exon_positions_transcript}" \
 "${output}"/Unique_Protein_Regions_genomic.bed

# ----- Get the genomic positions for all valid Split-ORFs ----- #
python "${script_dir}"/Genomic_scripts_18_10_24/create_orf_coords_bed_file.py \
 "${output}"/UniqueProteinORFPairs.txt \
  "${output}"/ORF_transcript_coords.bed


python "${script_dir}"/Genomic_scripts_18_10_24/genomic_dna_regions_polars.py \
 "${output}"/ORF_transcript_coords.bed \
 "${output}"/"${exon_positions_transcript}" \
 "${output}"/ORF_genomic_coords.bed


# ----- Test functionality of unique region to genomic coordinate conversion ----- #
python "${script_dir}"/Genomic_scripts_18_10_24/test/test_transcriptomic_to_genomic_coordinates.py \
 "${script_dir}"/Genomic_scripts_18_10_24/test/test_conversion_with_gen_trans.bed \
 "${output}"/Unique_DNA_Regions_genomic.bed

# ----- Subtract the genomic CDS coordinates from the unique regions ----- #
bedtools subtract -s -a "${output}"/Unique_DNA_Regions_genomic.bed \
 -b "${cds_coordinate_bed}" > "${output}"/Unique_DNA_Regions_genomic_CDS_subtraction.bed

# ----- CDS subtraction for transcriptomic regions and name correction genomic regions ----- #
python "${script_dir}"/Genomic_scripts_18_10_24/genomic_to_transcript_positions.py \
 --ur_dna_regions_genomic "${output}"/Unique_DNA_Regions_genomic_CDS_subtraction.bed \
 --exon_transcript_positions "${output}"/"${exon_positions_transcript}" \
 --out_unique_dna_gen "${output}"/Unique_DNA_Regions_genomic_final.bed \
 --out_unique_dna_trans "${output}"/Unique_DNA_Regions_transcriptomic.bed \
 --out_unique_dna_for_fasta "${output}"/Unique_DNA_Regions_transcriptomic_for_FASTA.bed

# ----- calculate overlap between unique DNA and protein regions genomic         ----- #
echo "calculate genomic overlap between unique DNA and protein regions"
bedtools intersect\
 -a "${output}"/Unique_DNA_Regions_genomic_CDS_subtraction.bed\
 -b "${output}"/Unique_Protein_Regions_genomic.bed\
 > "${output}"/Unique_Regions_Overlap_genomic.bed


# ----- calculate overlap between unique DNA and protein regions transcriptomic         ----- #
echo "calculate transcriptomic overlap between unique DNA and protein regions"
bedtools intersect\
 -a "${output}"/Unique_DNA_Regions_transcriptomic.bed \
 -b "${output}"/Unique_Protein_Regions_transcript_coords.bed \
 > "${output}"/Unique_Regions_Overlap_transcriptomic.bed


# ----- Use bedtools getfasta to extract the fasta sequences of the unique regions annotated in the produced bedfiles ----- #
bedtools getfasta\
 -fi "${output}"/ValidORF_DNA_Sequences.fa\
 -fo "${output}"/Unique_DNA_Regions.fa\
 -bed "${output}"/Unique_DNA_Regions_transcriptomic_for_FASTA.bed

bedtools getfasta\
 -fi "${output}"/OrfProteins.fa\
 -fo "${output}"/Unique_Protein_Regions.fa\
 -bed "${output}"/Unique_Protein_Regions_gt8_valid_filtered.bed



# ----- Create a report with basics statistics of the uniqueness scripts    ----- #
R -e "rmarkdown::render('"${script_dir}"/Extended_Pipeline_new.Rmd',output_file='"${output}"/Uniqueness_Report.html',
params=list(args = c('${output}/Unique_DNA_Regions.fa', 
'${output}/Unique_Protein_Regions.fa',
'${output}/Unique_Protein_Regions_gt8_valid_filtered.bed',
'${output}/UniqueProteinORFPairs.txt',
'${output}/Unique_DNA_Regions_transcriptomic.bed', 
'${output}/Unique_Protein_Regions_transcript_coords.bed',
'${output}/Unique_Regions_Overlap_transcriptomic.bed')))"



# ----- get Split-ORF genes for GO analysis       ----- #
python "${script_dir}"/SplitOrfs-master/get_so_genes_for_go_analysis.py \
 "${output}"/UniqueProteinORFPairs.txt \
 "${output}"/SplitOrfGeneFile.txt 


python "${script_dir}"/SplitOrfs-master/get_background_genes_for_go_analysis.py \
 "${transcripts}" \
 "${output}"/BackgroundGeneFile.txt

 python "${script_dir}"/Uniqueness_scripts/categorize_so_transcripts_by_urs.py \
 --so_results "${output}"/UniqueProteinORFPairs.txt \
 --ur_path "${output}"/Unique_DNA_Regions_genomic_final.bed


# ----- Remove intermediate Results ----- #
mv "${output}"/Unique_Protein_Regions_gt8_valid_filtered.bed "${output}"/Unique_Protein_Regions_protein_coordinates.bed
rm "${output}"/Unique_DNA_Regions_gt20_filtered.bed
rm "${output}"/Unique_DNA_Regions_reorganized.bed
rm "${output}"/Unique_DNA_Regions_gt20.bed
rm "${output}"/ProteinDatabase*
rm "${output}"/intersectResults.txt
rm "${output}"/Unique_Protein_Regions_gt8_valid.bed
rm "${output}"/Unique_Protein_Regions_gt8.bed
rm "${output}"/Unique_DNA_Regions_genomic_CDS_subtraction.bed
rm "${output}"/Unique_DNA_Regions_transcriptomic_for_FASTA.bed
rm "${script_dir}"/Venn*.log