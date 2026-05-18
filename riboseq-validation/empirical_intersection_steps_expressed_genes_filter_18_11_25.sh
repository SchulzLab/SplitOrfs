#!/bin/bash
#Help message:
usage="
Usage: 

options...

where:
-h			show this help
"

#Check that all arguments are provided and are in the correct file format. Give error messages according to errors in the call and exit the programm if something goes wrong
RED='\033[0;31m' #Red colour for error messages
NC='\033[0m'	 #No colour for normal messages


#available options for the programm
while getopts ':h' option; do
  case "$option" in
    h) echo "$usage"
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
  shift 2
done

if [[ $# -ne 7 ]]; then #check for right number of arguments
  echo -e "${RED}
  ERROR while executing the script!
  Wrong number of arguments.${NC}"
    echo "Supplied arguments: $@" >&2 
    echo "$usage" >&2
    exit 1
fi

################################################################################
# READ AND CHECK ARGUMENTS                                                     #
################################################################################

input_file=$1
unique_regions=$2
coordinates_3_prime=$3
outname=$4
random_region_path=$5
genome_fasta=$6
tmp_dir=$6
script_path=$7

# echo "$input_file"
# echo $unique_regions
# echo $coordinates_3_prime
# echo $outname
# echo $random_region_path



sample=$(basename "$outname")
out_path=$(dirname "$outname")
# echo $out_path

# mkdir $outname


################################################################################
# GENERATE BED FILE IF NECESSARY                                              #
################################################################################

if [[ "$input_file" == *.bed ]]; then
  sorted_bedfile="$input_file"
  present_chromosomes="$out_path/$(basename $outname)_chromosomes.txt"
  samtools view -H "$(dirname "$input_file")/$(basename "$sorted_bedfile" _chrom_sort.bed)_sorted.bam" | grep '@SQ' | cut -f 2 | cut -d ':' -f 2  | sort | uniq > "$present_chromosomes"
else
    echo "generating bed file"
    bed_file="$out_path/${sample}.bed"
    sorted_bedfile="$out_path/${sample}_chrom_sort.bed"
    echo "$bed_file"
    echo "$sorted_bedfile"
    bedtools bamtobed -i "$input_file" -split > "$bed_file"
    sort -k1,1 -k2,2n "$bed_file" > "$sorted_bedfile"
    rm "$bed_file"

    present_chromosomes="$out_path/$(basename $outname)_chromosomes.txt"
    samtools view -H "$input_file" | grep '@SQ' | cut -f 2 | cut -d ':' -f 2  | sort | uniq > "$present_chromosomes"
fi


################################################################################
# GENERATE REF FILES IF NECESSARY                                              #
################################################################################
if [ ! -s  "$random_region_path/genome_file.txt" ]; then
    cut -f1,2 "${genome_fasta}.fai" >\
    "$random_region_path/genome_chrom_ordering.txt"
    chmod 777 "$random_region_path/genome_chrom_ordering.txt"
fi


# chromosomes present in the unique regions
if [ ! -s  "$random_region_path/chromosomes_unique_regions.txt" ]; then
    cut -f1 "$unique_regions" | sort | uniq > "$random_region_path/chromosomes_unique_regions.txt"
fi

# sort unique regions according to chromosomes
sorted_unique_regions="$(dirname $unique_regions)/$(basename $unique_regions .bed)_chrom_sorted.bed"
if [ ! -s  "$sorted_unique_regions" ]; then
    sort -k1,1 -k2,2n "$unique_regions" > "$sorted_unique_regions"
fi


################################################################################
# INDEX, SORT, MAKE BED FROM GENOMIC ALIGNMENTS                                #
################################################################################
if [ ! -s  "$out_path/genome_chrom_ordering_$(basename $outname).txt" ]; then
  present_chromosomes="$out_path/$(basename $outname)_chromosomes.txt"
  grep -Fwf "$present_chromosomes" "$random_region_path"/genome_chrom_ordering.txt | sort -k1,1 -k2,2n > "$out_path/genome_chrom_ordering_$(basename $outname).txt"
fi


echo "intersecting with unique regions"
intersectBedfile="${outname}"_intersect_counts_sorted.bed




################################################################################
# INTERSECT WITH UNIQUE REGIONS                                                #
################################################################################
bedtools intersect\
   -s\
   -wao\
   -a "$sorted_unique_regions"\
   -b "$sorted_bedfile"\
   -sorted\
   -g "$out_path/genome_chrom_ordering_$(basename $outname).txt"\
   | sort -T "${tmp_dir}" -nr -k13,13\
   > "$intersectBedfile"

#/scratch/tmp/$USER



################################################################################
# INTERSECT WITH BACKGROUND REGIONS                                            #
################################################################################

for i in {1..20}; do
  randomfile="${outname}_random_background_regions_${i}.bed"
  sorted_randomfile="${outname}_random_background_regions_sorted_${i}.bed"
  if [ ! -s  $sorted_randomfile ]; then
    python "${script_path}"/BackgroundRegions_bed_genomic_fix.py\
    "$sorted_unique_regions"\
    "$coordinates_3_prime"\
    "$randomfile"\
    "$i"

    # sort the randomfile with dummy entries
    sort -T "${tmp_dir}" -k1,1 -k2,2n "$randomfile" > "$sorted_randomfile"
    echo "$sorted_randomfile"
    # rm $randomfile
  fi


  # debug mode on
  # set -x 
  randomintersectfile="${outname}_${i}_random_intersect_counts.bed"
  bedtools intersect\
    -s\
    -wao\
    -sorted\
    -g "$out_path/genome_chrom_ordering_$(basename $outname).txt"\
    -a "$sorted_randomfile" \
    -b "$sorted_bedfile"\
    | sort -T "${tmp_dir}" -nr -k13,13\
    > $randomintersectfile


    # set +x 
    # debug mode off
done





