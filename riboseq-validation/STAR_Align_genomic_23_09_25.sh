#----- Aligns the Ribo-seq reads to the genome using STAR ----- #
#----- in -i mode performs the indexing ----- #
#----- in -a mode performs the alignment ----- #
#!/bin/bash
# Help message:
usage="
Usage: ./STAR_Align_genomic_23_09_25.sh [-options] [arguments]

where:
-h			show this help
-i number_of_threads, star_index, genome_fasta and genome_gtf
-a number_of_threads, star_index, ribo_reads, bam_file, align_ends_type"


RED='\033[0;31m' # Red colour for error messages
NC='\033[0m'	 # No colour for normal messages


# available options for the programm
while getopts ':hi:a:' option; do
  case "$option" in
    h) 
        echo "$usage"
        exit 1
        ;;
    i)
        run_index=true
        if [ "$#" -lt 5 ]; then
            echo "Error: -i index requires 5 arguments: number_of_threads, star_index, genome_fasta and genome_gtf" >&2
            exit 1
        fi

        # Assign them
        number_of_threads=$2
        star_index=$3
        genome_fasta=$4
        genome_gtf=$5

        ;;

    a)
      run_align=true
      args=("$OPTARG" "${@:OPTIND:6}")   # get 4 args (1 from OPTARG, 3 more)
      if [ "${#args[@]}" -ne 6 ]; then
        echo "Error: -a requires 5 arguments: number_of_threads, star_index, ribo_reads, bam_file, align_ends_type" >&2
        exit 1
      fi

      number_of_threads=${args[0]}
      star_index=${args[1]}
      ribo_reads=${args[2]}
      bam_name=${args[3]}
      align_ends_type=${args[4]}
      genome_fasta=${args[5]}
      ;;
   \?) 
        printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1
        ;;
  esac
  shift 2
done


if [ "$run_index" = true ]; then
        # index generation
        STAR --runThreadN 50 --runMode genomeGenerate --genomeDir "$star_index" --genomeFastaFiles ${genome_fasta}\
         --sjdbGTFfile ${genome_gtf} --sjdbOverhang 49
        # --sjdbOverhang 49: maxreadlength - 1: 50bp length limit for Ribo-seq data Vlado
fi

if [ "$run_align" = true ]; then

    out_path=$(dirname $bam_name)


    if [ ! -s  $out_path/genome_file.txt ]; then
        cut -f1,2 ${genome_fasta}.fai >\
        $out_path/genome_chrom_ordering.txt
        chmod 777 $out_path/genome_chrom_ordering.txt
    fi

    if [ ! -e "${bam_name}_sorted.bam" ]; then
      sample=$(basename $bam_name)
      # align ribo_reads against the genome
      STAR\
      --runThreadN $number_of_threads\
      --alignEndsType ${align_ends_type} \
      --outSAMstrandField intronMotif\
      --alignIntronMin 20\
      --alignIntronMax 1000000\
      --genomeDir $star_index\
      --readFilesIn $ribo_reads\
      --twopassMode Basic\
      --seedSearchStartLmax 20\
      --seedSearchStartLmaxOverLread 0.5\
      --outFilterMatchNminOverLread 0.9\
      --outSAMattributes All\
      --outSAMtype BAM SortedByCoordinate\
      --outFileNamePrefix ${bam_name}

      bam_file=${bam_name}Aligned.sortedByCoord.out.bam

      filtered_bam_file=${bam_name}_filtered.bam
      samtools view -F 256 -F 2048 -q 10 -b $bam_file > $filtered_bam_file
      sorted_bam_file=${bam_name}_sorted.bam
      samtools sort -o $sorted_bam_file $filtered_bam_file
      samtools index -@ 10 $sorted_bam_file
      bed_file=${bam_name}.bed
    else 
      sample=$(basename $bam_name)
      sorted_bam_file=${bam_name}_sorted.bam
      bed_file=${bam_name}.bed
    fi

    present_chromosomes=${bam_name}_chromosomes.txt
    samtools view -H $sorted_bam_file | grep '@SQ' | cut -f 2 | cut -d ':' -f 2  | sort | uniq > $present_chromosomes

    if [ ! -e "${bam_name}_chrom_sort.bed" ]; then
      # echo "converting bam to bed"
      bedtools bamtobed -i $sorted_bam_file -split > $bed_file
      sorted_bed_file=$out_path/$(basename $bed_file .bed)_chrom_sort.bed
      sort -k1,1 -k2,2n $bed_file > $sorted_bed_file
    fi



    # subset the genome bed_file to the present genomes
    grep -Fwf $present_chromosomes $out_path/genome_chrom_ordering.txt | sort -k1,1 -k2,2n > $out_path/genome_chrom_ordering_${sample}.txt
    # -F: no regex, take chrs literally
    # -w: match the whole word, e.g. 1 does not match 10, 11 etc but only 1
    # -f: file input of the pattern that are searched for
fi









