#----- This script maps Ribo-seq data to the genome using the supplied annotation ----- #
# ----- then an intersection with unique regions in genome coords from the split-ORF pipeline ----- #
# ----- is performed as well as background regions of 3' UTRs, an empirical  ----- #
# ----- background distribution is used to determine which unique regions are validated ----- #
# ----- and this is summarized in an Rmd report ----- #

#!/bin/bash -l



# output_star="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter"
# unique_region_dir="/home/ckalk/tools/SplitORF_pipeline/Output/run_13.01.2026-11.44.16_NMD_CDS_subtraction_minAlignLength_15"
# ensembl_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf"
# genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
# The following are the paths to the riboseq reads
# input_fastq_path=(/projects/splitorfs/work/Riboseq/data/fastp/fastp_single_samples/*.fastq)
# file with the 3'UTR background regiosn
# three_primes="/projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic_merged_numbered.bed"

# please note: these are still aligned to Ens110
# would need to realign to respective index with TAMA GTF and also deduplicate
# hypo_ribo_dedup="/projects/splitorfs/work/Riboseq/Output/Michi_Vlado_round_1/alignment_genome/STAR/only_R1/deduplicated"
# upf10_dedup_dir="/projects/splitorfs/work/UPF1_deletion/Output/alignment_genome/STAR/deduplicated"

# cds_coordinates="/projects/splitorfs/work/Riboseq/data/region_input/genomic/Ens_110_CDS_coordinates_genomic_protein_coding_tsl_refseq_filtered.bed"
# script_path="/home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_validation/genomic/resample_random"
# ribo_pipe_path="/home/ckalk/scripts/SplitORFs/Riboseq/Riboseq_analysis_pipeline"

# input_name="NMD"
# region_type="NMD"


# Help message:
usage="
Usage: ./run_Riboseq_validation_pipeline.sh [-options] [arguments]

where:
-h			show this help
"

# available options for the programm
while getopts 'b:c:d:e:g:hi:n:o:p:r:s:t:u:' option; do
  case "$option" in
    b)
        bam_ending="$OPTARG"
        ;;
    c)
        cds_coordinates="$OPTARG"
        ;;
    d)
        dedup="$OPTARG"
        ;;
    e)
        ensembl_gtf="$OPTARG"
        ;;
    g)
        genome_fasta="$OPTARG"
        ;;
    i)
        input_fastq_path="$OPTARG"
        ;;
    n) 
        input_name="$OPTARG"
        ;;
    o)
        output_star="$OPTARG"
        ;;
    p)
        tmp_dir="$OPTARG"
        ;;
    r)
        region_type="$OPTARG"
        ;;
    s)
        script_path="$OPTARG"
        ;;
    t)
        three_primes="$OPTARG"
        ;;  
    u)
        unique_region_dir="$OPTARG"
        ;;   
    h) 
        echo "$usage"
        exit 1
        ;;
   \?) 
        printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1
        ;;
  esac
done

mkdir -p "$output_star"
file_dir="${output_star}/${input_name}_genome"
mkdir -p "$file_dir"
if [ -n "${input_fastq_path}" ]; then
    shopt -s nullglob  # enable nullglob
    input_fastq_files=(${input_fastq_path}/*.fastq)
    if [ ${#input_fastq_files[@]} -eq 0 ]; then
        input_fastq_files=(${input_fastq_path}/*.fq)
        if [ ${#input_fastq_files[@]} -eq 0 ]; then
            echo "Make sure that there are FASTQ files ending in .fq or .fastq in ${input_fastq_path}"  
            echo "Also make sure that the files are not gzipped!"  
            exit 1                       
        fi
    fi
fi



# if [[ -n "${filtered_gtf}" ]]; then
#     # filtering the 3' UTR regions
#     bash "${script_path}/region_handling/get_3prime_genomic_coords.sh" \
#     -c "${cds_coordinates}" \
#     -f "${filtered_gtf}" \
#     -r "${region_type}" \
#     -t "${three_primes}" \
#     -u "${unique_region_dir}"

#     # idea check for if running the filtering step, otherwise, they can be supplied
#     # ready to go, then they are required in BED format: can write tests for this!!!
#     # 3_prime_UTR_coords_genomic_Ensembl_110.txt
#     # and CDS: Ens_110_CDS_coordinates_genomic_all.txt
# fi


if [[ -n "$input_fastq_path" ]]; then
    # Create a Logfile for the alignments in the output directory
    exec > >(tee -i $output_star/AlignmentLogfile.txt)
    exec 2>&1

    echo "STAR index genome"

    if [ ! -d  "$output_star"/index ]; then
        bash "${script_path}/STAR_Align_genomic_23_09_25.sh" \
        -i 50 "$output_star"/index \
        "$genome_fasta" \
        "$ensembl_gtf"
    fi

    echo "Starting alignment against genome"


    # STAR alignment for the samples that are "normal", no UMI deduplication etc
    for i in "${input_fastq_files[@]}"; do
        if [[ "$i" =~ \.fastq$ ]]; then
            sample_name=$(basename "$i" .fastq)
        elif [[ "$i" =~ \.fq$ ]]; then
            sample_name=$(basename "$i" .fq)
        fi
        
        if [[ ! -e  ""${file_dir}"/${sample_name}/${sample_name}_"${input_name}"_chrom_sort.bed" ]]; then
            echo $i
            echo "${file_dir}"/"${sample_name}"
            mkdir -p "${file_dir}"/"${sample_name}"

            bash "${script_path}/STAR_Align_genomic_23_09_25.sh" \
            -a 16 "$output_star"/index "$i" \
            "${file_dir}/${sample_name}/${sample_name}_${input_name}" \
            EndToEnd \
            "$genome_fasta"

            rm "${file_dir}/${sample_name}/${sample_name}"*Aligned.sortedByCoord.out.bam
            rm "${file_dir}/${sample_name}/${sample_name}"*_filtered.bam

            echo "===================       Sample "$sample_name" mapped"
        fi

    done
fi



################################################################################
# Ribo-seq unique region and random region intersection                        #
################################################################################
# Intersection for the "normal" Ribo-seq files
# check in the file dir and subfolders for BAM files
for folder in "$file_dir" "$file_dir"/*/; do
    shopt -s nullglob
    bams=("$folder"/*"$bam_ending")
    if [[ -e "${bams[0]}" ]]; then
        for bam in "${bams[@]}"; do
            sample_name="$(basename "$bam" "$bam_ending")"
                if [ ! -e  "$output_star"/"${region_type}"_genome/"${sample_name}"/Unique_DNA_Regions_genomic_"${sample_name}".bed ]; then
                    mkdir -p "$output_star/"${region_type}"_genome/${sample_name}"
                    if [[ "${dedup}" == "true" ]]; then
                        bash "${script_path}/"filter_intersection_pipeline_region_type.sh \
                            -b "$bam" \
                            -c "$cds_coordinates" \
                            -e "$ensembl_gtf"\
                            -g "$genome_fasta" \
                            -i "${file_dir}/${sample_name}/${sample_name}_${input_name}_chrom_sort.bed" \
                            -n "${file_dir}/${sample_name}/${sample_name}_${input_name}" \
                            -o "$output_star" \
                            -p "$script_path"\
                            -s "$sample_name" \
                            -r "${region_type}" \
                            -t "$three_primes" \
                            -u "$unique_region_dir" \
                            -d
                    else
                        bash "${script_path}/"filter_intersection_pipeline_region_type.sh \
                            -b "$bam" \
                            -c "$cds_coordinates" \
                            -e "$ensembl_gtf"\
                            -g "$genome_fasta" \
                            -i "$bam" \
                            -n "${file_dir}/${sample_name}/${sample_name}_${input_name}" \
                            -o "$output_star" \
                            -p "$script_path"\
                            -s "$sample_name" \
                            -r "${region_type}" \
                            -t "$three_primes" \
                            -u "$unique_region_dir"
                    fi
                fi
        done
    fi
done



