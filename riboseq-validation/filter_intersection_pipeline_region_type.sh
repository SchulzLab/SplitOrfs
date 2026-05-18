#----- ----- #

#!/bin/bash -l

eval "$(conda shell.bash hook)"
source activate Riboseq


# Help message:
usage="
Usage: ./filter_intersection_pipeline_region_type.sh [-options] [arguments]

where:
-h			show this help
"

# available options for the programm
while getopts 'b:c:e:g:hi:n:o:p:r:s:t:u:d' option; do
  case "$option" in
    b)
        bam="$OPTARG"
        ;;
    c)
        cds_coordinates="$OPTARG"
        ;;
    d)
        dup=true
        ;;
    e)
        ensembl_gtf="$OPTARG"
        ;;
    g)
        genome_fasta="$OPTARG"
        ;;
    i)
        intersection_input="$OPTARG"
        ;;
    n) 
        name="$OPTARG"
        ;;
    o)
        output_star="$OPTARG"
        ;;
    p)
        script_path="$OPTARG"
        ;;
    r)
        region_type="$OPTARG"
        ;;
    s)
        sample_name="$OPTARG"
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

outdir="$output_star/${region_type}_genome/${sample_name}"
if [[ ! -d "${outdir}" ]];then
    mkdir "${outdir}"
fi



if [[ ! -e  "${outdir}"/"${sample_name}_${region_type}_intersect_counts_sorted.bed" ]]; then
    if [[ ! -e  "${outdir}/${sample_name}_${region_type}_htseq_counts.tsv" ]]; then
        htseq-count -f bam --secondary-alignments ignore \
        -c  "${outdir}/${sample_name}_${region_type}_htseq_counts.tsv" \
        --supplementary-alignments ignore "$bam" "$ensembl_gtf"
    fi


    threeprime_basename=$(basename "$three_primes" .bed)
    if [[ ! -e "${outdir}/${threeprime_basename}_${sample_name}.bed" ]]; then
        # echo "$three_primes"
        # echo "${outdir}/${sample_name}_${region_type}_htseq_counts.tsv"
        # echo "$ensembl_gtf"
        python "${script_path}"/filter_bed_file_for_expressed_genes_rnanrom.py \
            "$three_primes" \
            "${outdir}/${sample_name}_${region_type}_htseq_counts.tsv" \
            "$ensembl_gtf" \
            20
    fi

    if [[ ! -e  "${outdir}/Unique_DNA_Regions_genomic_final_${sample_name}.bed" ]]; then
    echo "$unique_region_dir/Unique_DNA_Regions_genomic_final.bed"
    echo "${outdir}/${sample_name}_${region_type}_htseq_counts.tsv"
    echo "$ensembl_gtf"
        python "${script_path}"/filter_bed_file_for_expressed_genes_rnanrom.py \
            "$unique_region_dir/Unique_DNA_Regions_genomic_final.bed" \
            "${outdir}/${sample_name}_${region_type}_htseq_counts.tsv" \
            "$ensembl_gtf" \
            20
    fi
        
    cds_coordinates_tpm_filtered="${outdir}/"$(basename "${cds_coordinates}" .bed)_"${sample_name}".bed
    if [[ ! -e "${cds_coordinates_tpm_filtered}" ]]; then
        python "${script_path}"/filter_bed_file_for_expressed_genes_rnanrom.py \
            "${cds_coordinates}" \
            "${outdir}/${sample_name}_${region_type}_htseq_counts.tsv" \
            "$ensembl_gtf" \
            20
    fi

    bash "${script_path}"/riboseq_coverage_3UTRs_vs_CDS_16_12_25.sh -b "${bam}" -c "${cds_coordinates_tpm_filtered}" \
    -s "${sample_name}" -t "${outdir}"/"${threeprime_basename}"_"${sample_name}".bed -p "${script_path}"

    # implement a filter that only conducts the empiricial intersection pipeline
    # if a certain read depth is found

    # this file should contain all reads mapping to mRNA
    # keep sample if the total count ampping to mRNA is larger than a certain number
    if [[ $dup == true ]]; then
        keep_sample=$(python "${script_path}"/check_seq_depth.py "${outdir}/${sample_name}_${region_type}_htseq_counts.tsv" 5)
        if [[ "$keep_sample" = "True" ]]; then
            keep_sample=true
        else
            keep_sample=false
        fi
    else
        keep_sample=true
    fi
    
    if [[ $keep_sample == true  ]]; then
        echo "$sample_name"
        "${script_path}"/empirical_intersection_steps_expressed_genes_filter_18_11_25.sh  \
            "$intersection_input" \
            "${outdir}"/Unique_DNA_Regions_genomic_final_"${sample_name}".bed \
            "${outdir}"/3_primes_filtered_for_CDS_distribution_"${sample_name}"_merged.bed \
            "${outdir}/${sample_name}_${region_type}" \
            "${outdir}" \
            "${genome_fasta}" \
            "${script_path}"

        echo "===================       Sample "$sample_name" intersected"
    fi
fi

