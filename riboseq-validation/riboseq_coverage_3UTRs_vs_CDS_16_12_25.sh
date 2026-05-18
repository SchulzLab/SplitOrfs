#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate Riboseq


# Help message:
usage="
Usage: ./filter_intersection_pipeline.sh [-options] [arguments]

where:
-h			show this help
"

# available options for the programm
while getopts 'b:c:hp:s:t:' option; do
  case "$option" in
    b)
        bam="$OPTARG"
        ;;
    c)
        cds_coordinates="$OPTARG"
        ;;  
    p)
        script_path="$OPTARG"
        ;;
    s)
        sample_name="$OPTARG"
        ;;
    t)
        three_primes="$OPTARG"
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

# ################################################################################
# # READ AND CHECK ARGUMENTS                                                     #
# ################################################################################

# bam="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome/ERR3367797/ERR3367797_NMD_sorted.bam"
# # this should be the CDS coordinates of the transcripts that I also have for the 3' UTRs, but should double check
# cds_coordinates="/projects/splitorfs/work/Riboseq/data/region_input/genomic/Ens_110_CDS_coordinates_genomic_protein_coding_tsl_refseq_filtered.bed"
# sample_name="ERR3367797"
# # these are three prime UTRs filtered for untranslated genes, but not yet for TPM
# three_primes="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/NMD_genome/ERR3367797/3_primes_genomic_merged_numbered_ERR3367797.bed"

output_dir=$(dirname "$three_primes")
window_bed_three_primes="${output_dir}/$(basename $three_primes .bed)_windows.bed"
window_bed_cds_coordinates=$(dirname "$cds_coordinates")/$(basename "$cds_coordinates" .bed)_windows.bed

bedtools_coverage_three_prime_outfile="${output_dir}/$(basename $window_bed_three_primes .bed)_coverage.tsv"
bedtools_coverage_cds_outfile="${output_dir}/$(basename $window_bed_cds_coordinates .bed)_coverage.tsv"


if [ ! -s "${window_bed_three_primes}" ]; then
    bedtools makewindows -b "${three_primes}" -w 20 -i srcwinnum > "${window_bed_three_primes}"
    bedtools coverage -F 0.33 -split -a "${window_bed_three_primes}" \
    -b "$bam" > "${bedtools_coverage_three_prime_outfile}"
fi


echo "${window_bed_cds_coordinates}"
# for the CDS all samples may use the same windows, unless we want to filter these as well for the untranslated genes
if [ ! -s "${window_bed_cds_coordinates}" ]; then
    bedtools makewindows -b "${cds_coordinates}" -w 20 -i srcwinnum > "${window_bed_cds_coordinates}"
    bedtools coverage -F 0.33 -split -a "${window_bed_cds_coordinates}" \
    -b $bam > "${bedtools_coverage_cds_outfile}"
fi


if [ ! -s "${output_dir}/3_primes_filtered_for_CDS_distribution_${sample_name}.bed" ]; then
    python "${script_path}/cds_vs_three_prime_coverage_windows_empirical.py" \
    "${bedtools_coverage_three_prime_outfile}" \
    "${bedtools_coverage_cds_outfile}" \
    "${sample_name}" \
    "${three_primes}" \
    0.05
fi

sort -k1,1 -k2,2n "${output_dir}/3_primes_filtered_for_CDS_distribution_${sample_name}.bed" \
> "${output_dir}/3_primes_filtered_for_CDS_distribution_${sample_name}_sorted.bed"

# merge the overlapping 3'UTR regions to not have duplicates: when smapling the same region can only be selected once
bedtools merge -i  "${output_dir}/3_primes_filtered_for_CDS_distribution_${sample_name}_sorted.bed" -s -c 4,5,6 -o collapse,min,distinct \
  > "${output_dir}/3_primes_filtered_for_CDS_distribution_${sample_name}_merged.bed"
 