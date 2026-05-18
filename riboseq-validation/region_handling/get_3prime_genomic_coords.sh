#!/bin/bash -l

# Help message:
usage="
Usage: ./get_3prime_coords.sh [-options] [arguments]

where:
-c      CDS coordinates downloaded from Ensembl (TXT)
-f      filtered GTF to only keep 3'UTRs of present transcripts
-r      region type (e.g. RI, NMD)
-t      three prime coordinates downloaded from Ensmebl (TXT)
-u      directory that contains the "Unique_DNA_Regions_genomic.bed" file from the Split-ORF pipeline
-h			show this help
"

# available options for the programm
while getopts 'c:f:hr:t:u:' option; do
  case "$option" in
    c)
        cds_coordinates="$OPTARG"
        ;;
    f)
        filtered_gtf="$OPTARG"
        ;;
    r)
        region_type="$OPTARG"
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

output_dir=$(dirname "${three_primes}")
if [ ! -e ""${output_dir}"/3_primes_genomic_merged_numbered.bed" ]; then
  # get BED file from the 3' UTRs downloaded from Ensembl
  python "${script_path}"/region_handling/get_3prime_genomic_coords.py\
  "${filtered_gtf}" \
  "${three_primes}" \
  "${output_dir}"/3_primes_genomic.bed
  # /projects/splitorfs/work/Riboseq/data/region_input/genomic/3_prime_UTR_coords_genomic_Ensembl_110.txt\
  # /projects/splitorfs/work/Riboseq/data/region_input/genomic/3_primes_genomic.bed

  # get BED file of the CDS coordinates downloaded from Ensembl
  python "${script_path}"/region_handling/CDS_coordinates_to_bedfile.py \
  "${cds_coordinates}" \
  "${output_dir}"/CDS_coordinates_genomic_all.bed
  # /projects/splitorfs/work/Riboseq/data/region_input/genomic/Ens_110_CDS_coordinates_genomic_all.txt \
  # /projects/splitorfs/work/Riboseq/data/region_input/genomic/Ens_110_CDS_coordinates_genomic_all.bed

  # sort bedfiles
  sort -k1,1 -k2,2n "${output_dir}"/CDS_coordinates_genomic_all.bed \
  > "${output_dir}"/CDS_coordinates_genomic_all_sorted.bed

  sort -k1,1 -k2,2n "${output_dir}"/3_primes_genomic.bed\
  > "${output_dir}"/3_primes_genomic_sorted.bed

  # change the strand for the 3' UTR bed file to + and -
  awk 'BEGIN{OFS="\t"} {if($6=="1") $6="+"; else if($6=="-1") $6="-"; print}' \
  "${output_dir}"/3_primes_genomic_sorted.bed \
    > "${output_dir}"/3_primes_genomic_sorted_strand.bed

  # subtract the CDS coordinates from the 3' UTRs
  bedtools subtract -s -a "${output_dir}"/3_primes_genomic_sorted_strand.bed \
  -b "${output_dir}"/CDS_coordinates_genomic_all_sorted.bed \
  > "${output_dir}"/3_primes_genomic_CDS_filtered.bed

  # subtract the unique regions from the 3' UTRs
  bedtools subtract -s -a "${output_dir}"/3_primes_genomic_CDS_filtered.bed \
  -b  "${unique_region_dir}"/Unique_DNA_Regions_genomic.bed \
  > "${output_dir}"/3_primes_genomic_CDS_"${region_type}"_filtered.bed


  sort -k1,1 -k2,2n "${output_dir}"/3_primes_genomic_CDS_"${region_type}"_filtered.bed \
  > "${output_dir}"/3_primes_genomic_CDS_"${region_type}"_filtered_sorted.bed

  # merge the overlapping 3'UTR regions to not have duplicates: when smapling the same region can only be selected once
  bedtools merge -i  "${output_dir}"/3_primes_genomic_CDS_"${region_type}"_filtered_sorted.bed -s -c 4,5,6 -o collapse,min,distinct \
  > "${output_dir}"/3_primes_genomic_merged.bed

  python "${script_path}"/region_handling/number_genomic_regions.py \
  "${output_dir}"/3_primes_genomic_merged.bed \
  "${output_dir}"/3_primes_genomic_merged_numbered.bed

 fi