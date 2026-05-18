# ------------------ IMPORTS ------------------ #

import pandas as pd
import argparse
import os
import os.path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Get union of SO validated NMD and RI genes as well as union of background all NMD and RI genes for enrichment analysis DISGENET."
    )

    # Required positional arguments
    parser.add_argument("--nmd_val_genes",
                        help="TXT file of NMD validated genes")
    parser.add_argument("--ri_val_genes",
                        help="TXT file of NMD validated genes")
    parser.add_argument("--nmd_so_genes",
                        help="TXT file of NMD SO genes")
    parser.add_argument("--ri_so_genes",
                        help="TXT file of NMD SO genes")
    parser.add_argument("--nmd_all_genes",
                        help="TXT file of NMD all genes")
    parser.add_argument("--ri_all_genes",
                        help="TXT file of NMD all genes")
    parser.add_argument("--outdir",
                        help="Outdir to put resulting TXT files")

    return parser.parse_args()


def main(nmd_val_genes, ri_val_genes, nmd_so_genes, ri_so_genes, nmd_all_genes, ri_all_genes, outdir):
    def perform_union_intersection_unique_analysis(nmd_genes, ri_genes, gene_type):
        nmd_genes_df = pd.read_csv(nmd_genes, header=None)
        ri_genes_df = pd.read_csv(ri_genes, header=None)

        # union val genes
        union_genes = pd.concat(
            [nmd_genes_df, ri_genes_df], ignore_index=True)
        union_genes = union_genes.drop_duplicates(ignore_index=True)
        union_genes.to_csv(os.path.join(
            outdir, f'Union_RI_NMD_{gene_type}_genes.txt'), index=False, header=False)
        print(f'Number of {gene_type} genes in union NMD and RI',
              len(union_genes.index))

        # intersection val genes
        intersection_genes = set(nmd_genes_df[0]) & set(ri_genes_df[0])
        with open(os.path.join(outdir, f"Intersection_RI_NMD_{gene_type}_genes.txt"), "w") as file:
            for item in intersection_genes:
                file.write(f"{item}\n")
        print(f'Number of {gene_type} genes in intersection NMD and RI',
              len(intersection_genes))

        # unique NMD and RI genes
        unique_nmd_genes = set(gene for gene in set(
            nmd_genes_df[0]) if gene not in set(ri_genes_df[0]))
        print(f'Number of {gene_type} genes uniquely in NMD',
              len(unique_nmd_genes))
        with open(os.path.join(outdir, f"Unique_NMD_{gene_type}_genes.txt"), "w") as file:
            for item in unique_nmd_genes:
                file.write(f"{item}\n")

        unique_ri_genes = set(gene for gene in set(
            ri_genes_df[0]) if gene not in set(nmd_genes_df[0]))
        print(f'Number of {gene_type} genes uniquely in RI',
              len(unique_ri_genes))
        with open(os.path.join(outdir, f"Unique_RI_{gene_type}_genes.txt"), "w") as file:
            for item in unique_ri_genes:
                file.write(f"{item}\n")

    ####### validated genes ####################################################
    perform_union_intersection_unique_analysis(
        nmd_val_genes, ri_val_genes, 'validated')

    ####### Split-ORF genes ####################################################
    perform_union_intersection_unique_analysis(
        nmd_so_genes, ri_so_genes, 'Split-ORF')

    ####### All genes ##########################################################
    perform_union_intersection_unique_analysis(
        nmd_all_genes, ri_all_genes, 'all')


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    nmd_val_genes = args.nmd_val_genes
    ri_val_genes = args.ri_val_genes
    nmd_so_genes = args.nmd_so_genes
    ri_so_genes = args.ri_so_genes
    nmd_all_genes = args.nmd_all_genes
    ri_all_genes = args.ri_all_genes
    outdir = args.outdir

    nmd_val_genes = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/NMD/SO_validated_genes_NMD.txt"
    ri_val_genes = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/RI/SO_validated_genes_RI.txt"
    nmd_all_genes = "/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-13.07.17_NMD_for_paper/BackgroundGeneFile.txt"
    ri_all_genes = "/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-14.09.13_RI_for_paper/BackgroundGeneFile.txt"
    outdir = "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis/SO_valdiation/union_intersection"

    main(nmd_val_genes, ri_val_genes, nmd_so_genes,
         ri_so_genes, nmd_all_genes, ri_all_genes, outdir)
