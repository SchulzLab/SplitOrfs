import os
import os.path
import argparse

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "."
        )
    )

    # Required positional arguments
    parser.add_argument(
        "background_genes",
        help="TXT file of background genes for GO enrichment.",
    )
    parser.add_argument(
        "reference_gtf",
        help="Reference GTF used to map each transcript to its chromosome.",
    )

    # Optional output path
    parser.add_argument(
        "-o", "--output",
        default="parnet_input.fa",
        help="Path for the output FASTA.",
    )

    return parser.parse_args()


def load_gene_to_chrom(gtf_path):
    """Parse a reference GTF and return {gene_id: chromosome}.

    Reads the GTF, keeps 'gene' feature rows, and extracts the
    gene_id from the attribute column (col 9) paired with the
    chromosome (col 1).
    """
    gn2chrom = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom = fields[0]
            attributes = fields[8]
            # attribute format: gene_id "ENSG..."; ...
            gn_id = None
            for attr in attributes.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_id"):
                    gn_id = attr.split('"')[1]
                    break
            if gn_id is not None:
                gn2chrom[gn_id] = chrom
    return gn2chrom


def main():
    args = parse_args()

    # Unpack arguments into local variables
    background_genes = args.background_genes
    reference_gtf = args.reference_gtf
    output = args.output

    # --- Load reference: transcript_id -> chromosome ---
    gn2chrom = load_gene_to_chrom(reference_gtf)

    bkg_gene_df = pd.read_csv(background_genes, header=None, names=['geneID'])
    bkg_gene_df['chr'] = bkg_gene_df['geneID'].map(gn2chrom).fillna('scaffold')
    bkg_gene_df = bkg_gene_df[bkg_gene_df['chr'] != 'scaffold']
    bkg_gene_df['geneID'].to_csv(output, index=False, header=False)
    # background_genes = '/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/BackgroundGeneFile.txt'
    # reference_gtf = '/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf'
    # output = '/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/BackgroundGeneFile_no_scaffold.txt'

    # background_genes = '/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/BackgroundGeneFile.txt'
    # reference_gtf = '/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf'
    # output = '/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.10.51_RI_contamination_subtraction/BackgroundGeneFile_no_scaffold.txt'
if __name__ == "__main__":
    main()
