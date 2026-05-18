import argparse
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Get list of Split-ORF genes wihtout duplicates.")

    parser.add_argument(
        "unique_protein_orf_pairs",
        type=str,
        help="Path to the UniqueProteinORFPairs.txt file from the Split-ORF pipeline"
    )

    parser.add_argument(
        "split_orf_gene_file",
        type=str,
        help="Path to the output file with the Split-ORF genes (txt)"
    )

    return parser.parse_args()


def main(unique_protein_orf_pairs, split_orf_gene_file):
    unique_protein_orf_pairs_df = pd.read_csv(
        unique_protein_orf_pairs, header=0, sep='\t')
    split_orf_genes = list(unique_protein_orf_pairs_df['geneID'].unique())
    with open(split_orf_gene_file, 'w') as f:
        for so in split_orf_genes:
            f.write(f"{so}\n")


if __name__ == "__main__":
    args = parse_arguments()
    main(args.unique_protein_orf_pairs, args.split_orf_gene_file)
