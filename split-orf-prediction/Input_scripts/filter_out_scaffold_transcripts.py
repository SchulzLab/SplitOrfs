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
        "fasta",
        help="Merged input FASTA file (transcript cDNA sequences).",
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


def read_fasta(fasta_path):
    """Yield (header, sequence) pairs from a FASTA file.

    header is returned without the leading '>' and without the trailing
    newline; sequence has all its lines concatenated.
    """
    header, seq_chunks = None, []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)


def main():
    args = parse_args()

    # Unpack arguments into local variables
    fasta = args.fasta
    reference_gtf = args.reference_gtf
    output = args.output

    # --- Load reference: transcript_id -> chromosome ---
    gn2chrom = load_gene_to_chrom(reference_gtf)

    chr_genes = []

    with open(output, "w") as out:
        for header, seq in read_fasta(fasta):
            gid = header.split('|')[0]
            chrom = gn2chrom.get(gid, "scaffold")
            if chrom != "scaffold":
                out.write(f">{header}\n")
                out.write(seq + "\n")
                chr_genes.append(gid)

    txt_file_name = output.removesuffix(".fa") + "gid.txt"
    with open(txt_file_name, "w") as f:
        f.write("\n".join(chr_genes))

    # fasta = '/home/ckalk/tools/SplitORF_pipeline/Input2023/NMD_transcripts_CDNA.fa'
    # reference_gtf = '/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf'
    # output = '/home/ckalk/tools/SplitORF_pipeline/Input2023/NMD_transcripts_CDNA_no_scaffold.fa'


    # fasta = '/home/ckalk/tools/SplitORF_pipeline/Input2023/RI_transcripts_CDNA.fa'
    # reference_gtf = '/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.110.chr.gtf'
    # output = '/home/ckalk/tools/SplitORF_pipeline/Input2023/RI_transcripts_CDNA_no_scaffold.fa'
if __name__ == "__main__":
    main()
