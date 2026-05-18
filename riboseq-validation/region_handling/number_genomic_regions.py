import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="."
    )

    # Required positional arguments
    parser.add_argument("input_bed_file",
                        help="Path to input BED file")

    parser.add_argument("output_bed_file",
                        help="Path to output BED file")

    return parser.parse_args()


def main(input_bed_file, output_bed_file):
    """
    Adds a number to the name field in the BED file to circumvent having the same
    name several times for different regions

    Args:
        input_bed_file (str): path to the input BED file that is to be numbered
        output_bed_file (str): path to the numbered output BED file
    """
    input_bed_df = pd.read_csv(input_bed_file, sep='\t', header=None)
    input_bed_df.iloc[:, 3] = input_bed_df.iloc[:, 3] + \
        '|' + list(map(str, range(len(input_bed_df.index))))
    input_bed_df.to_csv(output_bed_file, sep='\t', header=False, index=False)


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    input_bed_file = args.input_bed_file
    output_bed_file = args.output_bed_file

    main(input_bed_file, output_bed_file)
