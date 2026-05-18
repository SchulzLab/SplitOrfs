##############################################################################
# BackgroundRegions_adapted_genomic.py takes a bedfile and a fasta as input and
# generates random regions of the fasta file with the length distribution as
# given in the bed file and returns them in bed format
# for the genmic regions the selection of the random regions is made by the
# name, not by the chromsome (the chromsome was the name for transcriptomic
# alignment)
##############################################################################

# usage: python BackgroundRegions_adapted.py in.bed in.fasta out.bed
import sys
import random
from pybedtools import BedTool


def get_length_dist(bed_file):
    # read in bedfile and obtain length distribution
    UR_region_lengths = {}
    for line in bed_file:
        elems = line.split('\t')
        length = int(elems[2])-int(elems[1])
        if elems[3] in UR_region_lengths.keys():
            UR_region_lengths[elems[3]] = UR_region_lengths[elems[3]] + length
        else:
            UR_region_lengths[elems[3]] = length
    lengthdistribution = list(UR_region_lengths.values())
    return lengthdistribution


def get_random_bed_ids(lengthdistribution, UTR_bed_file, outname):
    # select random keys for the fasta file, as many as there are
    # entries in the bed file
    # do not allow for duplicate keys
    # Start time for Step 1
    three_prime_UTRs = BedTool(UTR_bed_file)
    random_regions_lengths = [(interval.name, interval.length, interval.chrom,
                               interval.start, interval.end, interval.strand) for interval in three_prime_UTRs]
    randomlist = []
    # random_choices = []
    with open(outname, 'w') as out:
        for i in range(len(lengthdistribution)):
            required_length = lengthdistribution[i]
            # Ensembl format is 1-based in the genomic coordinate calculation this was not accoutned for
            # need to convert to bed format here
            filtered_intervals = [(name, length, chrom, str(int(start) - 1), stop, strand) for name, length,
                                  chrom, start, stop, strand in random_regions_lengths if length >= required_length + 18]
            random_seq = random.choices(filtered_intervals, k=1)[0]
            random_key = random_seq[0]
            # if key has already been sampled, randomly choose a new key
            assert random_key not in randomlist
            randomlist.append(random_key)

            UTR_length = int(random_seq[4]) - int(random_seq[3])
            assert required_length + 18 <= UTR_length
            if required_length + 18 < UTR_length and random_seq[5] == '+':
                start = random.randint(
                    int(random_seq[3]) + 18, UTR_length - required_length + int(random_seq[3]))
            elif required_length + 18 < UTR_length and random_seq[5] == '-':
                start = random.randint(
                    int(random_seq[3]),  UTR_length - required_length + int(random_seq[3]) - 18)
            elif UTR_length == required_length + 18 and random_seq[5] == '+':
                start = int(random_seq[3]) + 18
            elif UTR_length == required_length + 18 and random_seq[5] == '-':
                start = int(random_seq[3])

            end = start + required_length
            if random_seq[5] == '+':
                out.write(random_seq[2] + '\t' + str(start) + '\t' + str(end) + '\t' + str(random_key)
                          + '\t' + str(0) + '\t' + '+' + '\n')
            elif random_seq[5] == '-':
                out.write(random_seq[2] + '\t' + str(start) + '\t' + str(end) + '\t' + str(random_key)
                          + '\t' + str(0) + '\t' + '-' + '\n')

            # filter the regions iteratively to prevent having the same name twice in the list
            random_regions_lengths = [(name, length, chrom, start, stop, strand) for name, length,
                                      chrom, start, stop, strand in random_regions_lengths if name != random_key]


if len(sys.argv) < 3:
    print('usage python BackgroundRegions_adapted_genomic.py in.bed UTR.bed out.bed')
else:
    random.seed(int(sys.argv[4]))
    file = open(sys.argv[1], 'r')
    lengthdistribution = get_length_dist(file)

    get_random_bed_ids(lengthdistribution, sys.argv[2], sys.argv[3])

    print(f'Random regions generated for iteration {sys.argv[4]}')
