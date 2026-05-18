# This script filters the unique regions smaller than the Mummer length parameter
# As Mummer identifies exact matches (non-unique regions) only longer than the parameter
# everything is labelled as unique, which might not be meaningfull for such small regions
import sys

# Bedfile to be merged
unique_regions_bed = open(sys.argv[1], "r")
unique_regions = unique_regions_bed.readlines()
# Parameter X ()should be set at the same value as the minimum length parameter of the MUMmer maxmatch call.
min_length = int(sys.argv[2])

# Resulting bedfile
with open(sys.argv[3], "w") as out:
    for unique_region in unique_regions:
        columns = unique_region.split()
        if (int(columns[2])-int(columns[1]) >= min_length):
                    out.write(columns[0] + "\t" + columns[1] + "\t" + columns[2] + "\n")