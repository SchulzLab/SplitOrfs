import sys
uniquefile = open(sys.argv[1], 'r')
uniquelist = []
for unique_region in uniquefile:
    unique_info = unique_region.strip()
    unique_info = unique_info.split("\t")
    unique_ids = unique_info[0].split(":")
    temp = [unique_ids[0], unique_info[1], unique_info[2],
            unique_ids[1], unique_ids[2], unique_ids[3]]
    uniquelist.append(temp)
with open(sys.argv[2], 'w') as f:
    for i in uniquelist:
        # if the protein ends at the end of the ORF, correct for the termination codon positions
        if 3 * int(i[2]) + int(i[4]) == int(i[5]) - 3:
            # give the unique regions in terms of transcript not of the orf, needed for intersection with riboseq data
            f.write(i[0] + "\t" + str(3 * int(i[1]) + int(i[4])) + "\t" + str(3 *
                    int(i[2]) + int(i[4]) + 3) + "\t" + i[3] + ":" + i[4] + ":" + i[5] + "\n")
        else:
            f.write(i[0] + "\t" + str(3 * int(i[1]) + int(i[4])) + "\t" + str(3 *
                    int(i[2]) + int(i[4])) + "\t" + i[3] + ":" + i[4] + ":" + i[5] + "\n")
