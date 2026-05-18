import sys
file = open(sys.argv[1], 'r')
with open(sys.argv[2], 'w') as f:
    previous_tid = None
    for exon_region in file:
        exon_region = exon_region.strip()
        exon_info = exon_region.split("\t")
        strand = exon_info[6]
        tid = exon_info[1]
        chromosome = exon_info[7]
        # to transform the coordinates into the bed format need to substract one
        # from the start of the Ensembl 1-based coordinates
        exon_start_coord = int(exon_info[2]) - 1
        exon_end_coord = int(exon_info[3])
        whole_transcript_start_coordinate = int(exon_info[4]) - 1
        whole_transcript_end_coordinate = int(exon_info[5])
        if strand == '1' or strand == '+':
            if tid != previous_tid:
                assert (exon_start_coord == whole_transcript_start_coordinate)
                transcript_start_coordinate = 0
                # in bed format there is no need to add plus one to the length
                # there is somehow, as I checked in the Ensembl webbrowser: e.g. exons: ENST00000482929
                # the calculation of end - start is always one too little than the reported length
                transcript_end_coordinate = exon_end_coord - exon_start_coord
                previous_tid = tid
            else:
                # bed files are half open this means that the end coordinate is not included!
                # Therefore, the previous end coordinate is the current start :)
                transcript_start_coordinate = previous_end_coordinate
                transcript_end_coordinate = previous_end_coordinate + \
                    exon_end_coord - exon_start_coord
        elif strand == '-1' or strand == '-':
            if tid != previous_tid:
                assert (exon_end_coord == whole_transcript_end_coordinate)
                transcript_start_coordinate = 0
                transcript_end_coordinate = exon_end_coord - exon_start_coord
                previous_tid = tid
            else:
                transcript_start_coordinate = previous_end_coordinate
                transcript_end_coordinate = previous_end_coordinate + \
                    exon_end_coord - exon_start_coord
        f.write(exon_info[0] + "|" + tid + "\t" + str(exon_start_coord) + "\t" + str(exon_end_coord) +
                '\t' + str(transcript_start_coordinate) + '\t' + str(transcript_end_coordinate) + '\t' + strand + '\t' + chromosome + "\n")
        previous_end_coordinate = transcript_end_coordinate
