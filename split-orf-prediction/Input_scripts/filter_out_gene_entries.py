#filter gene feature entries out of the gtf file
import sys
from typing import List
from pygtftk.gtf_interface import GTF


def main():


    databank_gtf = GTF(sys.argv[1], check_ensembl_format=False)
    genes = databank_gtf.get_gn_ids(nr= True)
    transcripts = databank_gtf.get_tx_ids(nr= True)
    filtered_gtf = databank_gtf.select_by_key('feature', 'gene', invert_match = True)
    
    print('Number of transcripts in normal gtf:', len(transcripts))
    print("Number of genes in normal gtf:", len(genes))

    filtered_gtf.write(sys.argv[2])
   
if __name__ == "__main__":
    main()




   