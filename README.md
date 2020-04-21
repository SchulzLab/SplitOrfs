# SplitOrfs
A workflow to computationally detect transcripts that could generate several distinct ORFs that match to the same protein. This can be used to find for example nonsense-mediated decay transcripts that make more than one ORF.

## Dependencies

1) *NCBI blast installation*. The workflow uses a local NCBI Blast installation. It uses the makedb and blastp executables. The latest version can be downloaded [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
2) *bedtools2* (version v2.27.1 or higher). Only needed if the predicted protein regions should be overlapped with functional annotations.
3) *Python 2.7* All python scripts are written for Python 2.7.

Both blast and bedtools should be added to the unix PATH variable.


## Usage 

The runSplitOrfs.sh is the masterscript that executes the workflow for generation of transcripts, alignment to proteins as well as overlap with functional annotation if desired. The script calls a number of python scripts or executable as part of the workflow. Exact information can be found in runSplitOrfs.sh. 

As input you need to give two fasta files. 

1.One is a Fasta file of the proteins that you want to use for the analysis. It has the following format for the header: 

**>ENSG00000001626|ENST00000003084**

Where the first entry is the gene identifier and the second one the protein ID/transcript ID from which the protein was made. They need to be separated by a | character

2.The second a Fasta file that contains transcript sequences. From these transcript sequences possible open reading frames (ORFs) are generated and translated to peptides and then used for alignment with BlastP against the proteins.

The format looks similar to above:

**>ENSG00000100767|ENST00000216658**

The only difference is that the second ID is the transcript ID that represents the transcript to be subjected to ORF generation.

If desired an additional **annotation bed file** can be supplied for checking with overlap of known protein domains. It has the following format (header only shown for illustration should not be in the file):

|Protein/Transcript ID   | Start   | End   | Identifier  |
|---|---|---|---|
|ENST00000308027| 21  |    274  |   PF07690|
|ENST00000574588|104 |    414   |  PF00038|

The first column denotes the Protein or Transcript ID representing the protein (here the Ensembl human transcipt ID of the protein). The second and third denote the start and end of the domain annotation in the protein. The last column is the identifier of the domain type (here PFAM domain). The annotation of human and mouse proteins can be found in the folder *annotations* in the repo.

The default function call is:

```javascript
runSplitOrfs.sh outputFolder proteins.fa transcripts.fa annotation.bed

```

If you want to execute the analysis without computing overlap to annotations use:


```javascript
runSplitOrfs.sh outputFolder proteins.fa transcripts.fa

```

## Output
The pipeline produces a number of files as output, some are just intermediates not of relevance. The relevant ones are:

1. OrfProteins.fa -  a fasta file of all generated proteins (in all three reading frames) from the transcripts in the transcripts.fa file supplied to runSplitOrfs.sh.
2. ProteinDatabase.phr, ProteinDatabase.psq, ProteinDatabase.pin - database files from Blast
3. UniqueProteinORFPairs.txt - the final set of transcripts, that have at least 2 ORFs matching to one of the proteins supplied in proteins.fa. Format explained below.
4. UniqueProteinORFPairs_annotated.txt - an extended file from above, when you also ad an annotation bed file to the pipeline.

The different columns of the UniqueProteinORFPairs_annotated.txt file are explained below.

|column name|explanation|
|---|---|
|**geneID** | Gene identifier|
|**targetTransID** | Identifier of the target protein/transcript that the ORFs have been aligned to.|
|**OrfTransID** | Identifier from which the ORFs have been generated.|
|**NumOrfs** | Number of ORFs matching to targetTransID.|
|**OrfIDs** | The unique ORF identifiers representing the ORFs that aligned (comma separated). These ORFs can be found in the file OrfProteins.fa in the output folder.|
|**OrfPos**| The nucleotide start-stop positions from which the ORF was generated (comma separated list for all matching ORFs).|
|**OrfLengths** | The nucleotide length of the matching ORFs (comma separated).|
|**OrfSeqIdents** | The sequence identity values of the ORF-protein alignments as reported by BlastP (comma separated).|
|**MinSeqIdent** | Minimal observed sequence identity of all the ORF-protein matches.|
|**MaxSeqIdent** | Maximal observed sequence identity of all the ORF-protein matches.|
|**protAlignPos** | Alignment start-stop positions of the ORF in the protein.|
|**ProtCoverage** | Number of amino acid positions covered of the original protein by alignment from all ORFs.|
|**ORF-DomainAnnot** | Identifiers of annotations that overlap with an ORF (comma separated list in order of the ORFs). NA means *not available*, when no ORF annotation existed.|
|**NumOrfAnnot** | Number of ORFs that have at least one overlapping annotation.|
|**AnnotPercent** | The ratio of ORFs that have at least one overlapping annotation.|

The last 3 columns only exist, when an annotation.bed file is supplied to runSplitORFs.sh
