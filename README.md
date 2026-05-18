# Split-ORF pipeline

Split Open Reading frames (Split-ORFs) exist on transcripts that encode two or more open reading frames.
These multiple open reading frames each encode a part of a full length protein, splitting it into pieces. 

The Split-ORF pipeline predicts possible Split-ORF transcripts in a user defined set of transcript sequences.
Additionally, DNA and protein unique regions are predicted for the Split-ORFs. These regions do not occur in any
other annotated protein coding transcript. 

## Conda package
We are currently working on a bioconda package for the Split-ORF pipeline. The package can be installed floowing the instructions below and will be made available via bioconda for easier installation within the next weeks (state: 18.05.2026).

## Usage 

The Split-ORF pipeline has two modules: **split-orf-prediction** and **ribo-cov**. The **split-orf-prediction** module allows the prediction of Split-ORFs for a user supplied set of transcripts together with their unique regions, while the **ribo-cov** module performs testing for significant Ribo-seq coverage in the Split-ORF unique regions of the split-orf-prediction module based on user supplied Ribo-seq data.<br>
The Split-ORF pipeline is implemented in a modular fashion within the bash framework. The different modules of the pipeline are written in python (v.3.11) and the reports are generated with R (v.4.1.2) with the Rmarkdown framework.

To install the Split-ORF pipeline conda package run the following commands:
```bash
git clone git@github.com:SchulzLab/SplitOrfs.git
cd SplitOrfs/conda-recipe
conda/mamba create -n splitorf_env -c local -c bioconda -c conda-forge
conda activate splitorf_env
conda build .
```
When the respective conda environment is acitvated, then the **split-orf-prediction** and **ribo-cov** can be run from the command line with the respective JSON files as input that are described in more detail below.


## Input Files split-orf-prediction
The input files are defined in JSON files and should all be in the same folder. Example JSON files of the inputs are supplied in the split-orf-prediction folder: e.g. split_orf_pipeline_input.json.

The following arguments need to be set within the JSON file:<br>
```json
{
    "file_path": path where the input files are located (all need to be in the same folder) (1),
    "proteins": Reference protein sequences (amino acid sequences, FASTA) (2),
    "transcripts": Transcripts sequences for Split-ORF prediction (FASTA) (3),
    "annotation": PFAM annotation (with columns: transcript ID, start position annotation, stop position annotation, PFAM ID, TSV) (4),
    "reference_transcripts": Reference protein coding transcript sequences (DNA sequences, FASTA) (5),
    "exon_positions": Exon coordinates (can be downloaded from Ensembl Biomart Strucutres with the following columns:
    Gene stable ID,	Transcript stable ID,	Exon region start (bp),	Exon region end (bp),	Transcript start (bp),
  	Transcript end (bp),	Strand,	Chromosome/scaffold name, TSV) (6),
    "align_method": Alignment algorithm (blast or diamond) (7),
    "cds_coordinate_bed": CDS and contamination coordinates to subtract from Split-ORF unique regions (Chromosome/scaffold name, Gene stable ID,  Transcript stable ID, Genomic coding start, Genomic coding end, Strand) (8),
    "output_dir": directory where outout files should be written to (9)
}
```

### How to get Input files
All of these input files were downloaded from Ensembl Biomart (Ensembl Genes 110, GRCh38p.14) using the following mart options and transformed with the respective custom scripts which can be found in the Input_scripts folder of the Split-Orfs github repository: <br>
(2) Attributes: Sequences: Peptide, header information: Gene stable ID, Transcript stable ID, Filters: Gene: Transcript type: protein_coding.<br>
(3) Attributes: Sequences: cDNA sequences, header information: Gene stable ID, Transcript stable ID, Filters: Gene: Transcript type: nonsense_mediated_decay or retained_intron. <br>
(4) Attributes: Features: Gene: Gene stable ID Transcript stable ID; Protein Domains and Families: Pfam ID, Pfam start, Pfam end. Remodelling scripts: convert_ensembl_output_to_bed.py.<br>
(5) Attributes: Sequences: cDNA sequences, header information: Gene stable ID, Transcript stable ID, Filters: Gene: Transcript type: protein_coding.<br>
(6) Attributes: Structures: Gene stable ID, Transcript stable ID, Exon region start, Exon region end, Transcript start, Transcript End, Strand. No filters were selected. Downloading all exonic positions includes those of the transcripts of interest.<br>
(8) for the contaminations: Attributes: Structures: Chromosome/scaffold name, Gene stable ID, Transcript stable ID, Exon region start, Exon region end, Strand. Filter: Gene: the respective contaminating RNA (snoRNA, miRNA, snRNA, scaRNA, rRNA, rRNA_psuedogene, Mt_tRNA, Mt_rRNA).<br>
for the CDS coordinates Attributes: Structures: Chromosome/scaffold name, Gene stable ID,  Transcript stable ID, Genomic coding start, Genomic coding end, Strand.  <br>
Protein coding transcript (5) and protein sequences (2) as well as CDS coordinates (8) were filtered for transcript support level 1 or 2 or the presence of their exact intron chain in the RefSeq annotation (v.GCF_000001405.40-RS_2023_10) with custom scripts (filter_Ensembl_GTF.sh, Filter_prot_coding_reference.sh). The filtered CDS coordinates were combined with the contamination coordiantes into a single BED file using the remodelling script: generate_CDS_contamination_subtraction_coordinates.sh.<br>
The input files are specified via a JSON file and the Input data should be located in the same directory. <br>
The Split-ORF pipeline creates an output folder with a timestamp of the run at a user specified location. All results as well as intermediate result files are written into this output directory. <br>
The final output files are a TSV file of the predicted Split-ORFs, BED files of the genomic coordinates of the unique Split-ORF regions and two HTML reports about the predicted Split-ORF candidates and about their unique regions. The steps of the Split-ORF pipeline produce intermediate results which are also included in the output of the pipeline.



### File structure

The FASTA files (2,3,5) used for the Split-ORF pipeline need to have the following format for the header: 

**>ENSG00000001626|ENST00000003084**

Where the first entry is the gene identifier and the second one the gene ID/transcript ID from which the protein was made. They need to be separated by a | character


The **annotation bed file** (4) can be supplied for checking with overlap of known protein domains. It has the following format (header only shown for illustration should not be in the file):

| Protein/Transcript ID | Start | End | Identifier |
| --------------------- | ----- | --- | ---------- |
| ENST00000308027       | 21    | 274 | PF07690    |
| ENST00000574588       | 104   | 414 | PF00038    |

The first column denotes the Protein or Transcript ID representing the protein (here the Ensembl human transcipt ID of the protein). The second and third denote the start and end of the domain annotation in the protein. The last column is the identifier of the domain type (here PFAM domain). The annotation of human and mouse proteins can be found in the folder *annotations* in the repo.


## Output
The pipeline produces a number of files as output, some are just intermediates not of relevance. The relevant ones are:

1. Split-ORF_Report.html - HTML report with summary statistics and plots about the Split-ORF predictions
2. Uniqueness_Report.html - HTML report with summary statistics and plots about the unique regions of the predicted Split-ORFs
3. **Unique_DNA_Regions_genomic_final.bed** - Unique region genomic coordinates to use for the ribo-cov module
4. OrfProteins.fa -  a fasta file of all generated proteins (in all three reading frames) from the transcripts in the transcripts.fa file supplied to runSplitOrfs.sh.
5. UniqueProteinORFPairs.txt - the final set of transcripts, that have at least 2 ORFs matching to one of the proteins supplied in proteins.fa. Format explained below.
6. UniqueProteinORFPairs_annotated.txt - an extended file from above, when you also ad an annotation bed file to the pipeline.

The different columns of the UniqueProteinORFPairs_annotated.txt file are explained below.

| column name         | explanation                                                                                                                                                |
| ------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **geneID**          | Gene identifier                                                                                                                                            |
| **targetTransID**   | Identifier of the target protein/transcript that the ORFs have been aligned to.                                                                            |
| **OrfTransID**      | Identifier from which the ORFs have been generated.                                                                                                        |
| **NumOrfs**         | Number of ORFs matching to targetTransID.                                                                                                                  |
| **OrfIDs**          | The unique ORF identifiers representing the ORFs that aligned (comma separated). These ORFs can be found in the file OrfProteins.fa in the output folder.  |
| **OrfPos**          | The nucleotide start-stop positions from which the ORF was generated (comma separated list for all matching ORFs).                                         |
| **OrfLengths**      | The nucleotide length of the matching ORFs (comma separated).                                                                                              |
| **OrfSeqIdents**    | The sequence identity values of the ORF-protein alignments as reported by BlastP (comma separated).                                                        |
| **MinSeqIdent**     | Minimal observed sequence identity of all the ORF-protein matches.                                                                                         |
| **MaxSeqIdent**     | Maximal observed sequence identity of all the ORF-protein matches.                                                                                         |
| **protAlignPos**    | Alignment start-stop positions of the ORF in the protein.                                                                                                  |
| **ProtCoverage**    | Number of amino acid positions covered of the original protein by alignment from all ORFs.                                                                 |
| **ORF-DomainAnnot** | Identifiers of annotations that overlap with an ORF (comma separated list in order of the ORFs). NA means *not available*, when no ORF annotation existed. |
| **NumOrfAnnot**     | Number of ORFs that have at least one overlapping annotation.                                                                                              |
| **AnnotPercent**    | The ratio of ORFs that have at least one overlapping annotation.                                                                                           |


