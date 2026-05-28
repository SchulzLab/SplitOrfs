# Split-ORF pipeline

Split Open Reading frames (Split-ORFs) exist on transcripts that encode two or more open reading frames.
These multiple open reading frames each encode a part of a full length protein, splitting it into pieces. 

The Split-ORF pipeline predicts possible Split-ORF transcripts in a user defined set of transcript sequences.
Additionally, DNA and protein unique regions are predicted for the Split-ORFs. These regions do not occur in any
other annotated protein coding transcript. 

These unique regions can then, in a second step, be tested for significantly more Ribo-seq coverage than background regions.

## Conda package
We are currently working on a bioconda package for the Split-ORF pipeline. The conda package can already be installed following the instructions below and will be made available via bioconda for easier installation within the next weeks (date: 18.05.2026).

## Installation 

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

# split-orf-prediction

## Input Files split-orf-prediction
The input files are defined in a JSON file which is the only argument of the **split-orf-prediction** module. All of the input files should be in the same folder. Example JSON files of the inputs are supplied in the split-orf-prediction folder: e.g. split_orf_pipeline_input.json.

The following arguments need to be set within the JSON file:<br>
```json
{
    "file_path": path where the input files are located
    (all need to be in the same folder),
    "proteins": Reference protein sequences (amino acid sequences, FASTA) (1),
    "transcripts": Transcripts sequences for Split-ORF prediction (FASTA) (2),
    "annotation": PFAM annotation (with columns: transcript ID, 
    start position annotation, stop position annotation, PFAM ID, TSV) (3),
    "reference_transcripts": Reference protein coding transcript sequences
     (DNA sequences, FASTA) (4),
    "exon_positions": Exon coordinates (can be downloaded from Ensembl Biomart
     Strucutres with the following columns:
    Gene stable ID,	Transcript stable ID,	Exon region start (bp),	
    Exon region end (bp),	Transcript start (bp),
  	Transcript end (bp),	Strand,	Chromosome/scaffold name, TSV) (5),
    "align_method": Alignment algorithm (blast or diamond),
    "cds_coordinate_bed": CDS and contamination coordinates to subtract 
    from Split-ORF unique regions (Chromosome/scaffold name, Gene stable ID,  
    Transcript stable ID, Genomic coding start, Genomic coding end, Strand) (6),
    "output_dir": directory where outout files should be written to
}
```

### How to get Input files
All of the Input and Output files which were used for the Split-ORF prediction in the pipeline paper (add link), can be downloaded via Zenodo (add link). This is the preferred option if looking for Split-ORFs in NMD or RI transcripts. If Split-ORFs should be predicted in other transcripts, follow the description of the Input files closely and make sure the expected format is met.


All of these input files were downloaded from Ensembl Biomart (Ensembl Genes 110, GRCh38p.14) using the following mart options and transformed with the respective custom scripts which can be found in the Input_scripts folder of the Split-Orfs github repository: 

1. Attributes: Sequences: Peptide, header information: Gene stable ID, Transcript stable ID, Filters: Gene: Transcript type: protein_coding.
2. Attributes: Sequences: cDNA sequences, header information: Gene stable ID, Transcript stable ID, Filters: Gene: Transcript type: nonsense_mediated_decay or retained_intron. <
3. Attributes: Features: Gene: Gene stable ID Transcript stable ID; Protein Domains and Families: Pfam ID, Pfam start, Pfam end. Remodelling scripts: convert_ensembl_output_to_bed.py.
4. Attributes: Sequences: cDNA sequences, header information: Gene stable ID, Transcript stable ID, Filters: Gene: Transcript type: protein_coding.
5. Attributes: Structures: Gene stable ID, Transcript stable ID, Exon region start, Exon region end, Transcript start, Transcript End, Strand. No filters were selected. Downloading all exonic positions includes those of the transcripts of interest.
6. for the contaminations: Attributes: Structures: Chromosome/scaffold name, Gene stable ID, Transcript stable ID, Exon region start, Exon region end, Strand. Filter: Gene: the respective contaminating RNA (snoRNA, miRNA, snRNA, scaRNA, rRNA, rRNA_psuedogene, Mt_tRNA, Mt_rRNA).
for the CDS coordinates Attributes: Structures: Chromosome/scaffold name, Gene stable ID,  Transcript stable ID, Genomic coding start, Genomic coding end, Strand.


Protein coding transcript (4) and protein sequences (1) as well as CDS coordinates (6) were filtered for transcript support level 1 or 2 or the presence of their exact intron chain in the RefSeq annotation (v.GCF_000001405.40-RS_2023_10) with custom scripts (filter_Ensembl_GTF.sh, Filter_prot_coding_reference.sh). The filtered CDS coordinates were combined with the contamination coordiantes into a single BED file using the remodelling script: generate_CDS_contamination_subtraction_coordinates.sh.<br>

The Split-ORF pipeline creates an output folder with a timestamp of the run at a user specified location. All results as well as intermediate result files are written into this output directory. <br>
The final output files are a TSV file of the predicted Split-ORFs, BED files of the genomic coordinates of the unique Split-ORF regions and two HTML reports about the predicted Split-ORF candidates and about their unique regions. The steps of the Split-ORF pipeline produce intermediate results which are also included in the output of the pipeline, but these files are not relevant as a result.



### File structure

The FASTA files (1,2,4) used for the Split-ORF pipeline need to have the following format for the header: 

**>ENSG00000001626|ENST00000003084**

Where the first entry is the gene identifier and the second one the gene ID/transcript ID from which the protein was made. They need to be separated by a | character


The **annotation bed file** (3) can be supplied for checking with overlap of known protein domains. It has the following format (header only shown for illustration should not be in the file):

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
5. so_categorization_df.csv - CSV file of all predicted Split-ORFs giving information on the transcript level about the number of unique regions, whether these are located in the first or second ORF and whether unqiue regions overlap within the transcript
6. UniqueProteinORFPairs.txt - the final set of transcripts, that have at least 2 ORFs matching to one of the proteins supplied in proteins.fa. Format explained below.
7. UniqueProteinORFPairs_annotated.txt - an extended file from above, when you also ad an annotation bed file to the pipeline.

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


# ribo-cov
The **ribo-cov** module performs a statistical test with BH multiple testing correction of the unique regions of the predicted Split-ORFs of the **split-orf-prediction** step. **Unique_DNA_Regions_genomic_final.bed** , the BED file with the genomic positions of the unique regions is required in order to run the **ribo-cov** module, as well as preprocessed FASTQ or BAM files of Ribo-seq data and the genomic annotation files. The annotation version used for the **split-orf-prediction** is required to match the annotation version of the **ribo-cov** module.

## Input Files ribo-cov
Example Input files are supplied in the riboseq-validation folder, e.g. riboseq_validation_input_server_NMD.json.
```json
{
    "output_star": Directory where the Ribo-seq BAM files of the STAR mapping step are stored or should be stored,
    "unique_region_dir": Directory of the split-orf-prediction Output containing Unique_DNA_Regions_genomic_final.bed,
    "ensembl_gtf": Ensembl GTF file to use for mapping,
    "genome_fasta": genome FASTA file,
    "input_fastq_path": Directory containing FASTQ files of the Ribo-seq data 
    (only indicate if data not yet mapped),
    "three_primes": Three prime UTR genomic coordinates to use as the background in BED format,
    "cds_coordinates": CDS coordinate BED file,
    "input_name": name to give the BAM files (NMD),
    "region_type": region type (NMD or RI),
    "bam_ending": name of the bam file endings, when FASTQ files are used should be _input_name_sorted.bam
    otherwise the ending of the supplied BAM files e.g. "deduplicated.bam",
    "tmp_dir": TMP directory for intersection of large BED files,
    "report": Path to PDF file which will be the PDF report of the ribo-cov module,
    "duplicated": boolean indicating whether reads are duplicated or whether UMIs were used 
    and reads are already deduplicated
}
```
The Ribo-seq data can make use of UMIs or not. If UMI containing Ribo-seq data is used, it is required to be already mapped to the genome (using the correct version!) and deduplicated. If the deduplicated reads are used then the "duplicated" flag should be set to false, otherwise it should be set to false. The ribo-cov module is only tested using STAR as the aligner, so it is recommended to use it for the alignment. Please note that if supplying BAM files, these need to be placed in a directory file_dir=output_star/input_name_genome/subdirectory". The subdirectory name can be chosen freely. This is required as the when the mapping is performed, the BAM files will be placed in subdirectories of the folder output_star/input_name_genome.

The results will be placed in a folder: output_star/region_type_genome.

The distinction between region_type and input_name is done in order to enable the reuse of BAM files, if the same Ribo-seq data is used for different runs of the Split-ORF pipeline. In our example we use it for the NMD and the RI transcripts. We first run it on the NMD transcripts and then reuse the BAM files to circumvent redundant mapping steps by setting input_name to NMD when the region_type is RI for the RI run and the BAM ending to _NMD_sorted.bam. The output_star directory is set to the parent directory of both NMD_genome and RI_genome, which are created by the ribo-cov module.

The **three_primes** should contain the genomic coordinates of 3' UTR regions that will be used for the background Ribo-seq coverage used for the statistical test. We decided to only consider 3' UTRs of protein coding transcripts filtered for transcript support level 1 or 2 or the presence of their exact intron chain in the RefSeq annotation (v.GCF_000001405.40-RS_2023_10). The coordinates were downloaded from Ensembl Biomart (v.110) and then filtered with the get_3prime_genomic_coords.sh in the region_handling directory (within the riboseq-validation directory). It is also possible to supply a different set of regions to use for the background, but we recommend to make sure that these only contain regions which actually are untranslated, i.e. 5' UTRs may have Ribo-seq coverge due to the scanning of the ribosome. We also subtract regions of the 3' UTRs that overlap with any CDS as annotated in Ensembl. It may happen that the 3' UTR of one transcript isoform overlaps with the CDS of another transcript isoform of the same gene. These regions should be considered as CDSs in terms of the Ribo-seq coverage and not included as regions for the background.

The cds_coordinates are downloaded from Ensembl Biomart (v.110) and are the same as the ones used for the split-orf-prediction module, but without the contaminations.

All Input and Output files can be downloaded from [Zenodo](https://zenodo.org/records/20340925).

## Output Files ribo-cov

A PDF report of the unique regions with significantly more Ribo-seq coverage than in the background regions with Bh multiple testing correction is generated automatically by the pipeline. The report shows the number of unique regions validated for each sample as well as the overlap of the validated unique regions among samples using Upsetplots.

Two CSV files with the endings of **genomicunique_regions.csv** and **unique_regions.csv** are generated per sample and contain the information about all unique regions with significant Ribo-seq coverage, such as the number of reads mapping (num_reads), the p- and q-value and the name of the transcript, ORF and exact position of the respective unqiue region (new_name).

In the Split-ORF pipeline paper we performed in depth analysis of the ribo-cov unique regions with the scripts placed within the riboseq-validation/SO_categorization_by_UR_coverage_results folder. This analysis resulted in the sunburst plots as shown in the paper and the **so_categorization_df.csv** file from the split-prf-prediction step together with the **unique_regions.csv** result files from the Ribo-seq pipeline.

## Contact
Corresponding author: kalk@med.uni-frankfurt.de

## Citation
Please cite our manuscripts on [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.05.22.727176v1).


## License

