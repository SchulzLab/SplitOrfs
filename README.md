# SplitOrfs
A workflow to computationally detect transcripts that could generate several distinct ORFs that match to the same protein. This can be used to find for example nonsense-mediated decay transcripts that make more than one ORF.

## Dependencies

1) *NCBI blast installation*. The workflow uses a local NCBI Blast installation. It uses the makedb and blastp executables. The latest version can be downloaded [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
2) *bedtools2* (version v2.27.1 or higher). Only needed if the predicted protein regions should be overlapped with functional annotations.
3) *Python 2.7* All python scripts are written for Python 2.7.

Both blast and bedtools should be added to the unix PATH variable.


## Usage 

The runSplitOrfs.sh is the masterscript that executes the workflow for generation of transcripts, alignment to proteins as well as overlap with functional annotation if desired. The script calls a number of python scripts or executable as part of the workflow. Exact information can be found in runSplitOrfs.sh. The default function call is:

```javascript
runSplitOrfs.sh outputFolder proteins.fa transcripts.fa annotation.bed

```

If you want to execute the analysis without computing overlap to annotations use:


```javascript
runSplitOrfs.sh outputFolder proteins.fa transcripts.fa

```
