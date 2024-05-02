# dWIZGenomics

The goal of dWIZGenomics is to provide the code used to reproduce the
figures containing the genomic data for the publication "A molecular glue 
degrader of the WIZ transcription factor for fetal hemoglobin induction"

## Instructions

The figures can be reproduced by running  the code in the vignettes. Install
the package by following the next steps:

1. Clone the repository

```
git clone https://github.com/Novartis/dWIZGenomics.git
```

2. Download the data files available [in the zenodo repository](https://doi.org/10.5281/zenodo.11085537) into the
directory `dWIZGenomics/inst/zenodoData`

3. Download the bigwig files from the [GSE247096](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247096) entry into the directory
`dWIZGenomics/inst/bigwigs/GSE247096`.

4. Load the package via `devtools::load_all()` to avoid copying of large 
files (bigwigs) into your `.libPaths()`

5. Run the code in the vignettes.
