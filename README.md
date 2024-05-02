# dWIZGenomics

The goal of dWIZGenomics is to provide the code used to reproduce the
figures containing the genomic data for the publication "A molecular glue 
degrader of the WIZ transcription factor for fetal hemoglobin induction"

## Installation

You can install the development version of dWIZGenomics like so:

``` r
devtools::install_package("Novartis/dWIZGenomics")
```

## Instructions

The figures can be reproduced by running  the code in the vignettes. Install
the package by following the next steps:

1. Clone the repository

```
git clone <URL to repo>
```

2. Download the data files available the zenodo repository <XX> into the
directory `dWIZGenomics/inst/zenodoData`

3. Download the data files available in the GEO entry <XX> into the directory
`dWIZGenomics/inst/bigwigs/GSE247096`.

4. Load the package via `devtools::load_all()` to avoid copying of large 
files (bigwigs) into your `.libPaths()`

5. Run the code in the vignettes.
