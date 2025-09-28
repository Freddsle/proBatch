# proBatch

## General Overview

The proBatch package facilitates batch effects analysis and correction high-throughput experiments. 
Although the package has primarily been developed for DIA (SWATH) proteomics data, 
it should also be applicable to most omic data with minor adaptations.
    
The package contains functions for diagnostics (proteome/genome-wide and feature-level), 
correction (normalization and batch effects correction) and quality control.

Diagnostics part of the package features unified color scheme for plotting, 
    that allows to produce publication-quality graphs.

Correction functions are convenient wrappers for common normalization and batch 
effects removal approaches such as quantile normalization and median centering. 
Furthermore, the package includes non-linear fitting based approaches to deal 
with complex, mass spectrometry-specific signal drifts.

Quality control step, mostly based on correlation analysis, allows to assess whether 
the correction improved the quality of the data.

All steps of batch effects analysis and correction are illustrated in the vignette,
    using the subset of real-world large-scale dataset.

Please use following manuscript for citation: 
Diagnostics and correction of batch effects in large-scale proteomic studies: a tutorial. Molecular Systems Biology 17, e10240 (2021).


## Installing

To install the latest version of proBatch package, you need `BiocManager`:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("proBatch") 

library(proBatch)
```


NOTE: You might need to also install the following linux packages:
`apt-get install libxml2-dev libz-dev`


## Exploring the package

The complete documentation:
```
help(proBatch)
```

Browse the vignette:
```
browseVignettes('proBatch')
browseVignettes('proBatchFeatures')
```


## Citation
Diagnostics and correction of batch effects in large-scale proteomic studies: a tutorial. Molecular Systems Biology 17, e10240 (2021). https://doi.org/10.15252/msb.202110240
