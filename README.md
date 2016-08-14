# LTRpred

### Perform _de novo_ Prediction, Annotation, Meta-Genomics, and Analytics of LTR retrotransposons with R


> Transposable genetic elements (TEs) comprise a vast array of DNA sequences, all having the ability to move to new sites in genomes either directly by a cut-and-paste mechanism (transposons) or indirectly through an RNA intermediate (retrotransposons). 
>
>  \- Nina V. Federoff, Science 2012 

Due to their enormous contribution to genome structure and genome evolution transposable elements allow us to study fundamental mechanisms of phenotypic adaptation, diversification,
and evolution. In particular, understanding the recognition and 
regulation of transposable elements by the genetic regulatory machinery will enable us to 
systematically identify the key players and key processes that enable niche adaptation and
species diversification on the genetic level.

The `LTRpred` package aims to provide an integrated software framework to 
predict LTR transposons _de novo_ in any genomic sequence of interest.
LTR transposons have the capacity to move to new sites in genomes
through a copy-and-paste mechanism and by doing so are able to contribute generatively 
to genome evolution and environmental sensing on the genetic level.
Hence, predicting the presence of LTR transposons within genomes as well as their
capacity to perform this copy-and-paste strategy enables us to quantify the extent 
to which transposons shape the adaptation and evolution of life in general.

In particular the following analyses can be performed with `LTRpred`:

### _De novo_ prediction and annotation

- _de novo_ prediction of LTR retrotransposons (nested, overlapping, or pure template)
- annotation of predicted LTR retrotransposons using [Dfam](http://dfam.org/) or [Repbase](http://www.girinst.org/repbase/) as reference
- solo LTR prediction
- copy number estimation of LTR elements and LTR retrotransposons
- open reading frame prediction in LTR retrotransposons
- LTR retrotransposon clustering based on DNA sequence
- age estimation of predicted LTR retrotransposons in Mya
- Methylation mark quantification in predicted LTR retrotransposons (CHH, CHG, CG, ... content)
- filtering for (potentially) active LTR retrotransposons  
- quality assesment of input genomes used to predict LTR retrotransposons

### Meta-Genomics Analyses

- run `LTRpred` on entire kingdoms of life
- perform meta genomics studies customized for LTR retrotransposons
- cluster LTR retrotransposons within and between species
- quantify the diversity space of LTR retrotransposons for entire kingdoms of life

### Visualization and Analytics Framework

- visualize the properties of predicted LTR retrotransposons
- visualize correlation between geome size and LTR retrotransposon abundance for entire kingdoms of life

## Tutorials

These tutorials introduce users to `LTRpred`:

- [Install Prerequisite Tools](https://github.com/HajkD/LTRpred/blob/master/vignettes/Installation.Rmd)
- [Introduction to LTRpred](https://github.com/HajkD/LTRpred/blob/master/vignettes/Introduction.Rmd)
- [De Novo Annotation with LTRpred](https://github.com/HajkD/LTRpred/blob/master/vignettes/Annotation.Rmd)
- [Analysis and Visualization of Predicted LTRs](https://github.com/HajkD/LTRpred/blob/master/vignettes/Analysis.Rmd)
- [Perform Meta-Genomics Studies with LTRpred]()

## NEWS

The current status of the package as well as a detailed history of the
functionality of each version of `LTRpred` can be found in the [NEWS](https://github.com/HajkD/LTRpred/blob/master/NEWS.md) section.


## Installation

Users can download `LTRpred` from [CRAN](https://cran.r-project.org/web/packages/LTRpred/index.html) :

```r
# install LTRpred 0.0.1 from CRAN
install.packages("LTRpred", dependencies = TRUE)
```

## Getting started with `LTRpred`

Users can also read the tutorials within ([RStudio](http://www.rstudio.com/)) :

```r
# source the LTRpred package
library(LTRpred)

# look for all tutorials (vignettes) available in the LTRpred package
# this will open your web browser
browseVignettes("LTRpred")

## or as single tutorials

# open tutorial: Install
 vignette("Install", package = "LTRpred")
```

In the `LTRpred` framework users can find:

#### De Novo Annotation Functions:

* `LTRpred()` : Major pipeline to predict LTR retrotransposons in a given genome
* `LTRpred.meta` : Perform Meta-Analyses with LTRpred
* `meta.summarize()` : Summarize (concatenate) all predictions of a `LTRpred.meta()` run
* `meta.apply()` : Apply functions to meta data generated by `LTRpred()`
* `LTRharvest()` : Run LTRharvest to predict putative LTR Retrotransposons
* `LTRdigest()` : Run LTRdigest to predict putative LTR Retrotransposons

#### Sequence Clustering and Similarity Computations
* `CLUSTpred()` : Cluster Sequences with VSEARCH
* `cluster.members()` : Select members of a specific cluster
* `clust2fasta()` : Export sequences of TEs belonging to the same cluster to fasta files
* `AllPairwiseAlign()` : Compute all pairwise (global) alignments with VSEARCH
* `filter.uc()` : Filter for cluster members
* `SimMatAbundance()` : Compute histogram shape similarity between species


#### LTR Copy Number Estimation

* `ltr.cn()` : Detect solo LTR copies of predicted LTR transposons
* `cn2bed()` : Write copy number estimation results to BED file format.

#### Filter Functions

* `filter.jumpers()` : Detect LTR retrotransposons that are potential jumpers
* `tidy.datasheet()` : Select most important columns of 'LTRpred' output for further analytics

#### Import the Output Files of the Prediction Tools:

* `read.prediction()` : Import the output of LTRharvest or LTRdigest
* `read.tabout()` : Import information sheet returned by LTRdigest
* `read.orfs()` : Read output of `ORFpred()`
* `read.seqs()` : Import sequences of predicted LTR transposons
* `read.ltrpred()` : Import the data sheet file generated by `LTRpred()`
* `read.uc()` : Read file in USEARCH cluster format
* `read.blast6out()` : Read file in blast6out format generated by USEARCH or VSEARCH

#### Export the Output Files of the Prediction Tools:

* `pred2bed()` : Format LTR prediction data to BED file format
* `pred2fasta()` : Save the sequence of the predicted LTR Transposons in a fasta file
* `pred2gff()` : Format LTR prediction data to GFF3 file format
* `pred2annotation()` : Match LTRharvest, LTRdigest, or LTRpred prediction with a given annotation file in GFF3 format
* `pred2csv()` : Format LTR prediction data to CSV file format

#### Visualization and Analytics Tools:

* `ORFpred()` : Open Reading Frame prediction in putative LTR transposons
* `PlotLTRAge()` : Plot the age distribution of predicted LTR transposons
* `PlotLTRWidth()` : Plot the width distribution of putative LTR transposons or LTRs
* `PlotLTRRange()` : Plot Genomic Ranges of putative LTR transposons
* `PlotSimCount()` : Plot LTR Similarity vs. predicted LTR count
* `PlotSizeCorrelation()` : Plot Genome size vs. LTR transposon count
* `PlotJumperSizeCorrelation()` : Plot Genome size vs. LTR transposon count for jumpers
* `PlotFamily()` : Visualize the Superfamily distribution of predicted LTR retrotransposons
* `PlotProteinDomain()` : Visualize the Protein Domain distribution of predicted LTR retrotransposons
* `PlotCopyNumber()` : Plot correlation between LTR copy number and methylation context  
* `PlotCluster()` : Plot correlation between Cluster Number and any other variable
* `PlotInterSpeciesCluster()` : Plot inter species similarity between TEs (for a specific cluster)
* `PlotMainInterSpeciesCluster()` : Plot inter species similarity between TEs (for the top n clusters)

#### Annotation and Validation:

* `dfam.query()` : Annotation of `de novo` predicted LTR transposons via [Dfam](http://dfam.org/help/tools) searches
* `read.dfam()` : Import Dfam Query Output
* `repbase.clean()` : Clean the initial Repbase database for BLAST
* `repbase.query()` : Query the RepBase to annotate putative LTRs
* `repbase.filter()` : Filter the Repbase query output

#### Methylation Context Estimation

* `motif.count()` : Low level function to detect motifs in strings


#### Minor helper functions

* `bcolor()` : Beautiful colors for plots
* `file.move()` : Move folders from one location to another
* `get.pred.filenames()` : Retrieve file names of files genereated by LTRpred
* `get.seqs()` : Quickly retrieve the sequences of a 'Biostrings' object
* `ws.wrap.path()` : Wrap whitespace in paths
* `rename.fasta()` : rename.fasta

## Developer Version of `LTRpred`

The developer version of `LTRpred` might include more functionality than the stable version on CRAN.
Hence users can download the current developer version of `LTRpred` by typing:

```r
# The developer version can be installed directly from github:

# install.packages("devtools")

# install developer version of LTRpred
library(devtools)
install_github("HajkD/LTRpred", build_vignettes = TRUE, dependencies = TRUE)

# On Windows, this won't work - see ?build_github_devtools
# install_github("HajkD/LTRpred", build_vignettes = TRUE, dependencies = TRUE)

# When working with Windows, first you need to install the
# R package: rtools -> http://cran.r-project.org/bin/windows/Rtools/
# or consult: http://www.rstudio.com/products/rpackages/devtools/

# Afterwards you can install devtools -> install.packages("devtools")
# and then you can run:

devtools::install_github("HajkD/LTRpred", build_vignettes = TRUE, dependencies = TRUE)

# and then call it from the library
library("LTRpred", lib.loc = "C:/Program Files/R/R-3.1.1/library")

```

## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions
provided in this package.

Furthermore, in case you find some bugs or need additional (more flexible) functionality of parts
of this package, please let me know:

https://github.com/HajkD/LTRpred/issues


## Acknowledgement

I would like to thank the [Paszkowski team](http://www.slcu.cam.ac.uk/research/paszkowski-group/group-members) for incredible support and motivating discussions that led to 
the realization of this project.



