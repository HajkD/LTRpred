# LTRpred

### De Novo LTR Prediction and Annotation with R

Transposons ... LTR transposons are a type of transposons that have the capacity ...
The `LTRpred` package aims to provide an integrated software framework to 
predict LTR transposons _de novo_ in any genomic sequence of interest.

## Tutorials

These tutorials introduce users to `LTRpred`:

- [Install Prerequisite Tools]()
- [Introduction to LTRpred]()
- [De Novo Annotation with LTRpred]()
- [Analysis and Visualization of Predicted LTRs]()


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

* `LTRharvest()` : Run LTRharvest to predict putative LTR Retrotransposons
* `LTRdigest()` : Run LTRdigest to predict putative LTR Retrotransposons


#### Import the Output Files of the Prediction Tools:

* `read.prediction()` : Import the output of LTRharvest or LTRdigest
* `read.tabout()` : Import information sheet returned by LTRdigest


#### Visualization and Analytics Tools:

* `PlotLTRAgeDistribution()` : Plot the age distribution of predicted LTR transposons
* `PlotLTRTransposonWidthDistribution()` : Plot the width distribution of putative LTR transposons
* `PlotLTRWidthDistribution()` : Plot the LTR width distribution
* `PlotRanges()` : Plot Ranges of an genomic feature

#### Annotation and Validation:
* `CleanRepBase()` : Clean the initial Repbase database for BLAST
* `QueryRepBase()` : Query the RepBase to annotate putative LTRs
* `FilterRepBaseQueryOutput()` : Filter the Repbase query output
* `MatchPredictionWithAnnotation()` : Match LTRharvest or LTRdigest prediction with a given Annotation file in GFF3 format


#### Sequence Innput/Output

* `ReadLTRharvestPredictionSeqs()` : Import sequences of LTRharvest predicted LTR transposons
* `WritePredictionToFastA()` : Save the sequence of the predicted LTR Transposons in a fasta file

#### Minor helper functions
* `bcolor()` : Beautiful colors for plots


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



