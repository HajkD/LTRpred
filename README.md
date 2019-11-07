## LTRpred(ict): a pipeline for automated functional annotation of LTR retrotransposons for comparative genomics studies  <img src="inst/LTRpred_logo.png" align="right" height="174" width="150" />


An easy way to perform __de novo__ functional annotation of `LTR retrotransposons` from any genome assembly in `fasta` format.

![](vignettes/LTRfeatures.png)

## Install

```r
# install the current version of LTRpred on your system
source("http://bioconductor.org/biocLite.R")
biocLite("devtools")
biocLite("HajkD/LTRpred")
```
## Tutorials

### Quick Start

The fastest way to generate a LTR retrotransposon prediction for a genome of interest (after [installing](https://hajkd.github.io/LTRpred/articles/Introduction.html) all prerequisite command line tools) is to use the
`LTRpred()` function and relying on the default parameters. In the following example,
a LTR transposon prediction is performed for parts of the Human Y chromosome.

```r
# load LTRpred package
library(LTRpred)
# de novo LTR transposon prediction for the Human Y chromosome
LTRpred(genome.file = system.file("Hsapiens_ChrY.fa", package = "LTRpred"))
```

When running your own genome, please specify `genome.file = "path/to/your/genome.fasta` instead of `system.file(..., package = "LTRpred")`. The command `system.file(..., package = "LTRpred")` merely references the path to the example file stored in the LTRpred package itself.


## Citation
The `LTRpred` package is not formally published yet, but a manuscript is in preparation. For now, please cite one of the the following paper when using `LTRpred` for your own research. `LTRpred` is part of these studies and helped to predict potentially active retrotransposons that were later confirmed experimentally.

> M Benoit, __HG Drost__, M Catoni, Q Gouil, S Lopez-Gomollon, DC Baulcombe, J Paszkowski. [__Environmental and epigenetic regulation of Rider retrotransposons in tomato__](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008370). _PloS Genetics_ (2019) (__in press__). 

or

> J Cho, M Benoit, M Catoni, __HG Drost__, A Brestovitsky, M Oosterbeek and J Paszkowski.  [__Sensitive detection of pre-integration intermediates of LTR retrotransposons in crop plants__](https://www.nature.com/articles/s41477-018-0320-9). _Nature Plants_, 5,  26-33 (2019).


This tutorial introduces users to `LTRpred`:

- [Introduction to LTRpred](https://hajkd.github.io/LTRpred/articles/Introduction.html)

Users can also read the tutorials within ([RStudio](http://www.rstudio.com/)) :

```r
library(LTRpred)
browseVignettes("LTRpred")
```

### Studies that successfully used `LTRpred` to annotate functional retrotransposons

> - J Cho, M Benoit, M Catoni, __HG Drost__, A Brestovitsky, M Oosterbeek and J Paszkowski.  [__Sensitive detection of pre-integration intermediates of LTR retrotransposons in crop plants__](https://www.nature.com/articles/s41477-018-0320-9). __Nature Plants__, 5,  26-33 (2019).
>
> - M Benoit, __HG Drost__, M Catoni, Q Gouil, S Lopez-Gomollon, DC Baulcombe, J Paszkowski. [__Environmental and epigenetic regulation of Rider retrotransposons in tomato__](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008370). __PloS Genetics__ (2019) (__in press__). 
>
> - E Cerruti, C Gisbert, __HG Drost__, D Valentino, E Portis, L Barchi, J Prohens, S Lanteri, C Comino,  M Catoni. [__Epigenetic bases of grafting-induced vigour in eggplant__](https://www.biorxiv.org/content/10.1101/831719v1). __bioaRxiv__ (2019).
>
> - Nguinkal _et al._ [__The First Highly Contiguous Genome Assembly of Pikeperch (Sander lucioperca), an Emerging Aquaculture Species in Europe__](https://www.mdpi.com/2073-4425/10/9/708/htm) __Genes__, 0(9), 708 (2019).

## Discussions and Bug Reports

I would be very happy to learn more about potential improvements of the concepts and functions
provided in this package.

Furthermore, in case you find some bugs or need additional (more flexible) functionality of parts
of this package, please let me know:

https://github.com/HajkD/LTRpred/issues

In the `LTRpred` framework users can find:

### _De novo_ prediction and annotation

- _de novo_ prediction of LTR retrotransposons (nested, overlapping, or pure template) using [LTRharvest](http://www.zbh.uni-hamburg.de/forschung/arbeitsgruppe-genominformatik/software/ltrharvest.html) and [LTRdigest](http://www.zbh.uni-hamburg.de/forschung/gi/software/ltrdigest.html)
- annotation of predicted LTR retrotransposons using [Dfam](http://dfam.org/) or [Repbase](http://www.girinst.org/repbase/) as reference
- solo LTR prediction based on specialized [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) searches
- LTR retrotransposons family clustering using [vsearch](https://github.com/torognes/vsearch)
- open reading frame prediction in LTR retrotransposons using [usearch](https://www.drive5.com/usearch/)
- age estimation of predicted LTR retrotransposons in Mya (not implemented yet, but soon to come..)
- CHH, CHG, CG, ... content quantification in predicted LTR retrotransposons
- filtering for (potentially) functional LTR retrotransposons  
- quality assesment of input genomes used to predict LTR retrotransposons

### Meta-Genomics Analyses

- run `LTRpred` on entire kingdoms of life using only one command (see `?LTRpred.meta`)
- perform meta genomics studies customized for LTR retrotransposons
- cluster LTR retrotransposons within and between species
- quantify the diversity space of LTR retrotransposons for entire kingdoms of life

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

#### Analytics Tools:

* `ORFpred()` : Open Reading Frame prediction in putative LTR transposons

#### Annotation and Validation:

* `dfam.query()` : Annotation of `de novo` predicted LTR transposons via [Dfam](http://dfam.org/help/tools) searches
* `read.dfam()` : Import Dfam Query Output
* `repbase.clean()` : Clean the initial Repbase database for BLAST
* `repbase.query()` : Query the RepBase to annotate putative LTRs
* `repbase.filter()` : Filter the Repbase query output

#### Methylation Context Estimation

* `motif.count()` : Low level function to detect motifs in strings

#### Visualization Framework

* `plot_ltrsim_individual()` : Plot the age distribution of predicted LTR transposons
* `plot_ltrwidth_individual()` : Plot the width distribution of putative LTR transposons or LTRs for individual species
* `plot_ltrwidth_species()` : Plot the width distribution of putative LTR transposons or LTRs for all species
* `plot_ltrwidth_kingdom()` : Plot the width distribution of putative LTR transposons or LTRs for all kingdoms
* `plot_copynumber_individual()` : Plot the copy number distribution of putative LTR transposons or LTRs for individual species
* `plot_copynumber_species()` : Plot the copy number distribution of putative LTR transposons or LTRs for all species
* `plot_copynumber_kingdom()` : Plot the copy number distribution of putative LTR transposons or LTRs for all kingdoms
* `plotLTRRange()` : Plot Genomic Ranges of putative LTR transposons
* `PlotSimCount()` : Plot LTR Similarity vs. predicted LTR count
* `plotSize()` : Plot Genome size vs. LTR transposon count
* `plotSizeJumpers()` : Plot Genome size vs. LTR transposon count for jumpers
* `plotFamily()` : Visualize the Superfamily distribution of predicted LTR retrotransposons
* `plotDomain()` : Visualize the Protein Domain distribution of predicted LTR retrotransposons
* `plotCN()` : Plot correlation between LTR copy number and methylation context  
* `plotCluster()` : Plot correlation between Cluster Number and any other variable
* `PlotInterSpeciesCluster()` : Plot inter species similarity between TEs (for a specific cluster)
* `PlotMainInterSpeciesCluster()` : Plot inter species similarity between TEs (for the top n clusters)


#### Minor helper functions

* `bcolor()` : Beautiful colors for plots
* `file.move()` : Move folders from one location to another
* `get.pred.filenames()` : Retrieve file names of files genereated by LTRpred
* `get.seqs()` : Quickly retrieve the sequences of a 'Biostrings' object
* `ws.wrap.path()` : Wrap whitespace in paths
* `rename.fasta()` : rename.fasta


## Acknowledgement

I would like to thank the [Paszkowski team](http://www.slcu.cam.ac.uk/research/paszkowski-group/group-members) for incredible support and motivating discussions that led to 
the realization of this project.



