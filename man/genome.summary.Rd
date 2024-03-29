% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genome.summary.R
\name{genome.summary}
\alias{genome.summary}
\title{Generating genome summary files for \code{LTRpred.meta} results}
\usage{
genome.summary(
  genome.folder,
  ltrpred.meta.folder,
  file.name,
  sim = 70,
  cut.range = 2,
  quality.filter = TRUE,
  n.orfs = 0,
  strategy = "default"
)
}
\arguments{
\item{genome.folder}{a file path to a folder storing the genome assembly files in fasta format that
were used to generate \code{\link{LTRpred}} annotations of diverse species from the same taxonomic kingdom.}

\item{ltrpred.meta.folder}{a file path to a folder storing \code{\link{LTRpred}} annotations of diverse species from the same taxonomic kingdom.}

\item{file.name}{name of the output file.}

\item{sim}{LTR similarity threshold. Only putative LTR transposons that fulfill this 
LTR similarity threshold will be retained.}

\item{cut.range}{similarity interval size.}

\item{quality.filter}{shall a quality filter to remove possible false positive predictions be applied?}

\item{n.orfs}{minimum number of open reading frames a predicted retroelement shall possess.}

\item{strategy}{quality filter strategy. Options are
\itemize{
\item \code{strategy = "default"} : see section \code{Quality Control} 
\item \code{strategy = "stringent"} : in addition to filter criteria specified in section \code{Quality Control},
the filter criteria \code{!is.na(protein_domain)) | (dfam_target_name != "unknown")} is applied
}}
}
\description{
Generating genome summary files for \code{LTRpred.meta} results
}
\details{
Generating genome summary files for \code{LTRpred.meta} results.
}
\author{
Hajk-Georg Drost
}
