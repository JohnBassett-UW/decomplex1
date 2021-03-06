% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dplexutilities.r
\name{barUMAP}
\alias{barUMAP}
\title{Performs dimension reduction on a formatted sample barcode table}
\usage{
barUMAP(bar.table, normalize = F, n_neighbors = 50, min_dist = 0.1)

barUMAP(bar.table, normalize = F, n_neighbors = 50, min_dist = 0.1)
}
\arguments{
\item{bar.table}{formatted sample barcode table}

\item{normalize}{logical, should the return table contain CLR normalized values. defaults to False.}

\item{n_neighbors}{n neighbors parameter for the umap function from package uwot}

\item{min_dist}{minimum distance parameter for the umap function from package uwot}
}
\value{
sample barcode table with UMAP dimensions UMAP1 & UMAP2 inserted as first two columns.

sample barcode table with UMAP dimensions UMAP1 & UMAP2 inserted as first two columns.
}
\description{
Performs dimension reduction on a formatted sample barcode table

Performs dimension reduction on a formatted sample barcode table
}
\examples{
bar.UMAP <- barUMAP(bar.table, normalize =F, n_neighbors =50, min_dist = 0.1)
bar.UMAP <- barUMAP(bar.table)
bar.UMAP <- barUMAP(bar.table, normalize =F, n_neighbors =50, min_dist = 0.1)
bar.UMAP <- barUMAP(bar.table)
}
