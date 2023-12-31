% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IC3.R
\name{IC3}
\alias{IC3}
\title{Title}
\usage{
IC3(A, cellinfo, lrinfo, alpha = 0.02)
}
\arguments{
\item{A}{Gene expression count matrix, each row represents a cell and each colume represents a gene.}

\item{cellinfo}{Cell information matrix, each row represents a cell.The first column represents the x-axis coordinates of the cell, and the second column represents the y-axis coordinates of the cell. The third column represents the cell type of the cell. The row name is the name of the cell, and its order should be consistent with the order of gene expression matrix A.}

\item{lrinfo}{lrinfo: Ligand-Receptor Database, each row represents a ligand-receptor pair. Ligand is in the first column.}

\item{alpha}{Distance proportional threshold, cells with a distance exceeding the threshold are considered to have no interaction with each other. Default value is 0.02.}
}
\value{
The first term: the communication probability with cell type level; The second term: the communication probability with single cell level. The 3-5th term: parameter estimation of lambda;beta;r=(r0,r1,r2).
}
\description{
Title
}
\examples{
IC3(A,cellinfo,lrinfo,alpha=0.02)
}
