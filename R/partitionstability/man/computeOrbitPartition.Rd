% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orbits.R
\name{computeOrbitPartition}
\alias{computeOrbitPartition}
\title{computeOrbitPartition}
\usage{
computeOrbitPartition(S, n)
}
\arguments{
\item{S}{Set of generators}

\item{n}{Length of the set the group \eqn{\langle S\rangle} acts on}
}
\description{
Implementation of Algorithms 2 and 3 (combined as one function), which computes the
orbit partition of a permutation group from a given set of generators.
}
\examples{
library("sets")
p1 <- c(2L, 1L, 3L, 4L, 6L, 5L, 7L, 8L, 9L, 10L)
p2 <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 10L, 9L)
p3 <- c(1L, 9L, 3L, 4L, 5L, 6L, 7L, 8L, 2L, 10L)
S <- list(p1, p2, p3)
computeOrbitPartition(S, 10)
# [1] 1 1 3 4 5 5 7 8 1 1

}
