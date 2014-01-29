\name{datacells.Latimer2006}

\alias{datacells.Latimer2006}

\docType{data}

\title{Data of presence-absence (from Latimer et al. 2006)}

\description{Data come from a small region including 476 one minute by
one minute grid cells. This region is a small corner of South
Africa's Cape Floristic Region, and includes very high plant species
diversity and a World Biosphere Reserve. The data frame can be used as an example
for several functions in the \code{hSDM} package.}

\usage{datacells.Latimer2006}

\format{
  
  \code{datacells.Latimer2006} is a data frame with 476 observations (cells) on the following 9 variables.
  \describe{
    \item{\code{y}}{the number of times the species was observed to be present in each cell}
    \item{\code{n}}{the number of visits or sample locations in each cell (which can be zero)}
    \item{\code{rough}}{elevational range or "roughness"}
    \item{\code{julmint}}{July minimum temperature}
    \item{\code{pptcv}}{interannual variation in precipitation}
    \item{\code{smdsum}}{summer soil moisture days}
    \item{\code{evi}}{enhanced vegetation or "greenness" index}
    \item{\code{ph1}}{percent acidic soil}
    \item{\code{num}}{number of neighbors of each cell, this is a sparse representation of the adjacency
      matrix for the subregion.}
  }
  
}

\source{Latimer et al. (2006) \emph{Ecological Applications}, Appendix B}

\references{

  Latimer, A. M.; Wu, S. S.; Gelfand, A. E. and Silander, J. A. (2006)
  Building statistical models to analyze species
  distributions. \emph{Ecological Applications}, 16, 33-50.

}

\keyword{datasets}
