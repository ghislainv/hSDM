\name{frogs}

\alias{frogs}

\docType{data}

\title{
  Counts of the number of frogs in a water body
}

\description{
  Counts of the number of frogs in ponds of the Canton Aargau, Switzerland.
}

\format{
  A data frame with 481 observations on the following 10 variables.
  \describe{
    \item{\code{count1}}{number of counted frogs during the first visit}
    \item{\code{count2}}{number of counted frogs uring the second visit}
    \item{\code{elevation}}{elevation, meters above sea level}
    \item{\code{year}}{year}
    \item{\code{fish}}{presence of fish (1 = present, 0 = absent)}
    \item{\code{waterarea}}{area of the water body in square meters}
    \item{\code{vegetation}}{indicator of vegetation (1 = vegetation present, 0 = no vegetation present)}
    \item{\code{pondid}}{name of the pond, corresponds to observation id}
    \item{\code{x}}{x coordinate}
    \item{\code{y}}{y coordinate}
  }
}

\details{
  The amphibian monitoring program started in 1999 and is mainly aimed
  at surveying population trends of endangered amphibian species. Every
  year, about 30 water bodies in two or three randomly selected priority
  areas (out of ten priority areas of high amphibian diversity) are
  surveyed. Additionally, a random selection of water bodies that
  potentially are suitable for one of the endangered amphibian species
  but that do not belong to the priority areas were surveyed. Each water
  body is surveyed by single trained volunteer during two nocturnal
  visits per year. Volunteers recorded anurans by walking along the
  water's edge with precise rules for the duration of a survey taking
  account of the size of the surveyed water body and noting visual
  encounters and calls. As fare as possible, encountered individuals of
  the Pelophylax-complex were identified as Marsh Frog (Pelophylax
  ridibundus), Pool Frog (P. lessonaea) or hybrids (P. esculentus) based
  on morphological characteristics or based on their calls. In the given
  data set, however, these three taxa are lumped together.
}

\source{The data is provided by Isabelle Floess, Landschaft und Gewaesser, Kanton Aargau.}

\references{
Schmidt, B. R., 2005: Monitoring the distribution of pond-breeding
amphibians, when species are detected imperfectly. - Aquatic
conservation: marine and freshwater ecosystems 15: 681-692.

Tanadini, L. G.; Schmidt, B. R., 2011: Population size influences
amphibian detection probability: implications for biodiversity
monitoring programs. - Plos One 6: e28244.
}

\examples{
  data(frogs)
}

\keyword{datasets}
