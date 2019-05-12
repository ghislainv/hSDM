
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hSDM R Package <img src="man/figures/logo.png" align="right" alt="" width="120" />

[![CRAN
Status](https://www.r-pkg.org/badges/version/hSDM)](https://cran.r-project.org/package=hSDM)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.594920.svg)](https://doi.org/10.5281/zenodo.594920)
[![Downloads](https://cranlogs.r-pkg.org/badges/hSDM)](https://cran.r-project.org/package=hSDM)

`hSDM` is an R package for estimating parameters of hierarchical
Bayesian species distribution models. Such models allow interpreting the
observations (occurrence and abundance of a species) as a result of
several hierarchical processes including ecological processes (habitat
suitability, spatial dependence and anthropogenic disturbance) and
observation processes (species detectability). Hierarchical species
distribution models are essential for accurately characterizing the
environmental response of species, predicting their probability of
occurrence, and assessing uncertainty in the model results.

## Installation

``` r
# Install release version from CRAN
install.packages("hSDM")

# Install development version from GitHub
devtools::install_github("ghislainv/hSDM")
```

## Vignettes and manual

  - Presentation at ISEC 2014:
    [hSDM-ISEC2014.pdf](https://bioscenemada.cirad.fr/FileTransfer/hSDM-ISEC2014.pdf)
  - Long vignette with several examples:
    [hSDM-vignette.pdf](https://bioscenemada.cirad.fr/FileTransfer/hSDM-vignette.pdf)
  - Manual:
    [hSDM-manual.pdf](https://CRAN.R-project.org/package=hSDM/hSDM.pdf)

## In the wild

  - Tutorial on using opportunistic species occurrence data for
    occupancy modelling by [Adam M.
    Wilson](https://github.com/adammwilson) on
    [GitHub](https://github.com/adammwilson/hSDM_Tutorial/blob/master/hSDM_Tutorial.md).
  - Tutorial by Adam M. Wilson adapted by [Marta A.
    Jarzyna](https://www.majarzyna.com/) on
    [spatial-ecology.net](http://spatial-ecology.net/dokuwiki/doku.php?id=wiki:spdistr2)
  - Tutorial on modelling spatial autocorrelation by [Jérôme
    Guélat](https://www.random-nature.net) on
    [Amazonaws](https://rstudio-pubs-static.s3.amazonaws.com/9687_cc323b60e5d542449563ff1142163f05.html).

## Related publications

**Diez J. M. and Pulliam H. R.** 2007. Hierarchical analysis of species
distributions and abundance across environmental gradients. *Ecology*.
**88**(12): 3144-3152.

**Gelfand A. E., Silander J. A., Wu S. S., Latimer A., Lewis P. O.,
Rebelo A. G. and Holder M.** 2006. Explaining species distribution
patterns through hierarchical modeling. *Bayesian Analysis*. **1**(1):
41-92.

**Latimer, A. M.; Wu, S. S.; Gelfand, A. E. & Silander, J. A.** 2006.
Building statistical models to analyze species distributions.
*Ecological Applications*. **16**(1): 33-50.

**MacKenzie, D. I.; Nichols, J. D.; Lachman, G. B.; Droege, S.; Andrew
Royle, J. and Langtimm, C. A.** 2002. Estimating site occupancy rates
when detection probabilities are less than one. *Ecology*. **83**:
2248-2255.

**Royle, J. A.** 2004. N-mixture models for estimating population size
from spatially replicated counts. *Biometrics*. **60**: 108-115.
