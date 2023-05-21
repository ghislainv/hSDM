# hSDM R Package <img src="man/figures/logo.png" align="right" alt="" width="120" />

[![R-CMD-check](https://github.com/ghislainv/hSDM/workflows/R-CMD-check/badge.svg)](https://github.com/ghislainv/hSDM/actions/workflows/check-standard.yaml)
[![CRAN Status](https://www.r-pkg.org/badges/version/hSDM)](https://cran.r-project.org/package=hSDM)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.594920.svg)](https://doi.org/10.5281/zenodo.594920)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/hSDM)](https://cran.r-project.org/package=hSDM)

`hSDM` is an R package for estimating parameters of hierarchical Bayesian species distribution models. Such models allow interpreting the observations (occurrence and abundance of a species) as a result of several hierarchical processes including ecological processes (habitat suitability, spatial dependence and anthropogenic disturbance) and observation processes (species detectability). Hierarchical species distribution models are essential for accurately characterizing the environmental response of species, predicting their probability of occurrence, and assessing uncertainty in the model results.

## Installation

Install the latest stable version of `hSDM` from [CRAN](https://cran.r-project.org/) with:


```r
install.packages("hSDM")
```

Or install the development version of `hSDM` from [GitHub](https://github.com/ghislainv/hSDM) with:


```r
devtools::install_github("ghislainv/hSDM")
```

## Vignettes and manual

- Presentation at ISEC 2014: [hSDM-ISEC2014.pdf](https://ecology.ghislainv.fr/hSDM/pdfs/hSDM-ISEC2014.pdf)
- Long vignette with several examples: [hSDM-vignette.pdf](https://ecology.ghislainv.fr/hSDM/pdfs/hSDM-vignette.pdf)
- Manual: [hSDM-manual.pdf](https://ecology.ghislainv.fr/hSDM/pdfs/hSDM-manual.pdf)

## In the wild

- Tutorial on using opportunistic species occurrence data for occupancy modelling by [Adam M. Wilson](https://github.com/adammwilson) on [GitHub](https://github.com/adammwilson/hSDM_Tutorial/blob/master/hSDM_Tutorial.md).
- Tutorial on modelling spatial autocorrelation by [Jérôme Guélat](https://www.random-nature.net) on [Amazonaws](https://rstudio-pubs-static.s3.amazonaws.com/9687_cc323b60e5d542449563ff1142163f05.html).

## Contributing

The `hSDM` R package is Open Source and released under the [GNU GPL version 3](https://www.gnu.org/licenses/gpl-3.0.en.html) license. Anybody who is interested can contribute to the package development following our [Contributing guide](https://ecology.ghislainv.fr/hSDM/CONTRIBUTING.html). Every contributor must agree to follow the project's [Code of conduct](https://ecology.ghislainv.fr/hSDM/CODE_OF_CONDUCT.html).

## References

**Diez J. M. and Pulliam H. R.** 2007. Hierarchical analysis of species distributions and abundance across environmental gradients. _Ecology_. **88**(12): 3144-3152.

**Gelfand A. E., Silander J. A., Wu S. S., Latimer A., Lewis P. O., Rebelo A. G. and Holder M.** 2006. Explaining species distribution patterns through hierarchical modeling. _Bayesian Analysis_. **1**(1): 41-92.

**Latimer, A. M.; Wu, S. S.; Gelfand, A. E. & Silander, J. A.** 2006. Building statistical models to analyze species distributions. _Ecological Applications_. **16**(1): 33-50.

**MacKenzie, D. I.; Nichols, J. D.; Lachman, G. B.; Droege, S.; Andrew Royle, J. and Langtimm, C. A.** 2002. Estimating site occupancy rates when detection probabilities are less than one. _Ecology_. **83**: 2248-2255.

**Royle, J. A.** 2004. N-mixture models for estimating population size from spatially replicated counts. _Biometrics_. **60**: 108-115.
