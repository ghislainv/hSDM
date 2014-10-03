\name{hSDM.ncWriteOutput}
\alias{hSDM.ncWriteOutput}

\title{Export hSDM output to a netCDF file for archiving and model comparison}

\description{The \code{hSDM.nc} function performs convergence and model performance diagnostics on hSDM output (currently only hSDM.ZIB).}

\usage{hSDM.nc(results,species,modelname,model,data,fdata,outputdir,autocor=T,verbose=T)}

\arguments{
  
  \item{results}{A list containing multiple chains of hSDM output}

  \item{species}{A character vector indicating the species name (used as metadata)}
  
  \item{modelname}{A character vector indicating the model formula (used as metadata)}

  \item{model}{A character vector indicating a model name (e.g. "full") to be used as metadata}

  \item{data}{The data frame used for model prediction}
    
  \item{fdata}{The data frame used for model fitting}
    
  \item{outputdir}{The target directory of the netCDF file} 

  \item{autocor}{Logical indicating whether to calculate the spatial autocorrelation of the predicted p(occurrence) surface}

  \item{verbose}{logical indicating whether to print status updates}
}

\value{
  
  \item{modeloutput.nc}{The function writes a netCDF file to \code{outputdir}.}
    
}

\details{

  
}

\references{

}

\author{
Adam M. Wilson (adam.wilson@yale.edu)
}

\seealso{
  \code{\link[hSDM]{hSDM.ZIB}}, \code{\link[coda]{summary.mcmc}}
}



\keyword{Poisson regression}
\keyword{biodiversity}
\keyword{species distribution models}
\keyword{hierarchical Bayesian models}
\keyword{MCMC}
\keyword{Markov Chains Monte Carlo}
\keyword{Metropolis algorithm}
