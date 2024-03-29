% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/ancom_functions.r
\name{ANCOM}
\alias{ANCOM}
\title{Run the ANCOM method}
\usage{
ANCOM(OTUdat, sig = 0.05, multcorr = 3, tau = 0.02, theta = 0.1,
  repeated = FALSE)
}
\arguments{
\item{OTUdat}{the OTU dataset. See Details for formatting instructions.}

\item{sig}{the significance level (or FDR) at which the tests are to be performed.}

\item{multcorr}{type of correction for multiple comparisons, see Details.}

\item{tau}{a tuning parameter in the Stepwise testing method. See Details.}

\item{theta}{a tuning parameter in the Stepwise testing method. See Details.}

\item{repeated}{logical determining whether the data have repeated measures (e.g., longitudinal design).}
}
\value{
The function produces a list with the following elements:
\itemize{
\item{ \code{W}: }{ values of the test statistics.}
\item{ \code{detected}: }{ names of OTUs detected.}
\item{ \code{dframe}: }{ the input dataframe.}
}
}
\description{
Runs ANCOM to test for differences in OTU abundances between treatment groups.
}
\details{
The ANCOM method was developed and tested with default values of the two tuning parameters
(\code{tau=0.02} and \code{theta=0.1}). For consistency, users are recommended to leave
these tuning parameters at their default values, unless they wish to explore the performance
 of ANCOM for different values of the tuning parameters.

Data should be formatted as follows: each row is a subject, and each column is an OTU.
The final column should contain the grouping variable.

To adjust for multiple testing, \code{multcorr} may take on the following values:
\itemize{
\item{ \code{1}: }{ A stringent correction}
\item{ \code{2}: }{ A less stringent correction}
\item{ \code{3}: }{ No correction (default)}
}
The more stringent correction is not available in the shiny application.
}
\note{
The function \code{\link{plot_ancom}} will produce plots for objects produced by \code{ANCOM}.
}
\examples{
\dontrun{
## Create and run a small example

nn <- 10
pp <- 20
sim_otu <- matrix( 0, nrow=nn, ncol=pp+1 )
sim_otu <- data.frame(sim_otu)
colnames(sim_otu) <- c( paste0("OTU_", letters[1:pp] ), "Group" )
sim_otu[,pp+1]    <- c( rep("Control",nn/2), rep("Treatment",nn/2)  )
idx_trt <- sim_otu$Group=="Treatment"

for( ii in 1:pp ){
  sim_otu[,ii] <- rpois( nn, 1 )
}

# Create some significance
sim_otu[idx_trt,3] <- rpois( nn/2, 8)
sim_otu[idx_trt,7] <- rpois( nn/2, 8)
sim_otu[idx_trt,9] <- rpois( nn/2, 8)

ancom.out <- ANCOM( OTUdat = sim_otu, sig = 0.20, multcorr = 2 )
ancom.out$W
ancom.out$detected
}
}

