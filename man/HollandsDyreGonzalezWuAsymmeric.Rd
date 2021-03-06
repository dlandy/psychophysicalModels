% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classicPsychophysics.R
\name{HollandsDyreGonzalezWuAsymmeric}
\alias{HollandsDyreGonzalezWuAsymmeric}
\title{Create Predicted data patterns from a 'middle out' asymmetric gonzalez & Wu model}
\usage{
HollandsDyreGonzalezWuAsymmeric(stimuli, beta = 1, delta = 1, scaling = 1,
  leftScaling = scaling, rightScaling = scaling, center = 0.5)
}
\arguments{
\item{stimuli}{a Vector of stimuli}

\item{beta}{The curvature of the space, corresponding to a stevens beta. <1 yields overestimation of small values.}

\item{delta}{The location of the prior / the elevation}

\item{scaling}{The overall distance to the priors, used when the scaling should be symmetric. The value should go from 0 to 1. 0 means that the reference point is at infinity.}

\item{leftScaling}{The scaling of the lefthand side, or lower values.}

\item{rightScaling}{The scaling of the righthand side, or higher values}
}
\description{
This model calculates, given a set of parameters, the predicted results of an 'asymmetric' HD GW model. This 
model uses a central (but maybe off-center) reference point, and applies the G&W model from that central point ot each end
Each end may be separately scaled.
}
\examples{
HollandsDyreGonzalezWuAsymmeric(c(1,2))
}
\keyword{proportionJudgments}
