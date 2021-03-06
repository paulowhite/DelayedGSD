% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FinalEstimate.R
\name{FinalEstimate}
\alias{FinalEstimate}
\title{calculate de-biased estimate}
\usage{
FinalEstimate(Info.d, Info.i, ck, lk, uk, sided = 1, kMax, estimate)
}
\arguments{
\item{Info.d}{Information at all decision analyses up to stage where study was stopped (should include information at final analysis if stopped at final analysis)}

\item{Info.i}{Information at all interim analyses up to stage where study was stopped}

\item{ck}{decision boundaries for all decision analyses up to stage where study was stopped (should include final boundary if stopped at final analysis)}

\item{lk}{lower bounds up to stage where study was stopped}

\item{uk}{upper bounds up to stage where study was stopped}

\item{sided}{one or two sided}

\item{kMax}{maximum number of analyses}

\item{estimate}{naive estimate (e.g. using  ML or REML).}
}
\description{
calculate de-biased estimate
}
