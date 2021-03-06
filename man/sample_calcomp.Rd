% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_calcomp.R
\name{sample_calcomp}
\alias{sample_calcomp}
\title{Sample conditional age-at-length (CAL) data and write to file for use by
the EM.}
\usage{
sample_calcomp(dat_list, outfile = NULL, fleets = c(1, 2), years,
  Nsamp)
}
\arguments{
\item{dat_list}{An SS data list object as read in from
\code{\link[r4ss]{SS_readdat}} in the \pkg{r4ss} package. Make sure
you select option \code{section=2}.}

\item{outfile}{A character string specifying the file name to use
when writing the information to the disk. The string must include
the proper file extension. No file is written using the default value
of \code{NULL}, which leads to increased speed because writing the 
file takes time and computing resources.}

\item{fleets}{*A vector of integers specifying which fleets to include.
The order of the fleets pertains to the input order of other arguments.
An entry of \code{fleets=NULL} leads to zero samples for any fleet.}

\item{years}{*A list the same length as \code{fleets} giving the years as
numeric vectors. If no fleet collected samples, keep the value to
\code{years=NULL}.}

\item{Nsamp}{*A numeric list of the same length as \code{fleets}.
Either single values or vectors of the same length as the number of
years can be passed through. Single values are repeated for all
years. If no fleet collected samples, keep the value to
\code{Nsamp=NULL}.}
}
\value{
A modified \code{.dat} file if \code{!is.null(outfile)}. A list object
containing the modified \code{.dat} file is returned invisibly.
}
\description{
Sample conditional age-at-length (CAL) data and write to file for use by
the EM.
}
\details{
Take a \code{data.SS_new} file containing expected values and
sample from true lengths, using length comp sample sizes, to get
realistic sample sizes for age bins given a length. Only the multinomial
distribution is currently implemented. xIf no fish are sampled then that
row is discarded. A value of NULL for fleets indicates to delete the
data so the EM If used with \code{\link{run_ss3sim}} the case file
should be named \code{calcomp}.
}
\note{
This function is only reliable when using multinomial length
compositions for the matching fleet. The real-valued length compositions
resulting from the Dirichlet distribution cause difficulties in the
sampling code. See the vignette for more.
}
\section{Which arguments to specifiy in case files}{

All function argument descriptions that start with an asterisk (*) will be passed
through the case files to \code{\link{run_ss3sim}}. If one of these arguments
is not specified in a case file, then a value of \code{NULL} will be passed,
which may or may not be an appropriate value. Other arguments will be ignored
if specified.
}

\seealso{
Other sampling functions: \code{\link{clean_data}},
  \code{\link{sample_agecomp}}, \code{\link{sample_index}},
  \code{\link{sample_lcomp}}, \code{\link{sample_mlacomp}},
  \code{\link{sample_wtatage}}
}
\author{
Cole Monnahan, Kotaro Ono
}
\concept{sampling functions}
