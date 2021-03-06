% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/change_data.r
\name{change_data}
\alias{change_data}
\title{Change the data that is available as output from an SS operating model.}
\usage{
change_data(dat_list, outfile = NULL, fleets, years, types,
  age_bins = NULL, len_bins = NULL, pop_binwidth = NULL,
  pop_minimum_size = NULL, pop_maximum_size = NULL,
  lcomp_constant = NULL, tail_compression = NULL, nsex = 1)
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

\item{fleets}{A numeric vector of fleets}

\item{years}{A numeric vector of years}

\item{types}{A vector that can take combinations of the following entries:
\code{"index"}, \code{"len"}, \code{"age"}, \code{"cal"}, \code{"mla"}.
\code{types} controls what data structures the function acts on, with
\code{"index"} changing indices/CPUE, \code{"len"} augmenting the length
composition data, \code{"age"} augmenting the age composition, \code{"cal"}
augmenting the conditional age at length, and \code{"mla"} augmenting the
mean length at age data.}

\item{age_bins}{*A numeric vector of age bins to use. If left as \code{NULL}
then the age bin structure will be taken from the OM.}

\item{len_bins}{*A numeric vector of length bins to use. If left as
\code{NULL} then the length bin structure will be taken from the OM.
For conditional age-at-length (CAAL) data, the last value provided to
\code{len_bins} will be used for Lbin_lo and -1 will be used for Lbin_hi
for the largest length bin category, i.e., row of CAAL data.}

\item{pop_binwidth}{*Population length bin width. Note that this value must
be smaller than the bin width specified in length composition data
\code{len_bins} or SS will fail (see notes in the SS manual).}

\item{pop_minimum_size}{*Population minimum length bin value.}

\item{pop_maximum_size}{*Population maximum length bin value.}

\item{lcomp_constant}{*A new robustification constant for length composition
data to be used. Must be a numeric value, as a proportion. For example 0.1
means 10 percent. See the SS manual for further information. A \code{NULL}
value indicates no action resulting in using the current value, and a value
of 0 will throw an error since that leads to an error when zeroes exist in
the data. Instead use a very small value like \code{1e-07}.}

\item{tail_compression}{*A new tail compression value to be used in SS. Must
be a numeric value, as a proportion. For example 0.1 means 10 percent. See
the SS manual for further information. A \code{NULL} value indicates no
action, a negative value indicates to SS to ignore it (not use that
feature).}

\item{nsex}{An integer value of 1 or 2 specifying the number of sexes in the
model. If 1, then females are the only included sex. This information can be
found in the data file for a given model and dictates how the composition
data are structured.}
}
\value{
An invisible data list, and a file is written to the disk if an
entry other than \code{NULL} is provided for \code{outfile}.
}
\description{
\code{change_data} alters the data structure for a data list as read in by
\code{\link[r4ss]{SS_readdat}}, for use in preparing the data file for an SS
operating model. Original data is removed and dummy data is added, as
specified, to the SS \code{.dat} file. This causes SS to produce expected
values (OM "truth") when the operating model is run, from which data can be
sampled.  For each data type altered, \code{change_data} will add data for
the fleets and years given; potentially adding many rows of redundant data.
Currently, \code{.dat} files with multiple sexes cannot be manipulated with
\code{change_data}. \code{\link{calculate_data_units}} is used internally in
\code{\link{ss3sim_base}} to create a superset of fleets and years from
sample arguments, and \code{\link{clean_data}} to strip out unused data after
\code{change_data} is called (see examples below). \code{change_data} is
called internally automatically, but can also be used by an \pkg{ss3sim} user
to manipulate data as a case, or to prepare a new OM for use in a simulation.
See the vignette for more details.
}
\details{
The robustification constant is added to both the observed and
  expected proportions of length composition data, before being normalized
  internally. It is designed to help stabilize the model, but is unclear how
  and when to use it for optimal effect. The same value is used for all
  length data.
}
\section{Which arguments to specifiy in case files}{

All function argument descriptions that start with an asterisk (*) will be passed
through the case files to \code{\link{run_ss3sim}}. If one of these arguments
is not specified in a case file, then a value of \code{NULL} will be passed,
which may or may not be an appropriate value. Other arguments will be ignored
if specified.
}

\examples{
d <- system.file("extdata", package = "ss3sim")
fleets <- 1:2
years <- c(5, 10, 15)
types <- c("len", "age")
file_in <- r4ss::SS_readdat(file.path(d, "models", "cod-om", "codOM.dat"),
  version = NULL, verbose = FALSE)

# Basic test with just length data, default bins:
out <- change_data(file_in, outfile = NULL, types = "len",
  years = years, fleets = fleets)
print(out$lbin_vector)
print(out$lencomp)

# Change the length bins:
out <- change_data(file_in, outfile = NULL, types = "len",
  years = years, fleets = fleets, len_bins = 3:6)
out$lbin_vector
out$lencomp

# Change the population length bins:
out <- change_data(file_in, outfile = NULL, types = "len",
  years = years, fleets = fleets, pop_binwidth = 1, pop_minimum_size = 5,
  pop_maximum_size = 210)
out$binwidth
out$maximum_size
out$minimum_size
}
\seealso{
\code{\link{sample_lcomp}}, \code{\link{sample_agecomp}}

Other change functions: \code{\link{change_em_binning}},
  \code{\link{change_e}}, \code{\link{change_f_par}},
  \code{\link{change_f}}, \code{\link{change_o}},
  \code{\link{change_retro}}, \code{\link{change_tv}}
}
\author{
Cole Monnahan, Ian Taylor, Sean Anderson, Kelli Johnson
}
\concept{change functions}
