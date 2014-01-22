#' Alter the fishing mortality (F) values in an \code{ss3.par} file
#'
#' Takes an SS3 \code{.par} file and changes the F values for specified years.
#'
#' @author Curry James Cunningham
#'
#' @param years Vector of years for which F values are specified
#' @param years_alter Vector of years for the which F values will be altered
#' @param fvals Vector of F values to be entered into \code{ss3.par} file
#' @param file_in Input SS3 par file.
#' @param file_out Output SS3 par file.
#' @return A modified SS3 \code{.par} file.
#' @examples
#' # Create a temporary folder for the output:
#' temp_path <- file.path(tempdir(), "ss3sim-f-example")
#' dir.create(temp_path, showWarnings = FALSE)
#'
#' # Find the example .par file in the package data:
#' d <- system.file("extdata", package = "ss3sim")
#' par_file <- paste0(d, "/change_f/ss3.par")
#'
#' change_f(years = 1:49, years_alter = 2, fvals = 9999, file_in =
#' par_file, file_out = paste0(temp_path, "/test.par"))
#' @export

change_f <- function(years, years_alter, fvals, file_in="ss3.par",
  file_out="ss3.par") {

  n.years_alter <- length(years_alter)

  # Check that sufficient F values are supplied
  if(n.years_alter != length(fvals)) {
    stop(paste('#ERROR: Number of years to alter:', n.years_alter,
        'Does NOT equal length of supplied Fvalues:', length(fvals)))
  }

  # Read in ss3.par file
  ss3.par <- readLines(file_in) # Original
  ss3.par.new <- ss3.par # New file

  for(y in 1:n.years_alter) {
  	temp.year <- which(years == years_alter[y])
    temp.loc <- which(ss3.par == paste('# F_rate[', temp.year, ']:', sep=''))
    ss3.par.new[temp.loc+1] <- fvals[y]
  }

  # Write new .par file
  writeLines(ss3.par.new, con=file_out)
  close(file(file_out))
  invisible(ss3.par.new)
}
