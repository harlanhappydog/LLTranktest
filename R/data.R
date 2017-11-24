#' toenail
#'
#'

#' The data are from a randomized trial comparing two oral treatments for toenail dermatophyte onychomycosis infection, see De Backer et al. (1996). Participants were evaluated at baseline (month 0) and at 1, 2, 3, 6, 9 and 12 months.  The response variable is the unaffected nail length expressed in mm.
#'
#' Most recently, Mian and Hasan (2012) and Mahabadi (2014) have analyzed the toenail data.
#'
#'
#'
#' Additional information is available at
#'
#' http://www.aliquote.org/articles/tech/MDLD/README.html
#'
#' De Backer, M., De Keyser, P., De Vroey, C. and Lesaffre, E. (1996). A 12-week treatment for dermatophyte toe onychomycosis: terbinafine 250mg/day vs. itraconazole 200mg/day--a double-blind comparative trial. British Journal of Dermatology, 134, 16-17.
#' http://www.jaad.org/article/S0190-9622(98)70486-4/abstract
#'
#' Mian, R. I. and Hasan, M. T. (2012). Two-part pattern-mixture model for longitudinal incom- plete semi-continuous toenail data. International Journal of Statistics in Medical Research, 1(2):120–127.
#'
#' Mahabadi, S. E. (2014). A bayesian shared parameter model for incomplete semicontinuous longitudinal data: An application to toenail dermatophyte onychomycosis study. Journal of Statistical Theory and Applications, 13(4):317–332.
#'

#'
#'
#' @format A data frame with 1854 rows and 4 variables:
#' \describe{
#'   \item{ID}{Unique subject identifier}
#'   \item{treat_group}{Treatment group identifier}
#'   \item{month}{month}
#'   \item{UNL_mm}{Unaffected Nail Length (UNL) in millimeters, measured from the nail bed to the infected part of the nail.}
#' }
"toenail"
