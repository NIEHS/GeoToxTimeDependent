#' GeoTox Data
#'
#' Sample data for use in vignettes and function examples. See the Package
#' Data vignette, \code{vignette("package_data", package = "GeoTox")}, for
#' details on how this data was gathered.
#'
#' @format A list with items:
#' \describe{
#'   \item{exposure}{2019 AirToxScreen exposure concentrations for a subset of
#'   chemicals in North Carolina counties.}
#'   \item{dose_response}{Subset of chemicals curated by ICE cHTS as active
#'   within a set of assays.}
#'   \item{age}{County population estimates for 7/1/2019 in North Carolina.}
#'   \item{obesity}{CDC PLACES obesity data for North Carolina counties in
#'   2020.}
#'   \item{simulated_css}{Simulated steady-state plasma concentrations for
#'   various age groups and obesity status combinations.}
#'   \item{boundaries}{County and state boundaries for North Carolina in 2019.}
#' }
"geo_tox_data"


#' Chemical API Server url
#'
#' A section of url used in Chemical API Endpoints
'chemical_api_server'

#' Hazard API Server url
#'
#' A section of url used in Hazard API Endpoints
'hazard_api_server'

#' Bioactivity API Server url
#'
#' A section of url used in Bioactivity API Endpoints
'bioactivity_api_server'

#' Exposure API Server url
#'
#' A section of url used in Exposure API Endpoints
'exposure_api_server'

utils::globalVariables(c('chemical_api_server', 'hazard_api_server',
                         'bioactivity_api_server', 'exposure_api_server'))
