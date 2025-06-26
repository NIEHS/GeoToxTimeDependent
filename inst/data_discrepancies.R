# ################################################################################
#
#
#
# # 1. First load libraries
# # 2. Run simulation of css values following the vignette
# # 3. Determine differences between stored data and locally simulated data
#
#
# # 1. Load libraries
# library(dplyr)
# library(sf)
# library(tidyr)
# library(readr)
# library(stringr)
# library(purrr)
# library(readxl)
# library(httk)
# library(httr2)
#
#
# # 2. Run simulations
# set.seed(2345)
# n_samples <- 500
#
# # Get CASN for which httk simulation is possible. Try using load_dawson2021,
# # load_sipes2017, or load_pradeep2020 to increase availability.
# load_sipes2017()
# casn <- intersect(unique(geo_tox_data$dose_response$casn),
#                   get_cheminfo(suppress.messages = TRUE))
#
# # Define population demographics for httk simulation
# pop_demo <- cross_join(
#   tibble(age_group = list(c(0, 2), c(3, 5), c(6, 10), c(11, 15),
#                           c(16, 20), c(21, 30), c(31, 40), c(41, 50),
#                           c(51, 60), c(61, 70), c(71, 100))),
#   tibble(weight = c("Normal", "Obese"))) |>
#   # Create column of lower age_group values
#   rowwise() |>
#   mutate(age_min = age_group[1]) |>
#   ungroup()
#
# # Create wrapper function around httk steps
# simulate_css_d <- function(chem.cas, agelim_years, weight_category,
#                          samples, verbose = TRUE) {
#
#   if (verbose) {
#     cat(chem.cas,
#         paste0("(", paste(agelim_years, collapse = ", "), ")"),
#         weight_category,
#         "\n")
#   }
#
#   httkpop <- list(method = "vi",
#                   gendernum = NULL,
#                   agelim_years = agelim_years,
#                   agelim_months = NULL,
#                   weight_category = weight_category,
#                   reths = c(
#                     "Mexican American",
#                     "Other Hispanic",
#                     "Non-Hispanic White",
#                     "Non-Hispanic Black",
#                     "Other"
#                   ))
#
#   css <- try(
#     suppressWarnings({
#       mcs <- create_mc_samples(chem.cas = chem.cas,
#                                samples = samples,
#                                httkpop.generate.arg.list = httkpop,
#                                suppress.messages = TRUE)
#
#       calc_analytic_css(chem.cas = chem.cas,
#                         parameters = mcs,
#                         model = "3compartmentss",
#                         suppress.messages = TRUE)
#     }),
#     silent = TRUE
#   )
#
#   # Return
#   if (is(css, "try-error")) {
#     warning(paste0("simulate_css failed to generate data for CASN ", chem.cas))
#     list(NA)
#   } else {
#     list(css)
#   }
# }
#
# # Simulate `C_ss` values
# simulated_css <- map(casn, function(chem.cas) {
#   pop_demo |>
#     rowwise() |>
#     mutate(
#       css = simulate_css_d(.env$chem.cas, age_group, weight, .env$n_samples)
#     ) |>
#     ungroup()
# })
# simulated_css <- setNames(simulated_css, casn)
#
# # Remove CASN that failed simulate_css
# casn_keep <- map_lgl(simulated_css, function(df) {
#   !(length(df$css[[1]]) == 1 && is.na(df$css[[1]]))
# })
# simulated_css <- simulated_css[casn_keep]
#
# # Get median `C_ss` values for each age_group
# simulated_css <- map(
#   simulated_css,
#   function(cas_df) {
#     cas_df |>
#       nest(.by = age_group) |>
#       mutate(
#         age_median_css = map_dbl(data, function(df) median(unlist(df$css),
#                                                            na.rm = TRUE))
#       ) |>
#       unnest(data)
#   }
# )
#
# # Get median `C_ss` values for each weight
# simulated_css <- map(
#   simulated_css,
#   function(cas_df) {
#     cas_df |>
#       nest(.by = weight) |>
#       mutate(
#         weight_median_css = map_dbl(data, function(df) median(unlist(df$css),
#                                                               na.rm = TRUE))
#       ) |>
#       unnest(data) |>
#       arrange(age_min, weight)
#   }
# )
#
# #geo_tox_data$simulated_css <- simulated_css
#
# ################################################################################
#
# # 3. Determine discrepancies
#
# # This looks at the sets of css values and determines which are missing. This is
# # a population level analysis rather than a value by value analysis
#
# difference <- lapply(1:21, function(name) {lapply(1:22, function(ix) {
#   list(redo = setdiff(unlist(simulated_css[[name]]$css[[ix]]), unlist(geo_tox_data$simulated_css[[name]]$css[[ix]])),
#        geotox = setdiff(unlist(geo_tox_data$simulated_css[[name]]$css[[ix]]), unlist(simulated_css[[name]]$css[[ix]])))}
#   )})
#
#
# # Determine the sets of chemicals, weight class, and age groups where there are discrepancies
# indices <- lapply(difference, function(chem) {
#   lapply(1:22, function(ix) {which(chem[ix][[1]][[1]]>0)}
#   )})
