#' Generate steady state concentration (CSS)
#'
#' @param data_input A list with non-empty `dose_response` entry. This entry
#'   should contain a data.frame with a column `casn`.
#' @param n_samples Number of sampled CSS values to create.
#' @param load_all_httk_chemical_data Boolean determining whether to load
#'   additional `httk` chemical data.
#' @param set_seed Boolean for setting a seed.
#' @param seed Seed, used only if `set_seed` has value `TRUE`.
#'
#' @returns A data.frame with generated CSS values for a combination of age and
#'   weight parameters, including median values grouped by weight and by age.
#' @export
#'
#' @examplesIf FALSE
generate_steady_state_css <- function(data_input = NULL,
                                      n_samples = 500,
                                      load_all_httk_chemical_data = FALSE,
                                      set_seed = FALSE,
                                      seed = 42){
  # Check that data is provided
  if (is.null(data_input)){
    stop("Please provide data to produce steady-state plasma concentration!")
  }

  # Check that dose response data are available
  if (!("dose_response" %in% names(data_input))){
    stop("Please include dose response data in the 'dose_response' slot of `data_input`!")
  }

  # Check that n_samples is numeric and positive
  if (is.numeric(n_samples)) {
    n_samples <- as.integer(n_samples)
    if (n_samples <= 0) {
      warning('Parameter `n_samples` is non-positive, setting `n_samples` to 500!')
      n_samples <- 500
    }
  } else {
    stop('Parameter `n_samples` must be a positive integer!')
  }


  # Ensure httk data is loaded
  if (load_all_httk_chemical_data & !httk_data_loaded()){
    load_httk_data()
  }
  #httk::load_sipes2017()

  # Set seed if desired
  if (set_seed){
    set.seed(seed)
  }


  # Grab casn values from input data
  casn <- intersect(unique(data_input$dose_response$casn),
                    httk::get_cheminfo(suppress.messages = TRUE))

  # Define population demographics for httk simulation
  pop_demo <- dplyr::cross_join(
    tibble::tibble(age_group = list(c(0, 2), c(3, 5), c(6, 10), c(11, 15),
                                    c(16, 20), c(21, 30), c(31, 40), c(41, 50),
                                    c(51, 60), c(61, 70), c(71, 100))),
    tibble::tibble(weight = c('Normal', 'Obese'))) |>
    # Create column of lower age_group values
    dplyr::rowwise() |>
    dplyr::mutate(age_min = .data$age_group[1]) |>
    dplyr::ungroup()

  # Simulate `C_ss` values
  simulated_css <- purrr::map(casn, function(chem.cas) {
    pop_demo |>
      dplyr::rowwise() |>
      dplyr::mutate(
        css = simulate_css(chem.cas, .data$age_group, .data$weight, n_samples)
      ) |>
      dplyr::ungroup()
  })
  simulated_css <- stats::setNames(simulated_css, casn)

  # Remove CASN that failed simulate_css
  casn_keep <- purrr::map_lgl(simulated_css, function(df) {
    !(length(df$css[[1]]) == 1 && is.na(df$css[[1]]))
  })
  simulated_css <- simulated_css[casn_keep]
  print(simulated_css)

  # Get median `C_ss` values for each age_group
  simulated_css <- purrr::map(
    simulated_css,
    function(cas_df) {
      cas_df |>
        tidyr::nest(.by = .data$age_group) |>
        dplyr::mutate(
          age_median_css = purrr::map_dbl(.data$data, function(df) median(unlist(df$css),
                                                             na.rm = TRUE))
        ) |>
        tidyr::unnest(.data$data)
    }
  )
  print(simulated_css)
  #return(simulated_css)

  # Get median `C_ss` values for each weight
  simulated_css <- purrr::map(
    simulated_css,
    function(cas_df) {
      cas_df |>
        tidyr::nest(.by = .data$weight) |>
        dplyr::mutate(
          weight_median_css = purrr::map_dbl(.data$data, function(df) median(unlist(df$css),
                                                                na.rm = TRUE))
        ) |>
        tidyr::unnest(.data$data) |>
        dplyr::arrange(.data$age_min, .data$weight)
    }
  )

  return(simulated_css)
}


#' Generate population CSS
#'
#' @param chem.cas Chemical identifier CASN.
#' @param agelim_years Vector of length 2 denoting age range of population.
#' @param weight_category String with values `Normal` or `Obese`.
#' @param samples Number of CSS values to generate.
#' @param verbose Boolean determining whether to display progress report of
#'   internal function operations.
#'
#' @returns A list of CSS values based on input parameter values.
#' @export
#'
#' @examplesIf FALSE
simulate_css <- function(chem.cas, agelim_years, weight_category,
                         samples, verbose = TRUE) {

  if (verbose) {
    cat(chem.cas,
        paste0("(", paste(agelim_years, collapse = ", "), ")"),
        weight_category,
        "\n")
  }

  httkpop <- list(method = "vi",
                  gendernum = NULL,
                  agelim_years = agelim_years,
                  agelim_months = NULL,
                  weight_category = weight_category,
                  reths = c(
                    "Mexican American",
                    "Other Hispanic",
                    "Non-Hispanic White",
                    "Non-Hispanic Black",
                    "Other"
                  ))

  css <- try(
    suppressWarnings({
      mcs <- httk::create_mc_samples(chem.cas = chem.cas,
                                     samples = samples,
                                     httkpop.generate.arg.list = httkpop,
                                     suppress.messages = TRUE)

      httk::calc_analytic_css(chem.cas = chem.cas,
                              parameters = mcs,
                              model = "3compartmentss",
                              suppress.messages = TRUE)
    }),
    silent = TRUE
  )

  # Return
  if (is(css, "try-error")) {
    warning(paste0("simulate_css failed to generate data for CASN ", chem.cas))
    list(NA)
  } else {
    list(css)
  }
}

#' Generate population parameters for time-dependent `httk` models
#'
#' @param chem.cas Chemical identifier CASN.
#' @param agelim_years Vector of length 2 denoting age range of population.
#' @param age_matrix A matrix with columns `age` and `samples` to determine
#'   number of individuals to sample from each age.
#' @param weight_category String with values `Normal` or `Obese`.
#' @param samples Number of CSS values to generate.
#' @param model This can be 'pbtk', '3compartment', '3compartmentss', or
#'   '1compartment'. Default for this function is 'pbtk', See documentation for
#'   `httk::create_mc_samples` for more detailed information.
#' @param verbose Boolean determining whether to display progress report of
#'   internal function operations.
#' @param set_seed Boolean for setting a seed.
#' @param seed Seed, used only if `set_seed` has value `TRUE`.
#'
#' @returns A table of population parameters to use in time-dependent `httk`
#'   models. If `age_matrix` parameter is supplied, an additional column `age`
#'   is included in the output.
#' @export
#'
#' @examplesIf FALSE
simulate_population_variability_parameters <- function(chem.cas,
                                                       agelim_years = NULL,
                                                       age_matrix = NULL,
                                                       weight_category,
                                                       samples,
                                                       model = 'pbtk',
                                                       verbose = TRUE,
                                                       set_seed = FALSE,
                                                       seed = 42) {

  if (verbose) {
    if (!is.null(agelim_years)){
    cat(chem.cas,
        paste0("(", paste(agelim_years, collapse = ", "), ")"),
        weight_category,
        "\n")
    }
  }

  if (set_seed){
    set.seed(seed)
  }

  if (!is.null(age_matrix) & is.numeric(age_matrix)){
    age_col <- which(colnames(age_matrix) %in% 'age')
    sample_col <- which(colnames(age_matrix) %in% 'samples')
    valid_age_idx <- age_matrix[, age_col] > 0
    #print(age_col)
    #print(sample_col)
    #print(valid_age_idx)
    age_matrix <- age_matrix[valid_age_idx, ]
    population <- list()
    n <- length(age_matrix[, 1])
    valid_tries <- rep(FALSE, n)
    for (i in 1:length(age_matrix[, 1])){
      httkpop <- list(method = "vi",
                      gendernum = NULL,
                      agelim_years = age_matrix[i, age_col],
                      agelim_months = NULL,
                      weight_category = weight_category,
                      reths = c(
                        "Mexican American",
                        "Other Hispanic",
                        "Non-Hispanic White",
                        "Non-Hispanic Black",
                        "Other"
                      ))

      mcs <- try(
        suppressWarnings({
          httk::create_mc_samples(chem.cas = chem.cas,
                                  samples = age_matrix[i, sample_col],
                                  model = model,
                                  httkpop.generate.arg.list = httkpop,
                                  suppress.messages = TRUE)
        }),
        silent = TRUE
      )

      if (!is(mcs, 'try-error')){
        mcs$age <- age_matrix[i, age_col]
        population[[i]] <- mcs
        valid_tries[[i]] <- TRUE
      }
    }
    population <- population[valid_tries]
    return(data.table::rbindlist(population))
  }

  httkpop <- list(method = "vi",
                  gendernum = NULL,
                  agelim_years = agelim_years,
                  agelim_months = NULL,
                  weight_category = weight_category,
                  reths = c(
                    "Mexican American",
                    "Other Hispanic",
                    "Non-Hispanic White",
                    "Non-Hispanic Black",
                    "Other"
                  ))

  mcs <- try(
    suppressWarnings({
      httk::create_mc_samples(chem.cas = chem.cas,
                              samples = samples,
                              httkpop.generate.arg.list = httkpop,
                              suppress.messages = TRUE)
    }),
    silent = TRUE
  )

  # Return
  if (is(mcs, "try-error")) {
    warning(paste0("simulate_population_variability_parameters failed to generate data for CASN ", chem.cas))
    list(NA)
  } else {
    mcs <- cbind(mcs, data.frame(chem.cas = rep(chem.cas, samples)))
    return(mcs)
  }
}

#' Prepare parameters for a httk models
#'
#' @param chem.cas The chemical identifier CASN
#' @param model For which `httk` model to prepare parameters, from the list
#'   c('pbpk', '1comp', '3comp', 'gas_pbtk', 'fetal_pbtk', 'schmitt') with
#'   default `pbtk`.
#' @param mcs data.frame output from `httk::create_mc_samples`.
#'
#' @returns Model parameters for each row of the input `mcs` population
#'   parameters
#' @export
#'
#' @examplesIf FALSE
parametrize_httk <- function(chem.cas,
                             model = 'pbtk',
                             mcs = NULL){

  n_entries <- dim(mcs)[[1]]

  httk_model <- switch(
    model,
    'pbtk' = httk::parameterize_pbtk,
    '1comp' = httk::parameterize_1comp,
    '3comp' = httk::parameterize_3comp,
    'gas_pbtk' = httk::parameterize_gas_pbtk,
    'fetal_pbtk' = httk::parameterize_fetal_pbtk,
    httk::parameterize_pbtk
  )

  model_params <- as.data.frame(lapply(httk_model(chem.cas = chem.cas), rep, n_entries))

  col_intersect <- intersect(names(mcs), names(model_params))
  print(col_intersect)

  #print(mcs[1,])

  model_params <- model_params |> dplyr::select(-col_intersect)

  if ('age' %in% names(mcs)){
    col_intersect <- c(col_intersect, 'age')
  }

  model_params <- cbind(model_params, mcs |> dplyr::select(c(col_intersect, 'chem.cas')))


  #model_params[[col_intersect]] <- as.data.frame(mcs)[[col_intersect]]

  return(model_params)

}

#' Solve time-dependant httk model
#'
#' @param chem.cas The chemical identifier CASN.
#' @param model Which `httk` model to solve, from the list c('pbpk', '1comp',
#'   '3comp', 'gas_pbtk', 'fetal_pbtk', 'schmitt') with default `pbtk`.
#' @param parameters Set of parameters to supply to the specified httk model.
#'   Each row corresponds to a unique individual.
#' @param plot Boolean to determine whether to create a plot of blood plasma
#'   concentration vs time.
#' @param data.matrix Alternate parameter to supply time and dose information,
#'   with columns `time` and `dose`. If the `standarize` parameter is set to
#'   true, then the dose should be in units of \eqn{\frac{mg}{m^3}} and the
#'   `input.units` parameter will be automatically set to \eqn{\frac{mg}{kg}}.
#'   Otherwise, if the `standardize` parameter is set to FALSE, the `dose` units
#'   should match the `input.unit` parameter.
#' @param input.units Set of units for input dose. See documentation for `httk`
#'   models for more information.
#' @param output.units Set of units for output dose. See documentation for
#'   `httk` models for more information.
#' @param timestep Duration of each dose in days. If `NULL`, duration of doses
#'   is determined by using differences in dose times from `C_ext` matrix, with
#'   last dose determined by the mean value of dose duration intervals.
#' @param IR A single value or vector of inhalation rate values corresponding to
#'   the individuals given by the `parameters` data.
#' @param standardize Standardize exposure by individual-specific `IR` and
#'   bodyweight.
#' @param plot_average Whether to compute and plot the 'average' individual from
#'   the set of `parameters` representing the entire simulated population.
#'
#' @returns A list of httk model output values for each iteration and an
#'   optional `ggplot2` plot of the Cplasma vs time for each individual
#'   represented by the `parameters` argument.
#' @import ggplot2
#' @export
#'
#' @examplesIf FALSE
solve_httk_model <- function(chem.cas,
                             model = 'pbtk',
                             parameters = NULL,
                             plot = FALSE,
                             data.matrix = NULL,
                             input.units = 'mg/kg',
                             output.units = 'uM',
                             timestep = NULL,
                             IR = NULL,
                             standardize = TRUE,
                             plot_average = FALSE){

  n_entries <- dim(parameters)[[1]]

  httk_model <- switch(
    model,
    'pbtk' = httk::solve_pbtk,
    '1comp' = httk::solve_1comp,
    '3comp' = httk::solve_3comp,
    'gas_pbtk' = httk::solve_gas_pbtk,
    'fetal_pbtk' = httk::solve_fetal_pbtk,
    httk::solve_pbtk
  )

  results <- list()

  if (!is.null(data.matrix)) {
    col_index <- which(colnames(data.matrix) == 'time')
    days <- max(data.matrix[, col_index])

    if (!is.null(IR) & is.numeric(IR)){
      if (length(IR) > 1){


      IR_internal <- IR[IR > 0]
      if (length(IR_internal) < n_entries){
        k <- n_entries - length(IR_internal)
        IR_internal <- c(IR_internal, rep(1, k))
        warning('Length of valid entries in `IR` parameter is less than number of samples, filling remaining entries with value 1!')
      } else if (length(IR_internal) > n_entries){
        warning(paste('Length of valid entries in `IR` parameter is greater than number of samples, using only the first', n_entries, 'values!'))
        IR_internal <- IR_internal[1:n_entries]
      }
      } else if (length(IR) == 1){
        if (IR <= 0){
          stop('The `IR` parameter must be positive!')
        }
        IR_internal <- rep(IR, n_entries)
      }
    } else {
      IR_internal <- rep(1, n_entries)
      warning('Setting `IR` parameter to value 1 for each sample!')
    }

    #print(IR_internal)

    for (i in 1:n_entries){

      internal_dose <- data.matrix
      if (standardize){
        internal_dose <- calc_internal_dose_td(C_ext = data.matrix,
                                               IR = IR_internal[[i]], # Placeholder until this can be determined
                                               BW = parameters$BW[[i]],
                                               scaling = 1, # Placeholder until this can be determined
                                               timestep = timestep)
        input.units = 'mg/kg'
      }

      #print(paste('Iteration', i, 'internal dose', internal_dose[25:29,2]))
      #print(paste('inhalation rate', IR_internal[[i]]))
      temp <- as.data.frame(httk_model(chem.cas = chem.cas,
                                       parameters = parameters[i,],
                                       input.units = input.units,
                                       output.units = output.units,
                                       dosing.matrix = internal_dose,
                                       suppress.messages = TRUE,
                                       days = days))
      temp$iteration <- i
      temp$iteration <- as.factor(temp$iteration)
      temp$person <- 'general person'
      temp$person <- as.factor(temp$person)
      results[[i]] <- temp
    }
  } else {
  for (i in 1:n_entries){
    temp <- as.data.frame(httk_model(parameters = parameters[i,], suppress.messages = TRUE))
    temp$iteration <- i
    temp$iteration <- as.factor(temp$iteration)
    temp$person <- 'general person'
    temp$person <- as.factor(temp$person)
    results[[i]] <- temp
  }
  }

  if (plot & plot_average){
    average_person_param <- get_population_statistics(parameters)
    if(!is.null(data.matrix)){
      internal_dose <- data.matrix
      if (standardize){
        internal_dose <- calc_internal_dose_td(C_ext = data.matrix,
                                               IR = mean(IR_internal),
                                               BW = average_person_param$BW,
                                               scaling = 1,
                                               timestep = timestep)
      }
      average_person <- as.data.frame(httk_model(chem.cas = chem.cas,
                                                 parameters = average_person_param,
                                                 input.units = input.units,
                                                 output.units = output.units,
                                                 dosing.matrix = internal_dose,
                                                 suppress.messages = TRUE,
                                                 days = days))
      average_person$iteration <- NA
      average_person$iteration <- as.factor(average_person$iteration)
      average_person$person <- 'average_person'
      average_person$person <- as.factor(average_person$person)
    } else {
      average_person <- as.data.frame(httk_model(parameters = average_person_param,
                                                 suppress.messages = TRUE))
      average_person$iteration <- NA
      average_person$iteration <- as.factor(average_person$iteration)
      average_person$person <- 'average_person'
      average_person$person <- as.factor(average_person$person)
    }

  }

  names(results) <- paste('Person', 1:n_entries)

  if (plot) {
    test_plot <- ggplot() + geom_line(data = data.table::rbindlist(results), aes(.data$time, .data$Cplasma, color = .data$iteration, linetype = .data$person))
    if (plot_average){
      test_plot <- test_plot +  geom_line(data = average_person, aes(.data$time, .data$Cplasma, linetype = .data$person), color = 'black', linewidth =  0.5) +
        scale_linetype_manual(values = c('general person' = 'solid', 'average_person' = 'longdash')) +
        guides(linetype = guide_legend(title = 'Person type'))
      return(list(numeric = c(results, list('Average Person' = average_person)), plasma_plot = test_plot))
      }
    return(list(numeric = results, plasma_plot = test_plot))
  }

  return(results)
}

#' Calculate time-dependent internal chemical dose
#'
#' @param C_ext Matrix of external exposure, with columns given by `time` and
#'   `dose`. This is in units of \eqn{\frac{mg}{m^3}}
#' @param IR Inhalation rate in \eqn{\frac{m^3}{day}}.
#' @param BW Body weight in \eqn{kg}.
#' @param scaling scaling factor encompassing any required unit adjustments.
#' @param timestep Duration of each dose in days. If `NULL`, duration of doses
#'   is determined by using differences in dose times from `C_ext` matrix, with
#'   last dose determined by the mean value of dose duration intervals.
#'
#' @returns Matrix of internal chemical doses in \eqn{\frac{mg}{kg}}.
#' @export
#'
#' @examplesIf FALSE
calc_internal_dose_td <- function(C_ext,
                                  IR,
                                  BW,
                                  scaling = 1,
                                  timestep = NULL){
  if (!('matrix' %in% class(C_ext) & is.numeric(C_ext))){
    stop('C_ext must be a matrix that is numeric!')
  } else if (!all(c('time', 'dose') %in% colnames(C_ext))){
    stop('C_ext must contain columns `time` and `dose`!')
  }

  col_t <- which(colnames(C_ext) %in% 'time')
  col_d <- which(colnames(C_ext) %in% 'dose')

  exp_matrix <- C_ext[, c(col_t, col_d)]
  colnames(exp_matrix) <- c('time', 'dose')

  #print(exp_matrix)

  if (!is.null(timestep) & is.numeric(timestep)){
    if (timestep <= 0){
      stop('The `timestep` parameter must be positive!')
    }
    times <- rep(timestep, length(exp_matrix[, 1]))
  } else {
    time_intervals <- diff(exp_matrix[, 1])
    last_interval <- mean(time_intervals)
    times <- c(time_intervals, last_interval)
  }


  doses <- exp_matrix[, 2]

  #print(times)
  #print(doses)

  internal_doses <- doses * IR * times/BW * scaling

  D_int <- matrix(c(exp_matrix[, 1], internal_doses), ncol = 2)
  colnames(D_int) <- c('time', 'dose')

  return(D_int)

}

#' Get population statistics
#'
#' @param params Table of httk model parameters, output from
#' @param quantiles A vector of quantile values, from 0 to 1. If not supplied,
#'   the average of each column is returned.
#'
#' @returns Single row data.table of statistics of columns
#' @import data.table
#' @export
#'
#' @examplesIf FALSE
get_population_statistics <- function(params = NULL,
                                   quantiles = NULL){

  # check quantile and computed based on quantile values
  if (!is.null(quantiles)){
    if (any(quantiles > 1) | any(quantiles < 0)){
      warning('Removing values outside [0,1] from `quantile` parameter')
      quantiles <- quantiles[!((quantiles < 0) | (quantiles > 1))]
    }
    quantiles <- sort(quantiles)

    # Cast params as numeric and compute quantiles
    params <- as.data.table(params)[, lapply(.SD, as.numeric), .SDcols = 1:length(params)]
    return(params[, lapply(.SD, function(t) {stats::quantile(t, na.rm = TRUE, probs = quantiles)}), .SDcols = 1:length(params)])
  }


  # Cast columns as numeric and return column means
  return(as.data.frame(t(colMeans(data.table::as.data.table(params)[, lapply(.SD, as.numeric), .SDcols = 1:(length(params))], na.rm = TRUE))))
}

#' Sensitivity analysis for time-dependent external concentration
#'
#' @param chem.cas A vector containing the CASRN for a single chemical or
#'   several chemicals in a mixture.
#' @param C_ext A matrix with external chemical concentration time-series data.
#'   The column names must either be `time` and `dose` in the case of a single
#'   chemical or `time` and the CASRN for each chemical in a mixture of several
#'   chemicals.
#' @param model_parameters A data.frame of model parameters supplied by the
#'   functions `parametrize_httk()` representing individuals in a population.
#' @param model A character string indicating which `httk` model to use, from
#'   the options 'pbtk', '1comp', '3comp', 'gas_pbtk', 'fetal_pbtk'. The default
#'   is 'pbtk'.
#' @param hill_params A data.frame output from the function `fit_hill()`, with
#'   data corresponding to the chemical(s) contained in the `chem.cas`
#'   parameter.
#' @param IR The inhalation rate for the average person of the population given
#'   by the `model_parameters` input, in \eqn{\frac{m^3}{day}}.
#' @param max_mult The upper bound multiplier for max response. See
#'   `calc_concentration_response()` for more details.
#'
#' @returns A data.frame consisting of mixture response data from one or several
#'   chemicals. See `calc_concentration_response()` for more details. For a
#'   mixture of chemicals, this will be the set of common chemicals contained in
#'   the `chem.cas`, `C_ext`, and `hill_params` inputs.
#' @export
#'
#' @examplesIf FALSE
compute_c_ext_sensitivity <- function(chem.cas = NULL,
                                      C_ext = NULL,
                                      model_parameters = NULL,
                                      model = 'pbtk',
                                      hill_params = NULL,
                                      IR = NULL,
                                      max_mult = 1){

  if (any(lapply(c(C_ext, model_parameters, hill_params), is.null))){
    stop('User must input `C_ext`, `model_parameters`, and `hill_params` values!')
  }

  httk_model <- switch(
    model,
    'pbtk' = httk::solve_pbtk,
    '1comp' = httk::solve_1comp,
    '3comp' = httk::solve_3comp,
    'gas_pbtk' = httk::solve_gas_pbtk,
    'fetal_pbtk' = httk::solve_fetal_pbtk,
    httk::solve_pbtk
  )

  print(chem.cas)
  print(colnames(C_ext))
  print(hill_params$chem)


  col_index <- which(colnames(C_ext) == 'time')

  if (length(col_index) != 1){
    stop('The input matrix `C_ext` must have a single column named `time`!')
  }

  days <- max(C_ext[, col_index])

  average_person_param <- get_population_statistics(params = model_parameters)

  chem_names <- intersect(intersect(chem.cas, colnames(C_ext)), hill_params$chem)
  print(chem_names)

  C_invitro <- c()

  if (length(chem_names)){
    # Case where 'C_ext' consists of column `time` and a column for each
    # chemical concentration given by the associated CASRN
    for (i in 1:length(chem_names)){
      current_chemical_parameters <- model_parameters |> dplyr::filter(chem.cas == chem_names[[i]])
      print(dim(current_chemical_parameters))

      chem_col_index <- which(colnames(C_ext) == chem_names[[i]])
      C_ext_temp <- C_ext[, c(col_index, chem_col_index)]
      colnames(C_ext_temp) <- c('time', 'dose')

      internal_dose <- calc_internal_dose_td(C_ext = C_ext_temp,
                                             IR = IR,
                                             BW = average_person_param$BW,
                                             scaling = 1)

      average_person <- as.data.frame(httk_model(chem.cas = chem_names[[i]],
                                                 parameters = average_person_param,
                                                 dosing.matrix = internal_dose,
                                                 suppress.messages = TRUE,
                                                 days = days))

      max_Cplasma <- max(average_person$Cplasma)
      min_Cplasma <- min(average_person$Cplasma)

      #Cplasma_values <- seq(from = min_Cplasma, to = max_Cplasma, length = 1E3)

      C_invitro <- cbind(C_invitro, matrix(average_person$Cplasma, ncol = 1, dimnames = list(NULL, chem_names[[i]])))
    }
  } else if (all(colnames(C_ext) == c('time', 'dose'))) {
    internal_dose <- calc_internal_dose_td(C_ext = C_ext,
                                           IR = IR,
                                           BW = average_person_param$BW,
                                           scaling = 1)

    average_person <- as.data.frame(httk_model(chem.cas = chem.cas,
                                               parameters = average_person_param,
                                               dosing.matrix = internal_dose,
                                               suppress.messages = TRUE,
                                               days = days))

    max_Cplasma <- max(average_person$Cplasma)
    min_Cplasma <- min(average_person$Cplasma)

    #Cplasma_values <- seq(from = min_Cplasma, to = max_Cplasma, length = 1E3)

    C_invitro <- cbind(C_invitro, matrix(average_person$Cplasma, ncol = 1, dimnames = list(NULL, chem.cas)))
  } else {
    warning('Please check input parameters `chem.cas`, `C_ext`, and `hill_params`! Returning an empty list.')
    return(list())
  }



  resp <- calc_concentration_response(C_invitro = C_invitro,
                                      hill_params = hill_params,
                                      max_mult = max_mult,
                                      fixed = TRUE)


}


#' Time resolution sensitivity analysis for time-dependent external
#' concentration
#'
#' @param chem.cas A vector containing the CASRN for a single chemical or
#'   several chemicals in a mixture.
#' @param C_ext A matrix with external chemical concentration time-series data.
#'   The column names must either be `time` and `dose` in the case of a single
#'   chemical or `time` and the CASRN for each chemical in a mixture of several
#'   chemicals.
#' @param model_parameters A data.frame of model parameters supplied by the
#'   functions `parametrize_httk()` representing individuals in a population.
#' @param model A character string indicating which `httk` model to use, from
#'   the options 'pbtk', '1comp', '3comp', 'gas_pbtk', 'fetal_pbtk'. The default
#'   is 'pbtk'.
#' @param hill_params A data.frame output from the function `fit_hill()`, with
#'   data corresponding to the chemical(s) contained in the `chem.cas`
#'   parameter.
#' @param IR The inhalation rate for the average person of the population given
#'   by the `model_parameters` input, in \eqn{\frac{m^3}{day}}.
#' @param max_mult The upper bound multiplier for max response. See
#'   `calc_concentration_response()` for more details.
#' @param resolutions Time resolutions to subdivide the time-series data. These
#'   are positive numeric values, with units in days.
#' @param sample The values in `resolutions` divide the entire time-series data
#'   and `sample` determine which time point to use for each interval. The valid
#'   options are 'left', 'right', 'center', with a default value of 'left'.
#'
#' @returns A named list with entries corresponding to a data.frame consisting
#'   of mixture response data from one or several chemicals. See
#'   `compute_c_ext_sensitivity()` for more details on the structure of each
#'   entry. The names are the valid time resolutions given by the `resolutions`
#'   parameter.
#' @export
#'
#' @examplesIf FALSE
compute_time_resolution_sensitivity <- function(chem.cas = NULL,
                                                C_ext = NULL,
                                                model_parameters = NULL,
                                                model = 'pbtk',
                                                hill_params = NULL,
                                                IR = NULL,
                                                max_mult = 1,
                                                resolutions = NULL,
                                                sample = 'left'){

  if(is.null(resolutions) | !all(sapply(resolutions, is.numeric))){
    stop('Please input numeric values for parameter `resolution`!')
  }

  resolutions <- resolutions[resolutions > 0]

  if (!length(resolutions)){
    stop('Values for `resolutions` parameter must be positive!')
  }

  # Isolate the column of time values
  col_t <- which(colnames(C_ext) == 'time')
  times <- C_ext[, col_t]

  min_timestep <- min(diff(times))

  # Determine the time range
  start_time <- min(times)
  end_time <- max(times)

  # Remove any resolutions that exceed the max time. Stop if none exist.
  if (any(resolutions < end_time)){
    resolutions <- resolutions[resolutions < end_time]
  } else {
    stop('Values in `resolutions` must not exceed final timepoint of time-series!')
  }

  # Remove any resolutions that are below the minimum timestep. Stop if none exist.
  if (any(resolutions > min_timestep)){
    resolutions <- resolutions[resolutions > min_timestep]
  } else {
    stop('Values in `resolutions` must be greater than minimum timestep of time-series!')
  }


  # Determine the length of each time-series based on the resolution Take the
  # total elapsed time, divide by the resolution, take the ceiling to ensure
  # integral number of repetitions
  repetitions <- ceiling((end_time - start_time)/resolutions) + 1

  trials <- length(repetitions)

  trial_data <- list()

  for (i in 1:trials){
    temp_times <- seq(from = start_time, by = resolutions[[i]], length = repetitions[[i]])
    temp_times <- temp_times[temp_times <= end_time]


    # Sample values
    if (sample == 'left'){
      sample_times <- times[sapply(temp_times, function(t) {which.min(times < t)})]
      sample_index <- which(times %in% sample_times)
      dose_matrix <- cbind(matrix(temp_times, ncol = 1, dimnames = list(NULL, c('time'))), C_ext[sample_index, -col_t])
    } else if (sample == 'right') {
      sample_times <- sapply(c(temp_times[-1], end_time + 1), function(t) {utils::tail(times[times < t], 1)})
      sample_index <- which(times %in% sample_times)
      dose_matrix <- cbind(matrix(temp_times, ncol = 1, dimnames = list(NULL, c('time'))), C_ext[sample_index, -col_t])
    } else if (sample == 'center') {
      left_times <- times[sapply(temp_times, function(t) {which.min(times < t)})]
      left_indices <- which(times %in% left_times)
      right_times <- sapply(c(temp_times[-1], end_time + 1), function(t) {utils::tail(times[times < t], 1)})
      right_indices <- which(times %in% right_times)
      central_indices <- ceiling((right_indices + left_indices)/2)
      dose_matrix <- cbind(matrix(temp_times, ncol = 1, dimnames = list(NULL, c('time'))), C_ext[central_indices, -col_t])
    } else {
      warning("Invalid option selected for `sample`, defaulting to 'left'!")
      sample_times <- times[sapply(temp_times, function(t) {which.min(times < t)})]
      sample_index <- which(times %in% sample_times)
      dose_matrix <- cbind(matrix(temp_times, ncol = 1, dimnames = list(NULL, c('time'))), C_ext[sample_index, -col_t])
    }

    print(colnames(dose_matrix))
    print(colnames(C_ext))

    colnames(dose_matrix) <- c('time', colnames(C_ext)[-col_t])
    print(colnames(dose_matrix))




    temp_trial <- compute_c_ext_sensitivity(chem.cas = chem.cas,
                                            C_ext = dose_matrix,
                                            model_parameters = model_parameters,
                                            model = model,
                                            hill_params = hill_params,
                                            IR = IR,
                                            max_mult = max_mult)
    trial_data[[i]] <- temp_trial
  }

  names(trial_data) <- resolutions
  return(trial_data)
}

time_resolution_variation <- function(chem.cas,
                                      model = 'pbtk',
                                      parameters = NULL,
                                      plot = TRUE,
                                      C_ext = NULL,
                                      input.units = 'mg/kg',
                                      output.units = 'uM',
                                      timestep = NULL,
                                      IR = NULL,
                                      standardize = TRUE,
                                      plot_average = FALSE,
                                      resolutions = NULL,
                                      sample = 'left'){
  if(is.null(resolutions) | !all(sapply(resolutions, is.numeric))){
    stop('Please input numeric values for parameter `resolution`!')
  }

  resolutions <- resolutions[resolutions > 0]

  if (!length(resolutions)){
    stop('Values for `resolutions` parameter must be positive!')
  }

  # Isolate the column of time values
  col_t <- which(colnames(C_ext) == 'time')
  times <- C_ext[, col_t]

  min_timestep <- min(diff(times))

  # Determine the time range
  start_time <- min(times)
  end_time <- max(times)

  # Remove any resolutions that exceed the max time. Stop if none exist.
  if (any(resolutions < end_time)){
    resolutions <- resolutions[resolutions < end_time]
  } else {
    stop('Values in `resolution` must not exceed final timepoint of time-series!')
  }

  # Remove any resolutions that are below the minimum timestep. Stop if none exist.
  if (any(resolutions > min_timestep)){
    resolutions <- resolutions[resolutions > min_timestep]
  } else {
    stop('Values in `resolutions` must be greater than minimum timestep of time-series!')
  }

  # Determine the length of each time-series based on the resolution Take the
  # total elapsed time, divide by the resolution, take the ceiling to ensure
  # integral number of repetitions
  repetitions <- ceiling((end_time - start_time)/resolutions) + 1

  trials <- length(repetitions)

  trial_data <- list()

  for (i in 1:trials){
    temp_times <- seq(from = start_time, by = resolutions[[i]], length = repetitions[[i]])
    temp_times <- temp_times[temp_times <= end_time]


    # Sample values
    if (sample == 'left'){
      sample_times <- times[sapply(temp_times, function(t) {which.min(times < t)})]
      sample_index <- which(times %in% sample_times)
      dose_matrix <- cbind(matrix(temp_times, ncol = 1, dimnames = list(NULL, c('time'))), C_ext[sample_index, -col_t])
    } else if (sample == 'right') {
      sample_times <- sapply(c(temp_times[-1], end_time + 1), function(t) {utils::tail(times[times < t], 1)})
      sample_index <- which(times %in% sample_times)
      dose_matrix <- cbind(matrix(temp_times, ncol = 1, dimnames = list(NULL, c('time'))), C_ext[sample_index, -col_t])
    } else if (sample == 'center') {
      left_times <- times[sapply(temp_times, function(t) {which.min(times < t)})]
      left_indices <- which(times %in% left_times)
      right_times <- sapply(c(temp_times[-1], end_time + 1), function(t) {utils::tail(times[times < t], 1)})
      right_indices <- which(times %in% right_times)
      central_indices <- ceiling((right_indices + left_indices)/2)
      dose_matrix <- cbind(matrix(temp_times, ncol = 1, dimnames = list(NULL, c('time'))), C_ext[central_indices, -col_t])
    } else {
      warning("Invalid option selected for `sample`, defaulting to 'left'!")
      sample_times <- times[sapply(temp_times, function(t) {which.min(times < t)})]
      sample_index <- which(times %in% sample_times)
      dose_matrix <- cbind(matrix(temp_times, ncol = 1, dimnames = list(NULL, c('time'))), C_ext[sample_index, -col_t])
    }

    print(colnames(dose_matrix))
    print(colnames(C_ext))

    colnames(dose_matrix) <- c('time', colnames(C_ext)[-col_t])
    print(colnames(dose_matrix))

    temp_trial <- solve_httk_model(chem.cas = chem.cas,
                                   model = model,
                                   parameters = parameters,
                                   plot = plot,
                                   data.matrix = dose_matrix,
                                   input.units = input.units,
                                   output.units = output.units,
                                   timestep = timestep,
                                   IR = IR,
                                   standardize = standardize,
                                   plot_average = plot_average)



    trial_data[[i]] <- temp_trial
  }

  names(trial_data) <- resolutions
  return(trial_data)

}


#' Combine dosing data from multiple pathways
#'
#' @param dosing_matrices List of matrices, each representing a time-series of
#'   chemical exposure.
#'
#' @returns A matrix with two columns, `time` and `dose`, which totals all the
#'   chemical exposures across all time points in the time-series exposure data.
#' @export
#'
#' @examplesIf FALSE
#' a <- matrix(c(seq(from = 0, to = 1, by = .1), 20*exp(-1*0:10)),
#'  ncol = 2, dimnames = list(NULL, c('time', 'dose')))
#' b <- matrix(c(seq(from = 0, to = 1, by = 1/6), 1:7),
#' ncol = 2, dimnames = list(NULL, c('time', 'dose')))
#' c <- combine_dosing_pathways(dosing_matrices = list(a,b))
#' c
combine_dosing_pathways <- function(dosing_matrices = NULL){
  # Check that input is not NULL
  if (is.null(dosing_matrices)){
    stop('Must input a list of matrices!')
  }

  # Cast input as list
  if (!is.list(dosing_matrices)){
    dosing_matrices <- list(dosing_matrices)
  }

  matrix_indices <- which(sapply(dosing_matrices, function(t) {'matrix' %in% class(t)}))


  times <- c()
  for (i in matrix_indices){
    temp_matrix <- dosing_matrices[[i]]
    time_col <- which(colnames(temp_matrix) == 'time')
    times <- unique(c(times, temp_matrix[, time_col]))
  }

  times <- sort(times)


  all_doses <- matrix(times, ncol = 1, dimnames = list(NULL, c('time')))

  doses <- rep(0, length(times))

  for (i in matrix_indices){
    temp_matrix <- dosing_matrices[[i]]
    time_col <- which(colnames(temp_matrix) == 'time')
    dose_col <- which(colnames(temp_matrix) == 'dose')
    dosing_index <- which(times %in% temp_matrix[, time_col])
    doses[dosing_index] <- temp_matrix[, dose_col]
    new_doses <- matrix(doses, ncol = 1, dimnames = list(NULL, paste0('route_', i)))
    all_doses <- cbind(all_doses, new_doses)
    doses <- rep(0, length(times))
  }

  row_sums <- apply(all_doses[, -1], 1, sum)

  final_dose <- matrix(c(times, row_sums), ncol = 2, dimnames = list(NULL, c('dose', 'time')))
  return(final_dose)
}



#' Load `httk` data
#'
#' This function loads in extended `httk` data and toggles the option
#' `httk_data_loaded` to `TRUE`.
#' @returns NULL
#' @export
#'
#' @examplesIf FALSE
load_httk_data <- function(){
  httk::load_sipes2017()
  httk::load_pradeep2020()
  httk::load_dawson2021()
  httk::load_honda2023()
  set_GeoToxTimeDependent_option('httk_data_loaded' = TRUE)
}

#' Reset `httk` data
#'
#' Helper function that resets loaded `httk` data and toggles the option
#' `httk_data_loaded` to `FALSE`.
#'
#' @returns NULL
#' @export
#'
#' @examplesIf FALSE
#'
reset_httk_data <- function(){
  httk::reset_httk()
  set_GeoToxTimeDependent_option('httk_data_loaded' = FALSE)
}
