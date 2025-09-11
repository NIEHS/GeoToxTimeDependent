#' Exponential decay of response in time-varying dose-response
#'
#' This function models time-varying dose-response given a time series of blood
#' plasma concentration, parameters associated to the dose-response curve, a
#' response decay rate, and optional weighting parameters for balancing the
#' decaying plasma levels with new plasma levels.
#'
#' For a given dose-response relationship, let \eqn{r = f(d)} denote the response
#' level corresponding to a dose, with \eqn{f} the dose-response function. In this
#' framework, we combine the response from the previous time step with the new
#' dose to determine the next response level.
#'
#' Let \eqn{r_t} denote the response at time \eqn{t}, \eqn{d_{t+\Delta t}} the
#' new plasma level at time \eqn{t + \Delta t}, and \eqn{r_{t+\Delta t}} the
#' response at time \eqn{t + \Delta t} and \eqn{t_{\frac{1}{2}}} the half-life
#' for the response decay. The decay at time \eqn{t + \Delta t} from the
#' response at time \eqn{t} is given by
#' \eqn{r^{d}_{t + \Delta t} = r_t\times\exp(-\frac{\ln(2)\times\Delta
#' t}{t_{\frac{1}{2}}})} and the corresponding dose is given by
#' \eqn{f^{-1}(r^d_{t + \Delta t})}. This is combined with the new plasma
#' reading, \eqn{d_{t + \Delta t}} as
#' \eqn{a\times d_{t + \Delta t} + b\times f^{-1}(r^d_{t + \Delta t})} and input
#' into \eqn{f} to determine the new response,
#' \eqn{r_{t + \Delta t} = f(a\times d_{t + \Delta t} + b\times f^{-1}(r^d_{t + \Delta t}))}.
#'
#' @param plasma_data A data.frame containing dose-response data. Must contain
#'   columns `time` and `Cplasma`.
#' @param max The max response level for assay.
#' @param AC50 AC50 level for assay.
#' @param n The slope factor for the assay response.
#' @param t_half Half-life in days of response decay.
#' @param a Control parameter for new plasma level.
#' @param b Control parameter for decaying plasma level.
#'
#' @returns A data.frame with a single column, `response`, that gives the
#'   response due to new plasma readings and decaying plasma levels.
#' @export
#'
#' @examplesIf FALSE
response_decay <- function(plasma_data = NULL,
                           max = NULL,
                           AC50 = NULL,
                           n = NULL,
                           t_half = NULL,
                           a = 1,
                           b = 1){
  if (!(all(c('time', 'Cplasma') %in% names(plasma_data)))){
    stop('`plasma_data` must include columns for "time" and "Cplasma"!')
  }

  if(any(sapply(c(max, AC50, n, t_half), is.null))){
    stop('Missing input value for one of `max`, `AC50`, `n`, `t_half`!')
  }

  if (!is.numeric(a) | a < 0){
    warning('The value for `a` must be positive! Setting to default value 1.')
    a <- 1
  }

  if (!is.numeric(b) | b < 0){
    warning('The value for `b` must be positive! Setting to default value 1.')
    b <- 1
  }

  t_idx <- which(names(plasma_data) == 'time')
  p_idx <- which(names(plasma_data) == 'Cplasma')

  times <- plasma_data[, t_idx]
  time_diffs <- diff(times)
  plasma <- plasma_data[, p_idx]

  response <- numeric(length(times))

  response[[1]] <- hill_val(conc = plasma[[1]], max = max, AC50 = AC50, n = n)

  for (i in 1:length(time_diffs)){
    p_new <- plasma[[1 + i]]
    p_old <- hill_conc(resp = response[[i]]*exp(-log(2)*time_diffs[[i]]/t_half), max = max, AC50 = AC50, n = n)
    p_combined <- a*p_new + b*p_old
    r_new <- hill_val(conc = p_combined, max = max, AC50 = AC50, n = n)
    response[[i + 1]] <- r_new
  }

  return(data.frame('time' = times,
                    'response' = response))
}

response_decay_exponential <- function(plasma_data = NULL,
                                       max = NULL,
                                       AC50 = NULL,
                                       n = NULL,
                                       K_decay = NULL,
                                       k_off = NULL,
                                       k_on = NULL){
  if (!(all(c('time', 'Cplasma') %in% names(plasma_data)))){
    stop('`plasma_data` must include columns for "time" and "Cplasma"!')
  }

  if(any(sapply(c(max, AC50, n), is.null))){
    stop('Missing input value for one of `max`, `AC50`, `n`!')
  }

  dose_dependent <- FALSE

  if (is.null(K_decay)){
    if (is.null(k_off) | is.null(k_on)){
    stop('Either `K_decay` or both `k_on` and `k_off` must be supplied!')
    } else {
      dose_dependent <- TRUE
    }
  } else if (!is.null(K_decay)){
    warning('Defaulting to decay constant that is not dose-dependent!')
  }

  t_idx <- which(names(plasma_data) == 'time')
  p_idx <- which(names(plasma_data) == 'Cplasma')

  times <- plasma_data[, t_idx]
  time_diffs <- diff(times)
  plasma <- plasma_data[, p_idx]


  response <- numeric(length(times))

  response[[1]] <- 0

  if (dose_dependent){
    for (i in 1:length(time_diffs)){
      current_response <- response[[i]]
      response_addition <- hill_val(conc = plasma[[i]], max = max, AC50 = AC50, n = n)
      decay_coefficient <- k_off + k_on*plasma[[i]]
      current_decay <- exp(-decay_coefficient*time_diffs[[i]])
      new_response <- response_addition + (current_response - response_addition)*current_decay
      response[[i+1]] <- new_response
    }
  } else {
    for (i in 1:length(time_diffs)){
      current_response <- response[[i]]
      response_addition <- hill_val(conc = plasma[[i]], max = max, AC50 = AC50, n = n)
      current_decay <- exp(-K_decay*time_diffs[[i]])
      new_response <- response_addition + (current_response - response_addition)*current_decay
      response[[i+1]] <- new_response
    }
  }

  return(data.frame('time' = times,
                    'response' = response))
}

#' ODE version of Response Decay
#'
#' Models time-dependent assay response given plasma concentration following the
#' ODE system \eqn{\frac{dR}{dt} = -KR + Kf(d(t))} where \eqn{R} is the
#' response, \eqn{d} is the dose function, \eqn{f} is the dose-response
#' function, and \eqn{K} is a decay rate coefficient. For constant dose regime
#' \eqn{d}, the steady-state response is given by \eqn{R = f(d)}. The value of
#' the decay rate constant \eqn{K} is determined by \eqn{K =
#' -\frac{\log(\text{tol})}{t_{ss}}} where \eqn{\text{tol}} is the tolerance
#' level and \eqn{t_{ss}} is the time it takes the response to converge to
#' \eqn{R_{ss}} under constant dosing regime within a specified tolerance level.
#'  See the description of `tol` below for more details.
#'
#' @param plasma_data A data.frame containing dose-response data. Must contain
#'   columns `time` and `Cplasma`.
#' @param max The max response level for assay.
#' @param min The min response level for assay.
#' @param AC50 AC50 level for assay.
#' @param n The slope factor for the assay response.
#' @param t_ss Time to steady-state of assay response within tolerance level
#'   `tol`.
#' @param tol A value strictly between 0 and 1, giving the level determining how
#'   close the response level is to steady-state value. More specifically, let
#'   \eqn{R_{ss}} is the steady-state response corresponding to a specific
#'   constant dose level and \eqn{t_{ss}} the time it takes to achieve
#'   \eqn{|\frac{R(t)}{R_{ss}}| < 1 - \text{tol}} for \eqn{t > t_{ss}}.
#'
#' @returns A `deSolve` object with columns 'time', 'response', and
#'   'new_response'. The column 'new_response' is the incremental response value
#'   from the dose at the time point while 'response' gives the overall
#'   response.
#' @export
#'
#' @examplesIf FALSE
response_decay_ode <- function(plasma_data = NULL,
                               max = NULL,
                               min = NULL,
                               AC50 = NULL,
                               n = NULL,
                               t_ss = NULL,
                               tol = 1E-4){
  if (!(all(c('time', 'Cplasma') %in% names(plasma_data)))){
    stop('`plasma_data` must include columns for "time" and "Cplasma"!')
  }

  if(any(sapply(c(max, min, AC50, n, t_ss), is.null))){
    stop('Missing input value for one of `max`, `min`, `AC50`, `n`, `t_ss`!')
  }

  t_idx <- which(names(plasma_data) == 'time')
  p_idx <- which(names(plasma_data) == 'Cplasma')

  times <- plasma_data[, t_idx]
  plasma <- plasma_data[, p_idx]
  response <- hill_val(conc = plasma, max = max, AC50 = AC50, n = n)
  # Normalize response to between 0 and 1
  normalized_response <- (response - min)/(max - min)

  dose_table <- data.frame(time = times,
                           plasma = plasma,
                           response = response,
                           normalized_response = normalized_response)



  input_ode_response <- stats::approxfun(dose_table[, -c(2,4)], rule = 2)


  response_ode <- function(t, x, parms){
    with(as.list(c(parms, x)),{
      response <- input_ode_response(t)
      dR <- -K*R + K*response
      list(dR, new_reponse = response)
    })
  }

  parms <- c(K = -log(tol)/t_ss)
  xstart <- c(R = hill_val(conc = 1E-10, max = max, AC50 = AC50, n=n))

  ode_response_solution <- deSolve::ode(y = xstart, times = times, func = response_ode, parms)

  # Undo normalization
  #ode_response_solution[, 2] <- (max-min)*ode_response_solution[, 2] + min

  colnames(ode_response_solution) <- c('time', 'response', 'new_response')

  return(ode_response_solution)

}

#' Determine the time and AUC of excess response
#'
#' Given a response curve and response threshold(s), it is possible to determine
#' for how long the response exceeds a particular threshold and the AUC of the
#' excess response above that threshold. This function takes in a named vector
#' or list of threshold values and a time-series of response data output from
#' the function [response_decay()].
#'
#' @param thresholds A named list or vector of threshold names and response
#'   threshold values.
#' @param response_data A time-series of response data, with columns `time` and
#'   `response`.
#'
#' @returns A data.frame consisting of a column of valid threshold names, the
#'   duration of time the response exceeded the threshold values, and the AUC
#'   between the threshold value and response curve.
#' @export
#'
#' @examplesIf FALSE
threshold_exceedance <- function(thresholds,
                                 response_data){
  if (is.null(thresholds) | !(is.numeric(thresholds) | is.list(thresholds))){
    stop('The argument `thresholds` must be a named list or vector of numeric values!')
  }

  numeric_idx <- which(unname(sapply(thresholds, is.numeric)))

  if (length(numeric_idx)){
    thresholds <- thresholds[numeric_idx]
  } else {
    stop('The parameter `thresholds` contains no numeric values!')
  }

  if (!(all(c('time', 'response') %in% names(response_data)))){
    stop('`response` must include columns for "time" and "response"!')
  }

  t_idx <- which(names(response_data) == 'time')
  r_idx <- which(names(response_data) == 'response')

  times <- response_data[, t_idx]
  time_diffs <- diff(times)

  response <- response_data[1:(length(times)-1), r_idx]

  time_exceeded <- numeric(length(numeric_idx))
  response_sum <- numeric(length(numeric_idx))

  for (i in 1:length(numeric_idx)){
    threshold_level <- thresholds[[i]]
    time_exceeded[[i]] <- sum(time_diffs[which(threshold_level <= response)])
    response_sum[[i]] <- sum((response[which(threshold_level <= response)] - threshold_level)*time_diffs[which(threshold_level <= response)])
  }

  exceedance_data <- data.frame('threshold_name' = names(thresholds),
                                'threshold_value' = unlist(unname(thresholds)),
                                'time_exceeded' = time_exceeded,
                                'response_sum' = response_sum)
  return(exceedance_data)


}
