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

  exceedance_data <- data.frame('threshold' = names(thresholds),
                                'time_exceeded' = time_exceeded,
                                'response_sum' = response_sum)
  return(exceedance_data)


}
