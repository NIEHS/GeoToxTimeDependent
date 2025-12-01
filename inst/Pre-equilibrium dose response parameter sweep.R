# Generate different exposure profiles
time_dr <- seq(from = 0, to = 30, by = 1/24)
time_2 <- seq(from = (30+1/24), to = 60, by = 1/24)

dose_a_low <- c(rep(0, 24), 1E-1*exp(-time_dr[1:697]))
dose_b_low <- 1E-1/19.15*abs(sin(time_dr))
dose_c_low <- rep(1E-2/3, 721)

acute_exposure_matrix_low <- matrix(c(time_dr, dose_a_low), ncol = 2, dimnames = list(NULL, c('time', 'dose')))
periodic_exposure_matrix_low <- matrix(c(time_dr, dose_b_low), ncol = 2, dimnames = list(NULL, c('time', 'dose')))
constant_exposure_matrix_low <- matrix(c(time_dr, dose_c_low), ncol = 2, dimnames = list(NULL, c('time', 'dose')))

dose_a_high <- c(rep(0, 24), 1E3*exp(-time_dr[1:697]))
dose_b_high <- 1E3/19.15*abs(sin(time_dr))
dose_c_high <- rep(1E2/3, 721)

acute_exposure_matrix_high <- matrix(c(time_dr, dose_a_high), ncol = 2, dimnames = list(NULL, c('time', 'dose')))
periodic_exposure_matrix_high <- matrix(c(time_dr, dose_b_high), ncol = 2, dimnames = list(NULL, c('time', 'dose')))
constant_exposure_matrix_high <- matrix(c(time_dr, dose_c_high), ncol = 2, dimnames = list(NULL, c('time', 'dose')))

dose_a_very_high <- c(rep(0, 24), 1E5*exp(-time_dr[1:697]))
dose_b_very_high <- 1E5/19.15*abs(sin(time_dr))
dose_c_very_high <- rep(1E4/3, 721)

acute_exposure_matrix_very_high <- matrix(c(time_dr, dose_a_very_high), ncol = 2, dimnames = list(NULL, c('time', 'dose')))
periodic_exposure_matrix_very_high <- matrix(c(time_dr, dose_b_very_high), ncol = 2, dimnames = list(NULL, c('time', 'dose')))
constant_exposure_matrix_very_high <- matrix(c(time_dr, dose_c_very_high), ncol = 2, dimnames = list(NULL, c('time', 'dose')))

# Set parameters
k_off <- -10:1
k_on <- 2:8
log_AC50 <- (-6:4)/2

log_kinetics <- expand.grid('k_off' = k_off,
                            'k_on' = k_on,
                            'logAC50' = log_AC50)
log_kinetics$KD <- log_kinetics$k_off - log_kinetics$k_on

# Set population and associated httk parameters
pop_simulator <- function(chem.cas = '',
                          samples = num_people){

  # Simulate populations of 10 year intervals, 'normal' weight category
  twenties_norm <- simulate_population_variability_parameters(chem.cas = chem.cas,
                                                              agelim_years = c(20, 29),
                                                              weight_category = 'Normal',
                                                              samples = samples,
                                                              set_seed = TRUE,
                                                              seed = 2345)
  thirties_norm <- simulate_population_variability_parameters(chem.cas = chem.cas,
                                                              agelim_years = c(30, 39),
                                                              weight_category = 'Normal',
                                                              samples = samples,
                                                              set_seed = TRUE,
                                                              seed = 2345)
  fourties_norm <- simulate_population_variability_parameters(chem.cas = chem.cas,
                                                              agelim_years = c(40, 49),
                                                              weight_category = 'Normal',
                                                              samples = samples,
                                                              set_seed = TRUE,
                                                              seed = 2345)
  fifties_norm <- simulate_population_variability_parameters(chem.cas = chem.cas,
                                                             agelim_years = c(50, 59),
                                                             weight_category = 'Normal',
                                                             samples = samples,
                                                             set_seed = TRUE,
                                                             seed = 2345)
  sixties_norm <- simulate_population_variability_parameters(chem.cas = chem.cas,
                                                             agelim_years = c(60, 69),
                                                             weight_category = 'Normal',
                                                             samples = samples,
                                                             set_seed = TRUE,
                                                             seed = 2345)
  seventies_norm <- simulate_population_variability_parameters(chem.cas = chem.cas,
                                                               agelim_years = c(70, 79),
                                                               weight_category = 'Normal',
                                                               samples = samples,
                                                               set_seed = TRUE,
                                                               seed = 2345)

  #simulate populations of 10 year intervals, 'obese' weight category

  twenties_ob <- simulate_population_variability_parameters(chem.cas = chem.cas,
                                                            agelim_years = c(20, 29),
                                                            weight_category = 'Obese',
                                                            samples = samples,
                                                            set_seed = TRUE,
                                                            seed = 2345)
  thirties_ob <- simulate_population_variability_parameters(chem.cas = chem.cas,
                                                            agelim_years = c(30, 39),
                                                            weight_category = 'Obese',
                                                            samples = samples,
                                                            set_seed = TRUE,
                                                            seed = 2345)
  fourties_ob <- simulate_population_variability_parameters(chem.cas = chem.cas,
                                                            agelim_years = c(40, 49),
                                                            weight_category = 'Obese',
                                                            samples = samples,
                                                            set_seed = TRUE,
                                                            seed = 2345)
  fifties_ob <- simulate_population_variability_parameters(chem.cas = chem.cas,
                                                           agelim_years = c(50, 59),
                                                           weight_category = 'Obese',
                                                           samples = samples,
                                                           set_seed = TRUE,
                                                           seed = 2345)
  sixties_ob <- simulate_population_variability_parameters(chem.cas = chem.cas,
                                                           agelim_years = c(60, 69),
                                                           weight_category = 'Obese',
                                                           samples = samples,
                                                           set_seed = TRUE,
                                                           seed = 2345)
  seventies_ob <- simulate_population_variability_parameters(chem.cas = chem.cas,
                                                             agelim_years = c(70, 79),
                                                             weight_category = 'Obese',
                                                             samples = samples,
                                                             set_seed = TRUE,
                                                             seed = 2345)
  return(list('normal' = list('twenties_norm' = twenties_norm,
                              'thirties_norm' = thirties_norm,
                              'fourties_norm' = fourties_norm,
                              'fifties_norm' = fifties_norm,
                              'sixties_norm' = sixties_norm,
                              'seventies_norm' = seventies_norm),
              'obese' = list('twenties_ob' = twenties_ob,
                             'thirties_ob' = thirties_ob,
                             'fourties_ob' = fourties_ob,
                             'fifties_ob' = fifties_ob,
                             'sixties_ob' = sixties_ob,
                             'seventies_ob' = seventies_ob)))
}

# Simulation parameters

simulate_parameters <- function(chem.cas = '',
                                mcs = list()){

  normal_pop <- mcs$normal
  obese_pop <- mcs$obese

  norm_20_param <- parametrize_httk(chem.cas = chem.cas, model = 'pbtk', mcs = normal_pop$twenties_norm)
  norm_30_param <- parametrize_httk(chem.cas = chem.cas, model = 'pbtk', mcs = normal_pop$thirties_norm)
  norm_40_param <- parametrize_httk(chem.cas = chem.cas, model = 'pbtk', mcs = normal_pop$fourties_norm)
  norm_50_param <- parametrize_httk(chem.cas = chem.cas, model = 'pbtk', mcs = normal_pop$fifties_norm)
  norm_60_param <- parametrize_httk(chem.cas = chem.cas, model = 'pbtk', mcs = normal_pop$sixties_norm)
  norm_70_param <- parametrize_httk(chem.cas = chem.cas, model = 'pbtk', mcs = normal_pop$seventies_norm)

  obese_20_param <- parametrize_httk(chem.cas = chem.cas, model = 'pbtk', mcs = obese_pop$twenties_ob)
  obese_30_param <- parametrize_httk(chem.cas = chem.cas, model = 'pbtk', mcs = obese_pop$thirties_ob)
  obese_40_param <- parametrize_httk(chem.cas = chem.cas, model = 'pbtk', mcs = obese_pop$fourties_ob)
  obese_50_param <- parametrize_httk(chem.cas = chem.cas, model = 'pbtk', mcs = obese_pop$fifties_ob)
  obese_60_param <- parametrize_httk(chem.cas = chem.cas, model = 'pbtk', mcs = obese_pop$sixties_ob)
  obese_70_param <- parametrize_httk(chem.cas = chem.cas, model = 'pbtk', mcs = obese_pop$seventies_ob)

  return(list('normal' = list('norm_20_param' = norm_20_param,
                              'norm_30_param' = norm_30_param,
                              'norm_40_param' = norm_40_param,
                              'norm_50_param' = norm_50_param,
                              'norm_60_param' = norm_60_param,
                              'norm_70_param' = norm_70_param),
              'obese' = list('obese_20_param' = obese_20_param,
                             'obese_30_param' = obese_30_param,
                             'obese_40_param' = obese_40_param,
                             'obese_50_param' = obese_50_param,
                             'obese_60_param' = obese_60_param,
                             'obese_70_param' = obese_70_param)))
}

# Generate parameters for Pentachlorophenol as an example
populations_test_pcp <- pop_simulator(chem.cas = '87-86-5', samples = 500)
parameters_test_pcp <- simulate_parameters(chem.cas = '87-86-5', mcs = populations_test_pcp)


# Simulate exposure
single_person_1_pcp <- solve_httk_model(chem.cas = '87-86-5', parameters = parameters_test_pcp$normal$norm_20_param[1,], plot = TRUE, data.matrix = rbind(periodic_exposure_matrix_high, matrix(c(time_2, rep(0, length(time_2))), ncol = 2, dimnames = list(NULL, c('time', 'dose')))), plot_average = TRUE, keep_all_columns = TRUE)
single_person_2_pcp <- solve_httk_model(chem.cas = '87-86-5', parameters = parameters_test_pcp$normal$norm_20_param[1,], plot = TRUE, data.matrix = rbind(acute_exposure_matrix_high, matrix(c(time_2, rep(0, length(time_2))), ncol = 2, dimnames = list(NULL, c('time', 'dose')))), plot_average = TRUE, keep_all_columns = TRUE)
single_person_3_pcp <- solve_httk_model(chem.cas = '87-86-5', parameters = parameters_test_pcp$normal$norm_20_param[1,], plot = TRUE, data.matrix = rbind(constant_exposure_matrix_high, matrix(c(time_2, rep(0, length(time_2))), ncol = 2, dimnames = list(NULL, c('time', 'dose')))), plot_average = TRUE, keep_all_columns = TRUE)
single_person_4_pcp <- solve_httk_model(chem.cas = '87-86-5', parameters = parameters_test_pcp$normal$norm_20_param[1,], plot = TRUE, data.matrix = rbind(periodic_exposure_matrix_very_high, matrix(c(time_2, rep(0, length(time_2))), ncol = 2, dimnames = list(NULL, c('time', 'dose')))), plot_average = TRUE, keep_all_columns = TRUE)
single_person_5_pcp <- solve_httk_model(chem.cas = '87-86-5', parameters = parameters_test_pcp$normal$norm_20_param[1,], plot = TRUE, data.matrix = rbind(acute_exposure_matrix_very_high, matrix(c(time_2, rep(0, length(time_2))), ncol = 2, dimnames = list(NULL, c('time', 'dose')))), plot_average = TRUE, keep_all_columns = TRUE)
single_person_6_pcp <- solve_httk_model(chem.cas = '87-86-5', parameters = parameters_test_pcp$normal$norm_20_param[1,], plot = TRUE, data.matrix = rbind(constant_exposure_matrix_very_high, matrix(c(time_2, rep(0, length(time_2))), ncol = 2, dimnames = list(NULL, c('time', 'dose')))), plot_average = TRUE, keep_all_columns = TRUE)
single_person_7_pcp <- solve_httk_model(chem.cas = '87-86-5', parameters = parameters_test_pcp$normal$norm_20_param[1,], plot = TRUE, data.matrix = rbind(periodic_exposure_matrix_low, matrix(c(time_2, rep(0, length(time_2))), ncol = 2, dimnames = list(NULL, c('time', 'dose')))), plot_average = TRUE, keep_all_columns = TRUE)
single_person_8_pcp <- solve_httk_model(chem.cas = '87-86-5', parameters = parameters_test_pcp$normal$norm_20_param[1,], plot = TRUE, data.matrix = rbind(acute_exposure_matrix_low, matrix(c(time_2, rep(0, length(time_2))), ncol = 2, dimnames = list(NULL, c('time', 'dose')))), plot_average = TRUE, keep_all_columns = TRUE)
single_person_9_pcp <- solve_httk_model(chem.cas = '87-86-5', parameters = parameters_test_pcp$normal$norm_20_param[1,], plot = TRUE, data.matrix = rbind(constant_exposure_matrix_low, matrix(c(time_2, rep(0, length(time_2))), ncol = 2, dimnames = list(NULL, c('time', 'dose')))), plot_average = TRUE, keep_all_columns = TRUE)

# Simulate dose-response
dose_response_pcp_plots_periodic <- list(dim(log_kinetics)[[1]])
dose_response_pcp_plots_acute <- list(dim(log_kinetics)[[1]])
dose_response_pcp_plots_constant <- list(dim(log_kinetics)[[1]])

for (i in 1:dim(log_kinetics)[[1]]){
  # Calculate response decay
  person <- single_person_7_pcp$numeric$`Person 1`
  max <- 100
  AC50 <- 10^log_kinetics$logAC50[[i]]
  n <- 1
  k_off <- 10^log_kinetics$k_off[[i]]
  k_on <- 10^log_kinetics$k_on[[i]]
  dose_response_pcp_plots_periodic[[i]] <- ggplot(response_decay_exponential(plasma_data = person, max = max, AC50 = AC50,
                                                                    n = n, k_off = k_off, k_on = k_on),
                                         aes(time, response)) + geom_point() + geom_line() +
    labs(title = paste('AC50 =', AC50, 'k_off =', k_off, 'k_on =', k_on, 'k_D =', k_off/k_on, 'Periodic'))
}

for (i in 1:dim(log_kinetics)[[1]]){
  # Calculate response decay
  person <- single_person_8_pcp$numeric$`Person 1`
  max <- 100
  AC50 <- 10^log_kinetics$logAC50[[i]]
  n <- 1
  k_off <- 10^log_kinetics$k_off[[i]]
  k_on <- 10^log_kinetics$k_on[[i]]
  dose_response_pcp_plots_acute[[i]] <- ggplot(response_decay_exponential(plasma_data = person, max = max, AC50 = AC50,
                                                                             n = n, k_off = k_off, k_on = k_on),
                                                  aes(time, response)) + geom_point() + geom_line() +
    labs(title = paste('AC50 =', AC50, 'k_off =', k_off, 'k_on =', k_on, 'k_D =', k_off/k_on, 'Acute'))
}

for (i in 1:dim(log_kinetics)[[1]]){
  # Calculate response decay
  person <- single_person_9_pcp$numeric$`Person 1`
  max <- 100
  AC50 <- 10^log_kinetics$logAC50[[i]]
  n <- 1
  k_off <- 10^log_kinetics$k_off[[i]]
  k_on <- 10^log_kinetics$k_on[[i]]
  dose_response_pcp_plots_constant[[i]] <- ggplot(response_decay_exponential(plasma_data = person, max = max, AC50 = AC50,
                                                                             n = n, k_off = k_off, k_on = k_on),
                                                  aes(time, response)) + geom_point() + geom_line() +
    labs(title = paste('AC50 =', AC50, 'k_off =', k_off, 'k_on =', k_on, 'k_D =', k_off/k_on, 'Constant'))
}
#ggplot(response_decay_exponential(plasma_data = single_person_7_pcp$numeric$`Person 1`, max = 100, AC50 = 1E-3, n = 1, k_off = 1E-9, k_on=1E2), aes(time, response)) + geom_point() + geom_line()
