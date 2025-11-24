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
k_off <- -5:1
k_on <- 2:8

log_kinetics <- expand.grid('k_off' = k_off,
                            'k_on' = k_on)
log_kinetics$KD <- log_kinetics$k_off - log_kinetics$k_on
