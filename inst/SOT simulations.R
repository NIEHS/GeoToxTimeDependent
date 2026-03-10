# Load libraries
library(httk)
library(dplyr)
library(tidyr)
library(data.table)
#library(GeoToxTimeDependent)
#library(gganimate)
#library(gifski)
library(ggridges)
library(cowplot)
#devtools::load_all()

#load_httk_data()

# Generate exposure scenarios. Same three across all chemicals, populations
time <- seq(from = 0, to = 30, by = 1/24)
dose_a <- c(rep(0, 24), 100*exp(-time[1:697]))
dose_b <- 100/19.15*abs(sin(time))
dose_c <- rep(10/3, length(time))

acute_exposure_matrix <- matrix(c(time, dose_a), ncol = 2, dimnames = list(NULL, c('time', 'dose')))
periodic_exposure_matrix <- matrix(c(time, dose_b), ncol = 2, dimnames = list(NULL, c('time', 'dose')))
constant_exposure_matrix <- matrix(c(time, dose_c), ncol = 2, dimnames = list(NULL, c('time', 'dose')))

# List of chemicals
chemical_casns <- c("53-96-3", "95-95-4", "88-06-2", "91-94-1", "119-90-4",
                    "60-11-7", "100-02-7", "101-14-4", "101-77-9", "92-87-5",
                    "133-06-2", "63-25-2", "120-80-9", "510-15-6", "84-74-2",
                    "123-31-9", "72-43-5", "106-50-3", "56-38-2", "87-86-5",
                    "95-80-7", "106-46-7", "107-21-1", "111-44-4", "1120-71-4",
                    "114-26-1", "117-81-7", "121-14-2", "121-69-7", "131-11-3",
                    "132-64-9", "133-90-4", "1582-09-8", "51-28-5", "51-79-6",
                    "532-27-4", "534-52-1", "57-74-9", "584-84-9", "64-67-5",
                    "76-44-8", "77-47-4", "77-78-1", "79-44-7", "82-68-8",
                    "822-06-0", "91-22-5", "92-52-4", "94-75-7", "96-09-3",
                    "98-86-2")


num_people <- 500

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

# Exposure simualtions

acute_exposure <- function(chem.cas = '',
                           sim_parms = list(),
                           acute_matrix = matrix()){

  norm_parms <- sim_parms$normal
  obese_parms <- sim_parms$obese

  acute_norm_20 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_20_param, plot = TRUE, data.matrix = acute_matrix, plot_average = TRUE)
  acute_norm_30 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_30_param, plot = TRUE, data.matrix = acute_matrix, plot_average = TRUE)
  acute_norm_40 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_40_param, plot = TRUE, data.matrix = acute_matrix, plot_average = TRUE)
  acute_norm_50 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_50_param, plot = TRUE, data.matrix = acute_matrix, plot_average = TRUE)
  acute_norm_60 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_60_param, plot = TRUE, data.matrix = acute_matrix, plot_average = TRUE)
  acute_norm_70 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_70_param, plot = TRUE, data.matrix = acute_matrix, plot_average = TRUE)

  acute_obese_20 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_20_param, plot = TRUE, data.matrix = acute_matrix, plot_average = TRUE)
  acute_obese_30 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_30_param, plot = TRUE, data.matrix = acute_matrix, plot_average = TRUE)
  acute_obese_40 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_40_param, plot = TRUE, data.matrix = acute_matrix, plot_average = TRUE)
  acute_obese_50 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_50_param, plot = TRUE, data.matrix = acute_matrix, plot_average = TRUE)
  acute_obese_60 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_60_param, plot = TRUE, data.matrix = acute_matrix, plot_average = TRUE)
  acute_obese_70 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_70_param, plot = TRUE, data.matrix = acute_matrix, plot_average = TRUE)

  return(list('normal' = list('acute_norm_20' = acute_norm_20,
                              'acute_norm_30' = acute_norm_30,
                              'acute_norm_40' = acute_norm_40,
                              'acute_norm_50' = acute_norm_50,
                              'acute_norm_60' = acute_norm_60,
                              'acute_norm_70' = acute_norm_70),
              'obese' = list('acute_obese_20' = acute_obese_20,
                             'acute_obese_30' = acute_obese_30,
                             'acute_obese_40' = acute_obese_40,
                             'acute_obese_50' = acute_obese_50,
                             'acute_obese_60' = acute_obese_60,
                             'acute_obese_70' = acute_obese_70)))
}

periodic_exposure <- function(chem.cas = '',
                              sim_parms = list(),
                              periodic_matrix = matrix()){

  norm_parms <- sim_parms$normal
  obese_parms <- sim_parms$obese

  periodic_norm_20 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_20_param, plot = TRUE, data.matrix = periodic_matrix, plot_average = TRUE)
  periodic_norm_30 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_30_param, plot = TRUE, data.matrix = periodic_matrix, plot_average = TRUE)
  periodic_norm_40 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_40_param, plot = TRUE, data.matrix = periodic_matrix, plot_average = TRUE)
  periodic_norm_50 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_50_param, plot = TRUE, data.matrix = periodic_matrix, plot_average = TRUE)
  periodic_norm_60 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_60_param, plot = TRUE, data.matrix = periodic_matrix, plot_average = TRUE)
  periodic_norm_70 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_70_param, plot = TRUE, data.matrix = periodic_matrix, plot_average = TRUE)

  periodic_obese_20 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_20_param, plot = TRUE, data.matrix = periodic_matrix, plot_average = TRUE)
  periodic_obese_30 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_30_param, plot = TRUE, data.matrix = periodic_matrix, plot_average = TRUE)
  periodic_obese_40 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_40_param, plot = TRUE, data.matrix = periodic_matrix, plot_average = TRUE)
  periodic_obese_50 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_50_param, plot = TRUE, data.matrix = periodic_matrix, plot_average = TRUE)
  periodic_obese_60 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_60_param, plot = TRUE, data.matrix = periodic_matrix, plot_average = TRUE)
  periodic_obese_70 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_70_param, plot = TRUE, data.matrix = periodic_matrix, plot_average = TRUE)

  return(list('normal' = list('periodic_norm_20' = periodic_norm_20,
                              'periodic_norm_30' = periodic_norm_30,
                              'periodic_norm_40' = periodic_norm_40,
                              'periodic_norm_50' = periodic_norm_50,
                              'periodic_norm_60' = periodic_norm_60,
                              'periodic_norm_70' = periodic_norm_70),
              'obese' = list('periodic_obese_20' = periodic_obese_20,
                             'periodic_obese_30' = periodic_obese_30,
                             'periodic_obese_40' = periodic_obese_40,
                             'periodic_obese_50' = periodic_obese_50,
                             'periodic_obese_60' = periodic_obese_60,
                             'periodic_obese_70' = periodic_obese_70)))
}

constant_exposure <- function(chem.cas = '',
                              sim_parms = list(),
                              constant_matrix = matrix()){

  norm_parms <- sim_parms$normal
  obese_parms <- sim_parms$obese

  constant_norm_20 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_20_param, plot = TRUE, data.matrix = constant_matrix, plot_average = TRUE)
  constant_norm_30 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_30_param, plot = TRUE, data.matrix = constant_matrix, plot_average = TRUE)
  constant_norm_40 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_40_param, plot = TRUE, data.matrix = constant_matrix, plot_average = TRUE)
  constant_norm_50 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_50_param, plot = TRUE, data.matrix = constant_matrix, plot_average = TRUE)
  constant_norm_60 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_60_param, plot = TRUE, data.matrix = constant_matrix, plot_average = TRUE)
  constant_norm_70 <- solve_httk_model(chem.cas = chem.cas, parameters = norm_parms$norm_70_param, plot = TRUE, data.matrix = constant_matrix, plot_average = TRUE)

  constant_obese_20 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_20_param, plot = TRUE, data.matrix = constant_matrix, plot_average = TRUE)
  constant_obese_30 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_30_param, plot = TRUE, data.matrix = constant_matrix, plot_average = TRUE)
  constant_obese_40 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_40_param, plot = TRUE, data.matrix = constant_matrix, plot_average = TRUE)
  constant_obese_50 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_50_param, plot = TRUE, data.matrix = constant_matrix, plot_average = TRUE)
  constant_obese_60 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_60_param, plot = TRUE, data.matrix = constant_matrix, plot_average = TRUE)
  constant_obese_70 <- solve_httk_model(chem.cas = chem.cas, parameters = obese_parms$obese_70_param, plot = TRUE, data.matrix = constant_matrix, plot_average = TRUE)

  return(list('normal' = list('constant_norm_20' = constant_norm_20,
                              'constant_norm_30' = constant_norm_30,
                              'constant_norm_40' = constant_norm_40,
                              'constant_norm_50' = constant_norm_50,
                              'constant_norm_60' = constant_norm_60,
                              'constant_norm_70' = constant_norm_70),
              'obese' = list('constant_obese_20' = constant_obese_20,
                             'constant_obese_30' = constant_obese_30,
                             'constant_obese_40' = constant_obese_40,
                             'constant_obese_50' = constant_obese_50,
                             'constant_obese_60' = constant_obese_60,
                             'constant_obese_70' = constant_obese_70)))
}

# AUC and Cplasma distributions

httk_distributions <- function(exposure_sims = list(),
                               exposure_params = list(),
                               Scenario = '',
                               num_people){


  normal_exposure <- exposure_sims$normal
  obese_exposure <- exposure_sims$obese

  normal_params <- exposure_params$normal
  obese_params <- exposure_params$obese


  normal_httk <- data.table(individual = rep(1:(num_people+1), length(normal_exposure)),
                           Cplasma_max = unname(unlist(lapply(normal_exposure, function(t) {sapply(t$numeric, function(j) {max(j$Cplasma)})}))),
                           AUC = unname(unlist(lapply(normal_exposure, function(t) {sapply(t$numeric, function(j) {max(j$AUC)})}))),
                           BW = c(unname(unlist(lapply(normal_params, function(t) {c(t$BW, mean(t$BW))})))),
                           Age = rep(c(10*(1+1:length(normal_exposure))), each = (num_people+1)),
                           Scenario = Scenario,
                           Weight = 'Normal')

  #print(head(normal_httk, 10))

  cplasma_max_normal <- ggplot(normal_httk[-c((num_people+1)*1:length(normal_exposure)),], aes(x = Cplasma_max, y = as.character(Age), color = Age, fill = Age)) + ggridges::geom_density_ridges(alpha = 0.5) + xlab('Max Cplasma') + ylab('Age cohort') + labs(fill = 'Age cohort', color = 'Age cohort')
  auc_normal <- ggplot(normal_httk[-c((num_people+1)*1:length(normal_params)),], aes(x = AUC, y = as.character(Age), color = Age, fill = Age)) + ggridges::geom_density_ridges(alpha = 0.5) + xlab('AUC') + ylab('Age cohort') + labs(fill = 'Age cohort', color = 'Age cohort')



  obese_httk <- data.table(individual = rep(1:(num_people+1), length(obese_exposure)),
                           Cplasma_max = unname(unlist(lapply(obese_exposure, function(t) {sapply(t$numeric, function(j) {max(j$Cplasma)})}))),
                           AUC = unname(unlist(lapply(obese_exposure, function(t) {sapply(t$numeric, function(j) {max(j$AUC)})}))),
                           BW = c(unname(unlist(lapply(obese_params, function(t) {c(t$BW, mean(t$BW))})))),
                           Age = rep(c(10*(1+1:length(obese_exposure))), each = (num_people+1)),
                           Scenario = Scenario,
                           Weight = 'Obese')

  cplasma_max_obese <- ggplot(obese_httk[-c((num_people+1)*1:length(obese_exposure)),], aes(x = Cplasma_max, y = as.character(Age), color = Age, fill = Age)) + ggridges::geom_density_ridges(alpha = 0.5) + xlab('Max Cplasma') + ylab('Age cohort') + labs(fill = 'Age cohort', color = 'Age cohort')
  auc_obese <- ggplot(obese_httk[-c((num_people+1)*1:length(obese_exposure)),], aes(x = AUC, y = as.character(Age), color = Age, fill = Age)) + ggridges::geom_density_ridges(alpha = 0.5) + xlab('AUC') + ylab('Age cohort') + labs(fill = 'Age cohort', color = 'Age cohort')

  return(list('normal' = list('normal_httk' = normal_httk,
                       'plots' = list('cplasma_max_normal' = cplasma_max_normal,
                                      'auc_normal' = auc_normal)),
       'obese' = list('obese_httk' = obese_httk,
                      plots = list('cplasma_max_obese' = cplasma_max_obese,
                                   'auc_obese' = auc_obese))))
}

run_simulations <- function(chem.cas = '',
                            num_people = num_people,
                            overwrite = FALSE){
  current_path <- paste0(getwd(), '/inst/SOT_simulations')
  if (dir.exists(paste0(current_path, '/', chem.cas)) & !overwrite){
    stop('Directory already exists! Either move contents of directory and change chemical!')
  }

  dir.create(paste0(current_path, '/', chem.cas))

populations_test <- pop_simulator(chem.cas = chem.cas, samples = num_people)
parameters_test <- simulate_parameters(chem.cas = chem.cas, mcs = populations_test)
acute_test <- acute_exposure(chem.cas = chem.cas, sim_parms = parameters_test, acute_matrix = acute_exposure_matrix)
periodic_test <- periodic_exposure(chem.cas = chem.cas, sim_parms = parameters_test, periodic_matrix = periodic_exposure_matrix)
constant_test <- constant_exposure(chem.cas = chem.cas, sim_parms = parameters_test, constant_matrix = constant_exposure_matrix)
acute_distributions <- httk_distributions(exposure_sims = acute_test, exposure_params = parameters_test, Scenario = 'Acute', num_people = num_people)
periodic_distributions <- httk_distributions(exposure_sims = periodic_test, exposure_params = parameters_test, Scenario = 'Periodic', num_people = num_people)
constant_distributions <- httk_distributions(exposure_sims = constant_test, exposure_params = parameters_test, Scenario = 'Constant', num_people = num_people)


# print(paste0(current_path, '/', chem.cas, '/', 'Cplasma_AUC_acute.pdf'))
# plot_grid(acute_distributions$normal$plots$auc_normal + scale_x_log10() + labs(title = paste('Distribution of Cplasma AUC values, normal weight, acute exposure', 'CASRN', chem.cas)),
#           acute_distributions$obese$plots$auc_obese + scale_x_log10() + labs(title = paste('Distribution of Cplasma AUC values, obese weight, acute exposure', 'CASRN', chem.cas)))
# Save plots as PDF
ggsave(filename = paste0(current_path, '/', chem.cas, '/', 'Cplasma_AUC_acute.pdf'),
       plot = plot_grid(acute_distributions$normal$plots$auc_normal + scale_x_log10() + labs(title = paste('Distribution of Cplasma AUC values, normal weight, acute exposure', 'CASRN', chem.cas)),
                        acute_distributions$obese$plots$auc_obese + scale_x_log10() + labs(title = paste('Distribution of Cplasma AUC values, obese weight, acute exposure', 'CASRN', chem.cas))),
       width = 20,
       height = 14,
       units = 'in')
# pdf(file = paste0(current_path, '/', chem.cas, '/', 'Cplasma_AUC_acute.pdf'), width = 20, height = 14)
# plot_grid(acute_distributions$normal$plots$auc_normal + scale_x_log10() + labs(title = paste('Distribution of Cplasma AUC values, normal weight, acute exposure', 'CASRN', chem.cas)),
#           acute_distributions$obese$plots$auc_obese + scale_x_log10() + labs(title = paste('Distribution of Cplasma AUC values, obese weight, acute exposure', 'CASRN', chem.cas)))
# dev.off()

ggsave(filename = paste0(current_path, '/', chem.cas, '/', 'Cplasma_AUC_periodic.pdf'),
       plot = plot_grid(periodic_distributions$normal$plots$auc_normal + scale_x_log10() + labs(title = paste('Distribution of Cplasma AUC values, normal weight, periodic exposure', 'CASRN', chem.cas)),
                        periodic_distributions$obese$plots$auc_obese + scale_x_log10() + labs(title = paste('Distribution of Cplasma AUC values, obese weight, periodic exposure', 'CASRN', chem.cas))),
       width = 20,
       height = 14,
       units = 'in'
)
# pdf(file = paste0(current_path, '/', chem.cas, '/', 'Cplasma_AUC_periodic.pdf'), width = 20, height = 14)
# plot_grid(periodic_distributions$normal$plots$auc_normal + scale_x_log10() + labs(title = paste('Distribution of Cplasma AUC values, normal weight, periodic exposure', 'CASRN', chem.cas)),
#           periodic_distributions$obese$plots$auc_obese + scale_x_log10() + labs(title = paste('Distribution of Cplasma AUC values, obese weight, periodic exposure', 'CASRN', chem.cas)))
# dev.off()

ggsave(filename = paste0(current_path, '/', chem.cas, '/', 'Cplasma_AUC_constant.pdf'),
       plot = plot_grid(constant_distributions$normal$plots$auc_normal + scale_x_log10() + labs(title = paste('Distribution of Cplasma AUC values, normal weight, constant exposure', 'CASRN', chem.cas)),
                        constant_distributions$obese$plots$auc_obese + scale_x_log10() + labs(title = paste('Distribution of Cplasma AUC values, obese weight, constant exposure', 'CASRN', chem.cas))),
       width = 20,
       height = 14,
       units = 'in')

# pdf(file = paste0(current_path, '/', chem.cas, '/', 'Cplasma_AUC_constant.pdf'), width = 20, height = 14)
# plot_grid(constant_distributions$normal$plots$auc_normal + scale_x_log10() + labs(title = paste('Distribution of Cplasma AUC values, normal weight, constant exposure', 'CASRN', chem.cas)),
#           constant_distributions$obese$plots$auc_obese + scale_x_log10() + labs(title = paste('Distribution of Cplasma AUC values, obese weight, constant exposure', 'CASRN', chem.cas)))
# dev.off()

ggsave(filename = paste0(current_path, '/', chem.cas, '/', 'Max_Cplasma_acute.pdf'),
       plot = plot_grid(acute_distributions$normal$plots$cplasma_max_normal + scale_x_log10() + labs(title = paste('Distribution of maximum Cplasma values, normal weight, acute exposure', 'CASRN', chem.cas)),
                        acute_distributions$obese$plots$cplasma_max_obese + scale_x_log10() + labs(title = paste('Distribution of maximum Cplasma values, obese weight, acute exposure', 'CASRN', chem.cas))),
       width = 20,
       height = 14,
       units = 'in')

# pdf(file = paste0(current_path, '/', chem.cas, '/', 'Max_Cplasma_acute.pdf'), width = 20, height = 14)
# plot_grid(acute_distributions$normal$plots$cplasma_max_normal + scale_x_log10() + labs(title = paste('Distribution of maximum Cplasma values, normal weight, acute exposure', 'CASRN', chem.cas)),
#           acute_distributions$obese$plots$cplasma_max_obese + scale_x_log10() + labs(title = paste('Distribution of maximum Cplasma values, obese weight, acute exposure', 'CASRN', chem.cas)))
# dev.off()

ggsave(filename = paste0(current_path, '/', chem.cas, '/', 'Max_Cplasma_periodic.pdf'),
       plot = plot_grid(periodic_distributions$normal$plots$cplasma_max_normal + scale_x_log10() + labs(title = paste('Distribution of maximum Cplasma values, normal weight, periodic exposure', 'CASRN', chem.cas)),
                        periodic_distributions$obese$plots$cplasma_max_obese + scale_x_log10() + labs(title = paste('Distribution of maximum Cplasma values, obese weight, periodic exposure', 'CASRN', chem.cas))),
       width = 20,
       height = 14,
       units = 'in')

# pdf(file = paste0(current_path, '/', chem.cas, '/', 'Max_Cplasma_periodic.pdf'), width = 20, height = 14)
# plot_grid(periodic_distributions$normal$plots$cplasma_max_normal + scale_x_log10() + labs(title = paste('Distribution of maximum Cplasma values, normal weight, periodic exposure', 'CASRN', chem.cas)),
#           periodic_distributions$obese$plots$cplasma_max_obese + scale_x_log10() + labs(title = paste('Distribution of maximum Cplasma values, obese weight, periodic exposure', 'CASRN', chem.cas)))
# dev.off()

ggsave(filename = paste0(current_path, '/', chem.cas, '/', 'Max_Cplasma_constant.pdf'),
       plot = plot_grid(constant_distributions$normal$plots$cplasma_max_normal + scale_x_log10() + labs(title = paste('Distribution of maximum Cplasma values, normal weight, constant exposure', 'CASRN', chem.cas)),
                        constant_distributions$obese$plots$cplasma_max_obese + scale_x_log10() + labs(title = paste('Distribution of maximum Cplasma values, obese weight, constant exposure', 'CASRN', chem.cas))),
       width = 20,
       height = 14,
       units = 'in')

# pdf(file = paste0(current_path, '/', chem.cas, '/', 'Max_Cplasma_constant.pdf'), width = 20, height = 14)
# plot_grid(constant_distributions$normal$plots$cplasma_max_normal + scale_x_log10() + labs(title = paste('Distribution of maximum Cplasma values, normal weight, constant exposure', 'CASRN', chem.cas)),
#           constant_distributions$obese$plots$cplasma_max_obese + scale_x_log10() + labs(title = paste('Distribution of maximum Cplasma values, obese weight, constant exposure', 'CASRN', chem.cas)))
# dev.off()

acute_distributions$normal$plots <- NULL
acute_distributions$obese$plots <- NULL
periodic_distributions$normal$plots <- NULL
periodic_distributions$obese$plots <- NULL
constant_distributions$normal$plots <- NULL
constant_distributions$obese$plots <- NULL

saveRDS(populations_test, file = paste0(current_path, '/', chem.cas, '/', 'population.RDS'))
saveRDS(parameters_test, file = paste0(current_path, '/', chem.cas, '/', 'parameters.RDS'))
saveRDS(acute_test, file = paste0(current_path, '/', chem.cas, '/', 'acute_exposure.RDS'))
saveRDS(periodic_test, file = paste0(current_path, '/', chem.cas, '/', 'periodic_exposure.RDS'))
saveRDS(constant_test, file = paste0(current_path, '/', chem.cas, '/', 'constant_exposure.RDS'))
saveRDS(acute_distributions, file = paste0(current_path, '/', chem.cas, '/', 'acute_distribution.RDS'))
saveRDS(periodic_distributions, file = paste0(current_path, '/', chem.cas, '/', 'periodic_distribution.RDS'))
saveRDS(constant_distributions, file = paste0(current_path, '/', chem.cas, '/', 'constant_distribution.RDS'))
#print(paste0(getwd(), '/inst/SOT_simulations'))
}

httk_normal_steady_state <- function(n_people,
                                    n_cohorts,
                                    chemical,
                                    parameters){
       httk_steady_state_simulation(n_people = n_people,
                                   n_cohorts = n_cohorts,
                                   chemical = chemical,
                                   parameters = parameters$normal,
                                   weight = 'Normal')    
                                    }



httk_steady_state_simulation <- function(n_people,
                                    n_cohorts,
                                    chemical,
                                    parameters,
                                    weight
                                    ){
       print(n_people)
       print(n_cohorts)
       num_cohorts <- length(n_cohorts)
       people_cohorts = n_people*num_cohorts
       df <- data.frame(avg = numeric(people_cohorts),
                                        frac = numeric(people_cohorts), 
                                        max = numeric(people_cohorts), 
                                        the.day = numeric(people_cohorts), 
                                        Age = numeric(people_cohorts), 
                                        individual = numeric(people_cohorts))
       cohort_min <- unlist(lapply(n_cohorts, min))
       print(cohort_min)
for (i in 1:n_people){
  for (j in 1:num_cohorts) {
       df[i + (j-1)*n_people,] <- cbind(as.data.frame(httk::calc_css(chem.cas = chemical, 
                                                                     parameters = parameters[[j]][i,])), Age = cohort_min[j])
}
  }
df$casn <- chemical
df$weight <- weight
df$individual <- rep(1:n_people, num_cohorts)
return(df)
}

# To merge tables, follow the procedure below for both the 'normal' and 'obese' weight classifications.
# normal_101_14_4 <- merge.data.table(merge.data.table(acute_distribution_101_14_4$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_101_14_4$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_101_14_4$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '101-14-4']
# normal_101_77_9 <- merge.data.table(merge.data.table(acute_distribution_101_77_9$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_101_77_9$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_101_77_9$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '101-77-9']
# normal_106_50_3 <- merge.data.table(merge.data.table(acute_distribution_106_50_3$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_106_50_3$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_106_50_3$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '106-50-3']
# normal_107_21_1 <- merge.data.table(merge.data.table(acute_distribution_107_21_1$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_107_21_1$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_107_21_1$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '107-21-1']
# normal_111_44_4 <- merge.data.table(merge.data.table(acute_distribution_111_44_4$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_111_44_4$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_111_44_4$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '111-44-4']
# normal_114_26_1 <- merge.data.table(merge.data.table(acute_distribution_114_26_1$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_114_26_1$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_114_26_1$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '114-26-1']
# normal_117_81_7 <- merge.data.table(merge.data.table(acute_distribution_117_81_7$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_117_81_7$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_117_81_7$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '117-81-7']
# normal_119_90_4 <- merge.data.table(merge.data.table(acute_distribution_119_90_4$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_119_90_4$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_119_90_4$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '119-90-4']
# normal_120_80_9 <- merge.data.table(merge.data.table(acute_distribution_120_80_9$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_120_80_9$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_120_80_9$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '120-80-9']
# normal_121_14_2 <- merge.data.table(merge.data.table(acute_distribution_121_14_2$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_121_14_2$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_121_14_2$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '121-14-2']
# normal_121_69_7 <- merge.data.table(merge.data.table(acute_distribution_121_69_7$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_121_69_7$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_121_69_7$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '121-69-7']
# normal_123_31_9 <- merge.data.table(merge.data.table(acute_distribution_123_31_9$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_123_31_9$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_123_31_9$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '123-31-9']
# normal_131_11_3 <- merge.data.table(merge.data.table(acute_distribution_131_11_3$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_131_11_3$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_131_11_3$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '131-11-3']
# normal_133_90_4 <- merge.data.table(merge.data.table(acute_distribution_133_90_4$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_133_90_4$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_133_90_4$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '133-90-4']
# normal_51_28_5 <- merge.data.table(merge.data.table(acute_distribution_51_28_5$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_51_28_5$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_51_28_5$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '51-28-5']
# normal_51_79_6 <- merge.data.table(merge.data.table(acute_distribution_51_79_6$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_51_79_6$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_51_79_6$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '51-79-6']
# normal_510_15_6 <- merge.data.table(merge.data.table(acute_distribution_510_15_6$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_510_15_6$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_510_15_6$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '510-15-6']
# normal_53_96_3 <- merge.data.table(merge.data.table(acute_distribution_53_96_3$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_53_96_3$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_53_96_3$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '53_96_3']
# normal_532_27_4 <- merge.data.table(merge.data.table(acute_distribution_532_27_4$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_532_27_4$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_532_27_4$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '532-27-4']
# normal_53_96_3 <- merge.data.table(merge.data.table(acute_distribution_53_96_3$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_53_96_3$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_53_96_3$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '53-96-3']
# normal_534_52_1 <- merge.data.table(merge.data.table(acute_distribution_534_52_1$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_534_52_1$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_534_52_1$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '534-52-1']
# normal_56_38_2 <- merge.data.table(merge.data.table(acute_distribution_56_38_2$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_56_38_2$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_56_38_2$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '56-38-2']
# normal_584_84_9 <- merge.data.table(merge.data.table(acute_distribution_584_84_9$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_584_84_9$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_584_84_9$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '584-84-9']
# normal_60_11_7 <- merge.data.table(merge.data.table(acute_distribution_60_11_7$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_60_11_7$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_60_11_7$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '60-11-7']
# normal_63_25_2 <- merge.data.table(merge.data.table(acute_distribution_63_25_2$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_63_25_2$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_63_25_2$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '63-25-2']
# normal_64_67_5 <- merge.data.table(merge.data.table(acute_distribution_64_67_5$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_64_67_5$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_64_67_5$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '64-67-5']
# normal_72_43_5 <- merge.data.table(merge.data.table(acute_distribution_72_43_5$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_72_43_5$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_72_43_5$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '72-43-5']
# normal_77_78_1 <- merge.data.table(merge.data.table(acute_distribution_77_78_1$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_77_78_1$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_77_78_1$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '77-78-1']
# normal_79_44_7 <- merge.data.table(merge.data.table(acute_distribution_79_44_7$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_79_44_7$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_79_44_7$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '79-44-7']
# normal_822_06_0 <- merge.data.table(merge.data.table(acute_distribution_822_06_0$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_822_06_0$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_822_06_0$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '822-06-0']
# normal_84_74_2 <- merge.data.table(merge.data.table(acute_distribution_84_74_2$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_84_74_2$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_84_74_2$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '84-74-2']
# normal_87_86_5 <- merge.data.table(merge.data.table(acute_distribution_87_86_5$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_87_86_5$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_87_86_5$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '87-86-5']
# normal_88_06_2 <- merge.data.table(merge.data.table(acute_distribution_88_06_2$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_88_06_2$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_88_06_2$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '88-06-2']
# normal_91_22_5 <- merge.data.table(merge.data.table(acute_distribution_91_22_5$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_91_22_5$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_91_22_5$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '91-22-5']
# normal_92_87_5 <- merge.data.table(merge.data.table(acute_distribution_92_87_5$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_92_87_5$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_92_87_5$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '92-87-5']
# normal_94_75_7 <- merge.data.table(merge.data.table(acute_distribution_94_75_7$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_94_75_7$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_94_75_7$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '94-75-7']
# normal_95_80_7 <- merge.data.table(merge.data.table(acute_distribution_95_80_7$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_95_80_7$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_95_80_7$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '95-80-7']
# normal_95_95_4 <- merge.data.table(merge.data.table(acute_distribution_95_95_4$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_95_95_4$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_95_95_4$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '95-95-4']
# normal_96_09_3 <- merge.data.table(merge.data.table(acute_distribution_96_09_3$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_96_09_3$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_96_09_3$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '96-09-3']
# normal_98_86_2 <- merge.data.table(merge.data.table(acute_distribution_98_86_2$normal$normal_httk[, index := paste(individual, Age)], constant_distribution_98_86_2$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')), periodic_distribution_98_86_2$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := '98-86-2']
# normal_table <- rbind(normal_101_14_4, normal_101_77_9, normal_106_50_3, normal_107_21_1, normal_111_44_4, normal_114_26_1, normal_117_81_7, normal_119_90_4, normal_120_80_9, normal_121_14_2, normal_121_69_7, normal_123_31_9, normal_131_11_3, normal_133_90_4, normal_51_28_5, normal_51_79_6, normal_510_15_6, normal_532_27_4, normal_534_52_1, normal_53_96_3, normal_56_38_2, normal_584_84_9, normal_60_11_7, normal_63_25_2, normal_64_67_5, normal_72_43_5, normal_77_78_1, normal_79_44_7, normal_822_06_0, normal_84_74_2, normal_87_86_5, normal_88_06_2, normal_91_22_5, normal_92_87_5, normal_94_75_7, normal_95_80_7, normal_95_95_4, normal_96_09_3, normal_98_86_2)
# normal_table[individual != 501, .(Cplasma_max_ratio_ap = mean(Cplasma_max.a/Cplasma_max), Cplasma_max_ratio_ac = mean(Cplasma_max.a/Cplasma_max.c), AUC_ratio_ap = mean(AUC.a/AUC), AUC_ratio_ac = mean(AUC.a/AUC.c), Cplasma_max_ap = sum(Cplasma_max.a > Cplasma_max), Cplasma_max_pc = sum(Cplasma_max > Cplasma_max.c), AUC_pc = sum(AUC > AUC.c)), by = .(casn, Age)]

# Prepare the time-to-steadystate values for 'normal' and 'obese'
# df_normal_532_27_4 <- data.frame(avg = numeric(3000), frac = numeric(3000), max = numeric(3000), the.day = numeric(3000), Age = numeric(3000), individual = numeric(3000))
# for (i in 1:500){ df_normal_532_27_4[i,] <- cbind(as.data.frame(httk::calc_css(chem.cas = '532-27-4', parameters = parameters_532_27_4$normal$norm_20_param[i,])), Age = 20)
#  df_normal_532_27_4[i + 500,] <- cbind(as.data.frame(httk::calc_css(chem.cas = '532-27-4', parameters = parameters_532_27_4$normal$norm_30_param[i,])), Age = 30)
#  df_normal_532_27_4[i + 1000,] <- cbind(as.data.frame(httk::calc_css(chem.cas = '532-27-4', parameters = parameters_532_27_4$normal$norm_40_param[i,])), Age = 40)
#  df_normal_532_27_4[i + 1500,] <- cbind(as.data.frame(httk::calc_css(chem.cas = '532-27-4', parameters = parameters_532_27_4$normal$norm_50_param[i,])), Age = 50)
#  df_normal_532_27_4[i + 2000,] <- cbind(as.data.frame(httk::calc_css(chem.cas = '532-27-4', parameters = parameters_532_27_4$normal$norm_60_param[i,])), Age = 60)
#  df_normal_532_27_4[i + 2500,] <- cbind(as.data.frame(httk::calc_css(chem.cas = '532-27-4', parameters = parameters_532_27_4$normal$norm_70_param[i,])), Age = 70)
#  if (i%%10 == 0){print(i)}
# }
# df_normal_532_27_4$casn <- '532-27-4'
# df_normal_532_27_4$weight <- 'Normal'
# df_normal_532_27_4$individual <- rep(1:500, 6)
# normal_weight_css_table <- rbind(df_normal_51_28_5, df_normal_51_79_6, df_normal_53_96_3, df_normal_56_38_2,
#                                  df_normal_60_11_7, df_normal_63_25_2, df_normal_64_67_5, df_normal_72_43_5,
#                                  df_normal_77_78_1, df_normal_79_44_7, df_normal_84_74_2, df_normal_87_86_5,
#                                  df_normal_88_06_2, df_normal_91_22_5, df_normal_92_87_5, df_normal_94_75_7,
#                                  df_normal_95_80_7, df_normal_95_95_4, df_normal_96_09_3, df_normal_98_86_2,
#                                  df_normal_101_14_4, df_normal_101_77_9, df_normal_106_50_3, df_normal_107_21_1,
#                                  df_normal_111_44_4, df_normal_114_26_1, df_normal_117_81_7, df_normal_119_90_4,
#                                  df_normal_120_80_9, df_normal_121_14_2, df_normal_121_69_7, df_normal_123_31_9,
#                                  df_normal_131_11_3, df_normal_133_90_4, df_normal_510_15_6, df_normal_532_27_4,
#                                  df_normal_534_52_1, df_normal_584_84_9, df_normal_822_06_0)

# Analysis
# setnames(normal_table, old = c('individual', 'Cplasma_max', 'AUC', 'BW', 'Age', 'Scenario', 'Weight'),
#          new = paste0(c('individual', 'Cplasma_max', 'AUC', 'BW', 'Age', 'Scenario', 'Weight'), '.p'))
# normal_table[individual.a != 501, .(Cplasma_max_ratio_ap = mean(Cplasma_max.a/Cplasma_max.p),
#                                     Cplasma_max_ratio_ac = mean(Cplasma_max.a/Cplasma_max.c),
#                                     AUC_ratio_ap = mean(AUC.a/AUC.p), AUC_ratio_ac = mean(AUC.a/AUC.c),
#                                     Cplasma_max_ap = sum(Cplasma_max.a > Cplasma_max.p),
#                                     Cplasma_max_pc = sum(Cplasma_max.p > Cplasma_max.c),
#                                     AUC_pc = sum(AUC.p > AUC.c)), by = .(casn)][, .(casn, Cplasma_max_ratio_ap, AUC_ratio_ap)][order(Cplasma_max_ratio_ap),]

# merge.data.table(normal_table, normal_weight_css_table,
#                  by.x = c('individual.a', 'Age.a', 'casn'),
#                  by.y = c('individual', 'Age', 'casn'))[, .(Ratio_cplasma_max_ap = mean(Cplasma_max.a/Cplasma_max.p),
#                                                             Ratio_AUC_ap = mean(AUC.a/AUC.p), css_day = mean(the.day)),
#                                                         by = .(casn)]

