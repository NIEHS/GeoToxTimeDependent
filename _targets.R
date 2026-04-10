# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(crew)
library(mirai)
library(nanonext)
library(crew.cluster)
library(GeoToxTimeDependent)
library(httk)
library(data.table)
library(ggplot2)
# library(tarchetypes) # Load other packages as needed.

load_httk_data()

controller_local <- crew::crew_controller_local(
  name = "controller_local",
  options_local = crew::crew_options_local(
    log_directory = "slurm/",
    log_join = FALSE
  ),
  workers = 1
)

controller_01 <- crew.cluster::crew_controller_slurm(
  name = "controller_01",
  workers = 50,
  # This controlls the number of workers that can be used simultaneously
  options_cluster = crew.cluster::crew_options_slurm(
    verbose = TRUE,
    memory_gigabytes_required = 16,
    log_output = "slurm/slurm-crew-%j.out",
    log_error = "slurm/slurm-crew-%j.err",
    script_lines = c(
      "#SBATCH --job-name=dispatch"
    ),
    #partition = "geo"
    #partition =  "normal"
    partition = "highmem"
  ),
  tasks_max = Inf
)

controller_02 <- crew.cluster::crew_controller_slurm(
  name = "controller_02",
  workers = 50,
  # This controlls the number of workers that can be used simultaneously
  options_cluster = crew.cluster::crew_options_slurm(
    verbose = TRUE,
    memory_gigabytes_per_cpu = 2,
    cpus_per_task = 32,
    log_output = "slurm/slurm-crew-%j.out",
    log_error = "slurm/slurm-crew-%j.err",
    script_lines = c(
      "#SBATCH --job-name=threading",
      "#SBATCH --hint=multithread",
      "export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK"
    ),
    #partition = "geo"
    #partition =  "normal"
    partition = "highmem"
  ),
  tasks_max = Inf
)


#Set target options:
tar_option_set(
  packages = c('GeoToxTimeDependent', 'httk', 'data.table', 'ggplot2'),
  imports = c('GeoToxTimeDependent'),
  controller = crew_controller_group(controller_local, controller_01, controller_02),
  resources = tar_resources(
    crew = tar_resources_crew(controller = "controller_local")
  ),
  format = 'qs',
  storage = "worker", 
  retrieval = "worker"
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source(files = 'inst/SOT simulations.R')
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  tar_target(
    name = chemicals,
    command = c('51-28-5', '51-79-6', '53-96-3', '56-38-2', '60-11-7', 
                '63-25-2', '64-67-5', '72-43-5', '77-78-1', '79-44-7', 
                '84-74-2', '87-86-5', '88-06-2', '91-22-5', '92-87-5', 
                '94-75-7', '95-80-7', '95-95-4', '96-09-3', '98-86-2', 
                '101-14-4', '101-77-9', '106-50-3', '107-21-1', '111-44-4', 
                '114-26-1', '117-81-7', '119-90-4', '120-80-9', '121-14-2', 
                '121-69-7', '123-31-9', '131-11-3', '133-90-4', '510-15-6', 
                '532-27-4', '534-52-1', '584-84-9', '822-06-0'),
    iteration = 'list'
    )
    # format = "qs" # Efficient storage for general data objects.
  ,
  tar_target(
    name = number_people,
    command = c(500),
    iteration = 'list'
    )
    ,
    tar_target(
    name = age_limits,
    list(c(20,29), c(30, 39), c(40,49), c(50, 59), c(60, 69), c(70, 79)),
    iteration = 'list'
    ),
    tar_target(
    name = population,
    command = {GeoToxTimeDependent::load_httk_data()
      pop_simulator(chem.cas = chemicals, samples = number_people)
      },
    pattern = cross(chemicals, number_people),
    iteration = 'list',
     resources = targets::tar_resources(
       crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
     )
     ),
    # tar_target(
    # name = population_pairs,
    # command = c(chemicals, number_people),
    # pattern = cross(chemicals, number_people),
    #  iteration = 'list'#,
    # # resources = targets::tar_resources(
    # #   crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
    # # )
    #  ),
     tar_target(
    name = simulate_params,
    command = {load_httk_data()
      simulate_parameters(chem.cas = chemicals, mcs = population)
    },
    pattern = map(chemicals, population),
    iteration = 'list',
    resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
    )
     ),
     tar_target(
    name = simulate_params_test,
    command = {#load_httk_data()
      simulate_parameters(chem.cas = chemicals, mcs = population)
    },
    pattern = map(chemicals, population),
    iteration = 'list',
    resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
    )
     ),
     tar_target(
      acute_exp, 
      command = acute_exposure_matrix,
      iteration = 'list'
     ),
      tar_target(
      periodic_exp, 
      command = periodic_exposure_matrix,
      iteration = 'list'
     ),
      tar_target(
      constant_exp, 
      command = constant_exposure_matrix,
      iteration = 'list'
     ),
     tar_target(
      acute_scenario,
      command = acute_exposure(chem.cas = chemicals,
                               sim_parms = simulate_params,
                               acute_matrix = acute_exp),
      pattern = map(chemicals, simulate_params),
      iteration = 'list',
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
     )
     ),
     tar_target(
      periodic_scenario,
      command = periodic_exposure(chem.cas = chemicals,
                                  sim_parms = simulate_params,
                                  periodic_matrix = periodic_exp),
      pattern = map(chemicals, simulate_params),
      iteration = 'list',
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
     )
     ),
     tar_target(
      constant_scenario,
      command = constant_exposure(chem.cas = chemicals,
                                  sim_parms = simulate_params,
                                  constant_matrix = constant_exp),
      pattern = map(chemicals, simulate_params),
      iteration = 'list',
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
     )
     ),
     tar_target(
      acute_distribution,
      command = httk_distributions(exposure_sims = acute_scenario, 
                                   exposure_params = simulate_params, 
                                   Scenario = 'Acute', 
                                   num_people = number_people),
      pattern = cross(map(acute_scenario, simulate_params), number_people),
      iteration = 'list',
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
     )
     ),
     tar_target(
      periodic_distribution,
      command = httk_distributions(exposure_sims = periodic_scenario, 
                                   exposure_params = simulate_params, 
                                   Scenario = 'Periodic', 
                                   num_people = number_people),
      pattern = cross(map(periodic_scenario, simulate_params), number_people),
      iteration = 'list',
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
     )
     ),
     tar_target(
      constant_distribution,
      command = httk_distributions(exposure_sims = constant_scenario, 
                                   exposure_params = simulate_params, 
                                   Scenario = 'Constant', 
                                   num_people = number_people),
      pattern = cross(map(constant_scenario, simulate_params), number_people),
      iteration = 'list',
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
     )
     ),
     tar_target(
      analytic_css_acute,
      command = {
        # Take total exposure over thirty days
        total_intake <- sum(acute_exp[, 2])
        css_analytic <- analytic_css(population_parameters = simulate_params, 
                                                chemical = chemicals,
                                               # model = 'pbtk',
                                                dose = total_intake)
        css_analytic
      },
      pattern = map(chemicals, simulate_params),
      iteration = 'list',
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
     )
     ),
     tar_target(
      analytic_css_periodic,
      command = {
        # Take total exposure over thirty days
        total_intake <- sum(periodic_exp[, 2])
        css_analytic <- analytic_css(population_parameters = simulate_params, 
                                                chemical = chemicals,
                                               # model = 'pbtk',
                                                dose = total_intake)
        css_analytic
      },
      pattern = map(chemicals, simulate_params),
      iteration = 'list',
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
     )
     ),
     tar_target(
      analytic_css_constant,
      command = {
        # Take total exposure over thirty days
        total_intake <- sum(constant_exp[, 2])
        css_analytic <- analytic_css(population_parameters = simulate_params, 
                                                chemical = chemicals,
                                               # model = 'pbtk',
                                                dose = total_intake)
        css_analytic
      },
      pattern = map(chemicals, simulate_params),
      iteration = 'list',
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
     )
     ),
     tar_target(
      normal_distribution,
      command = {
        casn <- chemicals
        dt <- merge.data.table(merge.data.table(acute_distribution$normal$normal_httk[, index := paste(individual, Age)],
                                          constant_distribution$normal$normal_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')),
                                          periodic_distribution$normal$normal_httk[, index := paste(individual, Age)], by = 'index')[, casn := casn]
      setnames(dt, old = c('individual', 'Cplasma_max', 'AUC', 'BW', 'Age', 'Scenario', 'Weight'),
          new = paste0(c('individual', 'Cplasma_max', 'AUC', 'BW', 'Age', 'Scenario', 'Weight'), '.p'))
      
      dt
      },
      pattern = map(map(map(acute_distribution, constant_distribution),periodic_distribution), chemicals),
      iteration = 'list',
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
     )
     ),
     tar_target(
      obese_distribution,
      command = {
        casn <- chemicals
        dt <- merge.data.table(merge.data.table(acute_distribution$obese$obese_httk[, index := paste(individual, Age)],
                                          constant_distribution$obese$obese_httk[, index := paste(individual, Age)], by = 'index', suffixes = c('.a', '.c')),
                                          periodic_distribution$obese$obese_httk[, index := paste(individual, Age)], by = 'index')[, casn := casn]
      setnames(dt, old = c('individual', 'Cplasma_max', 'AUC', 'BW', 'Age', 'Scenario', 'Weight'),
          new = paste0(c('individual', 'Cplasma_max', 'AUC', 'BW', 'Age', 'Scenario', 'Weight'), '.p'))
      
      dt
      },
      pattern = map(map(map(acute_distribution, constant_distribution),periodic_distribution), chemicals),
      iteration = 'list',
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
     )
     ),
     tar_target(
      name = normal_dist_processed,
      command = {
        dt <- rbindlist(normal_distribution)
        aggregate_person <- number_people + 1
        dt[individual.a != aggregate_person, .(Cplasma_max_ratio_ap = mean(Cplasma_max.a/Cplasma_max.p),
                                              Cplasma_max_ratio_ac = mean(Cplasma_max.a/Cplasma_max.c),
                                              AUC_ratio_ap = mean(AUC.a/AUC.p), 
                                              AUC_ratio_ac = mean(AUC.a/AUC.c),
                                              Cplasma_max_ap = sum(Cplasma_max.a > Cplasma_max.p), 
                                              Cplasma_max_pc = sum(Cplasma_max.p > Cplasma_max.c), 
                                              AUC_pc = sum(AUC.p > AUC.c)), by = .(casn, Age.a)]
      }
     ),
     tar_target(
      name = obese_dist_processed,
      command = {
        dt <- rbindlist(obese_distribution)
        aggregate_person <- number_people + 1
        dt[individual.a != aggregate_person, .(Cplasma_max_ratio_ap = mean(Cplasma_max.a/Cplasma_max.p),
                                              Cplasma_max_ratio_ac = mean(Cplasma_max.a/Cplasma_max.c),
                                              AUC_ratio_ap = mean(AUC.a/AUC.p), 
                                              AUC_ratio_ac = mean(AUC.a/AUC.c),
                                              Cplasma_max_ap = sum(Cplasma_max.a > Cplasma_max.p), 
                                              Cplasma_max_pc = sum(Cplasma_max.p > Cplasma_max.c), 
                                              AUC_pc = sum(AUC.p > AUC.c)), by = .(casn, Age.a)]
      }
     ),
     tar_target(
      name = normal_steady_state,
      command = {
        parameters <- simulate_params
        normal_parameters <- parameters$normal
        #ages <- as.vector(age_limits)
        print(age_limits)
        print(chemicals)
        load_httk_data()
        httk_steady_state_simulation(n_people = number_people,
                                     n_cohorts = age_limits,
                                     chemical = chemicals,
                                     parameters = normal_parameters,
                                     weight = 'normal')
      },
      iteration = 'list',
      pattern = map(chemicals, simulate_params),
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_02") # Specify the SLURM controller
     )
     ),
     tar_target(
      name = normal_steady_state_test,
      command = {
        parameters <- simulate_params
        normal_parameters <- parameters$normal
        #ages <- as.vector(age_limits)
        print(age_limits)
        print(chemicals)
        #load_httk_data()
        httk_steady_state_simulation(n_people = number_people,
                                     n_cohorts = age_limits,
                                     chemical = chemicals,
                                     parameters = normal_parameters,
                                     weight = 'normal')
      },
      iteration = 'list',
      pattern = map(chemicals, simulate_params),
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_02") # Specify the SLURM controller
     )
     ),
     tar_target(
      name = obese_steady_state,
      command = {
        parameters <- simulate_params
        obese_parameters <- parameters$obese
        #ages <- as.vector(age_limits)
        print(age_limits)
        print(chemicals)
        load_httk_data()
        httk_steady_state_simulation(n_people = number_people,
                                     n_cohorts = age_limits,
                                     chemical = chemicals,
                                     parameters = obese_parameters,
                                     weight = 'obese')
      },
      iteration = 'list',
      pattern = map(chemicals, simulate_params),
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_02") # Specify the SLURM controller
     )
     ),
     tar_target(
      normal_plasma_steady_state_analysis,
      command = {
        normal_table <- rbindlist(normal_distribution)
        normal_weight_css_table <- rbindlist(normal_steady_state)
        merge.data.table(normal_table, normal_weight_css_table,
                 by.x = c('individual.a', 'Age.a', 'casn'),
                 by.y = c('individual', 'Age', 'casn'))[, .(Ratio_cplasma_max_ap = mean(Cplasma_max.a/Cplasma_max.p),
                                                            Ratio_cplasma_max_ac = mean(Cplasma_max.a/Cplasma_max.c),
                                                            Ratio_AUC_ap = mean(AUC.a/AUC.p), Ratio_AUC_ac = mean(AUC.a/AUC.c),
                                                            AUC_a = mean(AUC.a), Cplasma_max_a = mean(Cplasma_max.a),
                                                            css_day = mean(the.day)),
                                                        by = .(casn)]
      }
     ),
     tar_target(
      obese_plasma_steady_state_analysis,
      command = {
        obese_table <- rbindlist(obese_distribution)
        obese_weight_css_table <- rbindlist(obese_steady_state)
        merge.data.table(obese_table, obese_weight_css_table,
                 by.x = c('individual.a', 'Age.a', 'casn'),
                 by.y = c('individual', 'Age', 'casn'))[, .(Ratio_cplasma_max_ap = mean(Cplasma_max.a/Cplasma_max.p), 
                                                            Ratio_cplasma_max_ac = mean(Cplasma_max.a/Cplasma_max.c),
                                                            Ratio_AUC_ap = mean(AUC.a/AUC.p), Ratio_AUC_ac = mean(AUC.a/AUC.c),
                                                            AUC_a = mean(AUC.a), Cplasma_max_a = mean(Cplasma_max.a),
                                                            css_day = mean(the.day)),
                                                        by = .(casn)]
      }
     ),
    #  tar_target(
    #   name = test_parameters,
    #   command = {
    #     #parameters <- simulate_params
    #     #print(chemicals)
    #     #normal_parameters <- parameters$normal
    #     #list('cohorts' = age_limits, 'par' = normal_parameters, 'chem' = chemicals)
    #     age_limits

    #   },
    #   iteration = 'list',
    #    pattern = map(chemicals, simulate_params),
    #   resources = targets::tar_resources(
    #   crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
    #  )
    #  ),
     tar_target(
      k_off,
      command = -8:-3#,
      # iteration = 'list'
     ),
     tar_target(
      k_on,
      command = 0:6#,
      # iteration = 'list'
     ),
     tar_target(
      log_AC50,
      command = -2:2#,
      # iteration = 'list'
     ),
     tar_target(
      log_kinetics,
      command = {log_kinetics <- expand.grid('k_off' = k_off,
                                             'k_on' = k_on,
                                             'logAC50' = log_AC50)
                 log_kinetics$KD <- log_kinetics$k_off - log_kinetics$k_on
                 asplit(log_kinetics, 1)
                 #sapply(dim(log_kinetics)[[1]], function(i) {log_kinetics[i,]})
                 },
      iteration = 'list'#,
    #   resources = targets::tar_resources(
    #   crew = targets::tar_resources_crew(controller = "controller_01") # Specify the SLURM controller
    #  )
     ),
  
  #   tar_target(
  #    name = dr_chems,
  #    command = c('120-80-9', '79-44-7', '117-81-7'),
  #    iteration = 'list'
  #   ),
  #   tar_target(
  #    name = dr_population,
  #    command = pop_simulator(chem.cas = dr_chems, samples = 10),
  #    iteration = 'list',
  #    pattern = map(dr_chems),
  #    resources = targets::tar_resources(
  #    crew = targets::tar_resources_crew(controller = "controller_01")
  #   )
  #   ),
  #
  #   tar_target(
  #    name = dr_parameter_sweep,
  #    command = {},
  #    pattern = cross(log_kinetics, dr_population),
  #    resources = targets::tar_resources(
  #    crew = targets::tar_resources_crew(controller = "controller_01")
  #   )
  #   ),

    tar_target(
      acute_dr_sweep,
      command = {
        max <- 100
        AC50_ <- log_kinetics[3]
        AC50 <- 10^AC50_
        n <- 1
        k_off <- 10^log_kinetics[1]
        k_on <- 10^log_kinetics[2]
        dose_response_sweep(exposure_sims = acute_scenario,
                            max = max,
                            AC50 = AC50,
                            n = n,
                            k_off = k_off,
                            k_on = k_on,
                            chemical = chemicals)
      },
      iteration = 'list',
      pattern = cross(log_kinetics, map(chemicals, acute_scenario)),
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_01")
     )
    ),

    tar_target(
      periodic_dr_sweep,
      command = {
        max <- 100
        AC50_ <- log_kinetics[3]
        AC50 <- 10^AC50_
        n <- 1
        k_off <- 10^log_kinetics[1]
        k_on <- 10^log_kinetics[2]
        dose_response_sweep(exposure_sims = periodic_scenario,
                            max = max,
                            AC50 = AC50,
                            n = n,
                            k_off = k_off,
                            k_on = k_on,
                            chemical = chemicals)
      },
      iteration = 'list',
      pattern = cross(log_kinetics, map(chemicals, periodic_scenario)),
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_01")
     )
    ),
    tar_target(
      constant_dr_sweep,
      command = {
        max <- 100
        AC50_ <- log_kinetics[3]
        AC50 <- 10^AC50_
        n <- 1
        k_off <- 10^log_kinetics[1]
        k_on <- 10^log_kinetics[2]
        dose_response_sweep(exposure_sims = constant_scenario,
                            max = max,
                            AC50 = AC50,
                            n = n,
                            k_off = k_off,
                            k_on = k_on,
                            chemical = chemicals)
      },
      iteration = 'list',
      pattern = cross(log_kinetics, map(chemicals, constant_scenario)),
      resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "controller_01")
     )
    )
 
     )#,
#  tar_target(
#    name = model,
#    command = coefficients(lm(y ~ x, data = data))
#  )

