# Figures for posters
library(ggplot2)
library(Polychrome)
library(kableExtra)
library(ggnewscale)
library(cowplot)
library(data.table)
library(targets)



#Working on the ddn
# Generate tables
normal_cplasma_auc_table <- kableExtra::kbl(tar_read(normal_plasma_steady_state_analysis)[order(Ratio_cplasma_max_ap), ] |> setnames(old = 'casn', new = 'CASN'), caption = '`Normal` weight class blood plasma ratios and time to steady-state', digits = 2) |> kableExtra::kable_classic_2(latex_options = "striped", html_font = 'Calibri', full_width = FALSE)
##normal_cplasma_auc_table |> kableExtra::save_kable(file = './NCPostdoc2026_normal_blood_plasma_table_expanded_classic_capitalized.png')
obese_cplasma_auc_table <- kableExtra::kbl(tar_read(obese_plasma_steady_state_analysis)[order(Ratio_cplasma_max_ap), ] |> setnames(old = 'casn', new = 'CASN'), caption = '`Obese` weight class blood plasma ratios and time to steady-state', digits = 2) |> kableExtra::kable_classic_2(latex_options = "striped", html_font = 'Calibri', full_width = FALSE)
##obese_cplasma_auc_table |> kableExtra::save_kable(file = './NCPostdoc2026_obese_blood_plasma_table_expanded_classic_capitalized.png')

# Determine line of best fit
 normal_auc_cplasma_fit <- nls(Ratio_AUC_ap ~ A/Ratio_cplasma_max_ap + B, tar_read(normal_plasma_steady_state_analysis)[order(Ratio_cplasma_max_ap), ], start = list(A = 1, B = 1))
 obese_auc_cplasma_fit <- nls(Ratio_AUC_ap ~ A/Ratio_cplasma_max_ap + B, tar_read(obese_plasma_steady_state_analysis)[order(Ratio_cplasma_max_ap), ], start = list(A = 1, B = 1))

# Generate Cplasma and AUC plots
# Generate blood plasma stat plots
P40 <- Polychrome::createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
P40 <- Polychrome::sortByHue(P40)
P40 <- as.vector(t(matrix(P40, ncol=4)))
names(P40) <- NULL

cplasma_max_auc_plot_normal <- ggplot(tar_read(normal_plasma_steady_state_analysis),
       aes(x = Ratio_cplasma_max_ap, y = Ratio_AUC_ap, color = as.factor(casn))) +
  geom_point() + scale_color_manual(name = 'CASN', values = P40) +
  geom_function(fun = \(x) coef(normal_auc_cplasma_fit)[[1]]/x + coef(normal_auc_cplasma_fit)[[2]], color = 'black') + theme(axis.title = element_text(size = 14))

cplasma_max_auc_plot_obese <- ggplot(tar_read(obese_plasma_steady_state_analysis),
                                      aes(x = Ratio_cplasma_max_ap, y = Ratio_AUC_ap, color = as.factor(casn))) +
  geom_point() + scale_color_manual(name = 'CASN', values = P40) +
  geom_function(fun = \(x) coef(obese_auc_cplasma_fit)[[1]]/x + coef(normal_auc_cplasma_fit)[[2]], color = 'black') + theme(axis.title = element_text(size = 14))

## ggsave('./inst/NCPostdoc2026_cplasma_max_AUC_ap_normal_plot_large_font.png',
##       plot = cplasma_max_auc_plot_normal, device = 'png', width = 7680, height = 4110, units = 'px')
## ggsave('./inst/NCPostdoc2026_cplasma_max_AUC_ap_obese_plot_large_font.png',
##       plot = cplasma_max_auc_plot_obese, device = 'png', width = 7680, height = 4110, units = 'px')

#Generate plasma plots
# The chemicals are ordered c('584-84-9', '79-44-7', '117-81-7') and the plasma data is stored
# in the targets name acute_scenario_figures, periodic_scenario_figures, constant_scenario_figures 

acute_117_81_7 <- data.table::rbindlist(tar_read(acute_scenario_figures)[[3]]$normal$acute_norm_20$numeric, fill = TRUE)
acute_117_81_7_new <- acute_117_81_7[, .(AUC = max(AUC), Cplasma = max(Cplasma)), by = iteration]
acute_117_81_7_new[, Cplasma_mod := Cplasma - Cplasma[501]]
acute_117_81_7_new[, AUC_mod := AUC - AUC[501]]
acute_117_81_7_percentile <- acute_117_81_7_new[order(Cplasma_mod), iteration][c(25, 475)]

acute_117_81_7_plot <- ggplot() + geom_line(data = merge.data.table(acute_117_81_7, acute_117_81_7_new[, .(iteration, Cplasma_mod)], by = 'iteration'),
                                           aes(.data$time, .data$Cplasma, color = as.factor(abs(.data$Cplasma_mod)), linetype = .data$person)) +
  scale_colour_grey(guide = 'none', start = 0.75, end = .85) + ggnewscale::new_scale_color() +
  geom_line(data = acute_117_81_7[is.na(iteration),], aes(.data$time, .data$Cplasma, linetype = .data$person), color = 'black', linewidth = 0.5) +
  geom_line(data = acute_117_81_7[iteration == acute_117_81_7_percentile[[1]],], aes(.data$time, .data$Cplasma, linetype = '5th percentile'), color = 'black') +
  geom_line(data = acute_117_81_7[iteration == acute_117_81_7_percentile[[2]],], aes(.data$time, .data$Cplasma, linetype = '95th percentile'), color = 'black') +
  scale_linetype_manual(values = c('general person' = 'solid', 'average_person' = 'longdash', '5th percentile' = '1F', '95th percentile' = 'dotted'),
                        labels = c('general person' = 'General Person', 'average_person' = 'Average Person'), guide = 'none') +
  labs(title = 'Acute Exposure to 117-81-7', x = 'Time (d)', y = expression("Cplasma"~(mu~"M"))) + theme(axis.title = element_text(size = 14))


periodic_117_81_7 <- data.table::rbindlist(tar_read(periodic_scenario_figures)[[3]]$normal$periodic_norm_20$numeric, fill = TRUE)
periodic_117_81_7_new <- periodic_117_81_7[, .(AUC = max(AUC), Cplasma = max(Cplasma)), by = iteration]
periodic_117_81_7_new[, Cplasma_mod := Cplasma - Cplasma[501]]
periodic_117_81_7_new[, AUC_mod := AUC - AUC[501]]
periodic_117_81_7_percentile <- periodic_117_81_7_new[order(Cplasma_mod), iteration][c(25, 475)]

periodic_117_81_7_plot <- ggplot() + geom_line(data = merge.data.table(periodic_117_81_7, periodic_117_81_7_new[, .(iteration, Cplasma_mod)], by = 'iteration'),
                                              aes(.data$time, .data$Cplasma, color = as.factor(abs(.data$Cplasma_mod)), linetype = .data$person)) +
  scale_colour_grey(guide = 'none', start = 0.75, end = .85) + ggnewscale::new_scale_color() +
  geom_line(data = periodic_117_81_7[is.na(iteration),], aes(.data$time, .data$Cplasma, linetype = .data$person), color = 'black', linewidth = 0.5) +
  geom_line(data = periodic_117_81_7[iteration == acute_117_81_7_percentile[[1]],], aes(.data$time, .data$Cplasma, linetype = '5th percentile'), color = 'black') +
  geom_line(data = periodic_117_81_7[iteration == acute_117_81_7_percentile[[2]],], aes(.data$time, .data$Cplasma, linetype = '95th percentile'), color = 'black') +
  scale_linetype_manual(values = c('general person' = 'solid', 'average_person' = 'longdash', '5th percentile' = '1F', '95th percentile' = 'dotted'),
                        labels = c('general person' = 'General Person', 'average_person' = 'Average Person'), guide = 'none') +
  labs(title = 'Periodic Exposure to 117-81-7', x = 'Time (d)', y = expression("Cplasma"~(mu~"M"))) + theme(axis.title = element_text(size = 14))

constant_117_81_7 <- data.table::rbindlist(tar_read(constant_scenario_figures)[[3]]$normal$constant_norm_20$numeric, fill = TRUE)
constant_117_81_7_new <- constant_117_81_7[, .(AUC = max(AUC), Cplasma = max(Cplasma)), by = iteration]
constant_117_81_7_new[, Cplasma_mod := Cplasma - Cplasma[501]]
constant_117_81_7_new[, AUC_mod := AUC - AUC[501]]
constant_117_81_7_percentile <- constant_117_81_7_new[order(Cplasma_mod), iteration][c(25, 475)]

constant_117_81_7_plot <- ggplot() + geom_line(data = merge.data.table(constant_117_81_7, constant_117_81_7_new[, .(iteration, Cplasma_mod)], by = 'iteration'),
                                              aes(.data$time, .data$Cplasma, color = as.factor(abs(.data$Cplasma_mod)), linetype = .data$person)) +
  scale_colour_grey(guide = 'none', start = 0.75, end = .85) + ggnewscale::new_scale_color() +
  geom_line(data = constant_117_81_7[is.na(iteration),], aes(.data$time, .data$Cplasma, linetype = .data$person), color = 'black', linewidth = 0.5) +
  geom_line(data = constant_117_81_7[iteration == acute_117_81_7_percentile[[1]],], aes(.data$time, .data$Cplasma, linetype = '5th percentile'), color = 'black') +
  geom_line(data = constant_117_81_7[iteration == acute_117_81_7_percentile[[2]],], aes(.data$time, .data$Cplasma, linetype = '95th percentile'), color = 'black') +
  scale_linetype_manual(values = c('general person' = 'solid', 'average_person' = 'longdash', '5th percentile' = '1F', '95th percentile' = 'dotted'),
                        labels = c('general person' = 'General Person', 'average_person' = 'Average Person')) +
  guides(linetype = guide_legend(title = 'Person type')) +
  labs(title = 'Constant Exposure to 117-81-7', x = 'Time (d)', y = expression("Cplasma"~(mu~"M"))) + theme(axis.title = element_text(size = 14))

row_117_81_7 <- cowplot::plot_grid(acute_117_81_7_plot, periodic_117_81_7_plot, constant_117_81_7_plot, align = 'h', nrow = 1)


acute_79_44_7 <- data.table::rbindlist(tar_read(acute_scenario_figures)[[2]]$normal$acute_norm_20$numeric, fill = TRUE)
acute_79_44_7_new <- acute_79_44_7[, .(AUC = max(AUC), Cplasma = max(Cplasma)), by = iteration]
acute_79_44_7_new[, Cplasma_mod := Cplasma - Cplasma[501]]
acute_79_44_7_new[, AUC_mod := AUC - AUC[501]]
acute_79_44_7_percentile <- acute_79_44_7_new[order(Cplasma_mod), iteration][c(25, 475)]

acute_79_44_7_plot <- ggplot() + geom_line(data = merge.data.table(acute_79_44_7, acute_79_44_7_new[, .(iteration, Cplasma_mod)], by = 'iteration'),
                                           aes(.data$time, .data$Cplasma, color = as.factor(abs(.data$Cplasma_mod)), linetype = .data$person)) +
                      scale_colour_grey(guide = 'none', start = 0.75, end = .85) + ggnewscale::new_scale_color() +
                      geom_line(data = acute_79_44_7[is.na(iteration),], aes(.data$time, .data$Cplasma, linetype = .data$person), color = 'black', linewidth = 0.5) +
                      geom_line(data = acute_79_44_7[iteration == acute_79_44_7_percentile[[1]],], aes(.data$time, .data$Cplasma, linetype = '5th percentile'), color = 'black') +
                      geom_line(data = acute_79_44_7[iteration == acute_79_44_7_percentile[[2]],], aes(.data$time, .data$Cplasma, linetype = '95th percentile'), color = 'black') +
                      scale_linetype_manual(values = c('general person' = 'solid', 'average_person' = 'longdash', '5th percentile' = '1F', '95th percentile' = 'dotted'),
                                            labels = c('general person' = 'General Person', 'average_person' = 'Average Person'), guide = 'none') +
                      labs(title = 'Acute Exposure to 79-44-7', x = 'Time (d)', y = expression("Cplasma"~(mu~"M"))) + theme(axis.title = element_text(size = 14))


periodic_79_44_7 <- data.table::rbindlist(tar_read(periodic_scenario_figures)[[2]]$normal$periodic_norm_20$numeric, fill = TRUE)
periodic_79_44_7_new <- periodic_79_44_7[, .(AUC = max(AUC), Cplasma = max(Cplasma)), by = iteration]
periodic_79_44_7_new[, Cplasma_mod := Cplasma - Cplasma[501]]
periodic_79_44_7_new[, AUC_mod := AUC - AUC[501]]
periodic_79_44_7_percentile <- periodic_79_44_7_new[order(Cplasma_mod), iteration][c(25, 475)]

periodic_79_44_7_plot <- ggplot() + geom_line(data = merge.data.table(periodic_79_44_7, periodic_79_44_7_new[, .(iteration, Cplasma_mod)], by = 'iteration'),
                                           aes(.data$time, .data$Cplasma, color = as.factor(abs(.data$Cplasma_mod)), linetype = .data$person)) +
                        scale_colour_grey(guide = 'none', start = 0.75, end = .85) + ggnewscale::new_scale_color() +
                        geom_line(data = periodic_79_44_7[is.na(iteration),], aes(.data$time, .data$Cplasma, linetype = .data$person), color = 'black', linewidth = 0.5) +
                        geom_line(data = periodic_79_44_7[iteration == acute_79_44_7_percentile[[1]],], aes(.data$time, .data$Cplasma, linetype = '5th percentile'), color = 'black') +
                        geom_line(data = periodic_79_44_7[iteration == acute_79_44_7_percentile[[2]],], aes(.data$time, .data$Cplasma, linetype = '95th percentile'), color = 'black') +
                        scale_linetype_manual(values = c('general person' = 'solid', 'average_person' = 'longdash', '5th percentile' = '1F', '95th percentile' = 'dotted'),
                                              labels = c('general person' = 'General Person', 'average_person' = 'Average Person'), guide = 'none') +
                        labs(title = 'Periodic Exposure to 79-44-7', x = 'Time (d)', y = expression("Cplasma"~(mu~"M"))) + theme(axis.title = element_text(size = 14))

constant_79_44_7 <- data.table::rbindlist(tar_read(constant_scenario_figures)[[2]]$normal$constant_norm_20$numeric, fill = TRUE)
constant_79_44_7_new <- periodic_79_44_7[, .(AUC = max(AUC), Cplasma = max(Cplasma)), by = iteration]
constant_79_44_7_new[, Cplasma_mod := Cplasma - Cplasma[501]]
constant_79_44_7_new[, AUC_mod := AUC - AUC[501]]
constant_79_44_7_percentile <- constant_79_44_7_new[order(Cplasma_mod), iteration][c(25, 475)]

constant_79_44_7_plot <- ggplot() + geom_line(data = merge.data.table(constant_79_44_7, constant_79_44_7_new[, .(iteration, Cplasma_mod)], by = 'iteration'),
                                              aes(.data$time, .data$Cplasma, color = as.factor(abs(.data$Cplasma_mod)), linetype = .data$person)) +
                        scale_colour_grey(guide = 'none', start = 0.75, end = .85) + ggnewscale::new_scale_color() +
                        geom_line(data = constant_79_44_7[is.na(iteration),], aes(.data$time, .data$Cplasma, linetype = .data$person), color = 'black', linewidth = 0.5) +
                        geom_line(data = constant_79_44_7[iteration == acute_79_44_7_percentile[[1]],], aes(.data$time, .data$Cplasma, linetype = '5th percentile'), color = 'black') +
                        geom_line(data = constant_79_44_7[iteration == acute_79_44_7_percentile[[2]],], aes(.data$time, .data$Cplasma, linetype = '95th percentile'), color = 'black') +
                        scale_linetype_manual(values = c('general person' = 'solid', 'average_person' = 'longdash', '5th percentile' = '1F', '95th percentile' = 'dotted'),
                                              labels = c('general person' = 'General Person', 'average_person' = 'Average Person')) +
                        guides(linetype = guide_legend(title = 'Person type')) +
                        labs(title = 'Constant Exposure to 79-44-7', x = 'Time (d)', y = expression("Cplasma"~(mu~"M"))) + theme(axis.title = element_text(size = 14))

row_79_44_7 <- cowplot::plot_grid(acute_79_44_7_plot, periodic_79_44_7_plot, constant_79_44_7_plot, align = 'h', nrow = 1)

acute_584_84_9 <- data.table::rbindlist(tar_read(acute_scenario_figures)[[1]]$normal$acute_norm_20$numeric, fill = TRUE)
acute_584_84_9_new <- acute_584_84_9[, .(AUC = max(AUC), Cplasma = max(Cplasma)), by = iteration]
acute_584_84_9_new[, Cplasma_mod := Cplasma - Cplasma[501]]
acute_584_84_9_new[, AUC_mod := AUC - AUC[501]]
acute_584_84_9_percentile <- acute_584_84_9_new[order(Cplasma_mod), iteration][c(25, 475)]

acute_584_84_9_plot <- ggplot() + geom_line(data = merge.data.table(acute_584_84_9, acute_584_84_9_new[, .(iteration, Cplasma_mod)], by = 'iteration'),
                                           aes(.data$time, .data$Cplasma, color = as.factor(abs(.data$Cplasma_mod)), linetype = .data$person)) +
  scale_colour_grey(guide = 'none', start = 0.75, end = .85) + ggnewscale::new_scale_color() +
  geom_line(data = acute_584_84_9[is.na(iteration),], aes(.data$time, .data$Cplasma, linetype = .data$person), color = 'black', linewidth = 0.5) +
  geom_line(data = acute_584_84_9[iteration == acute_584_84_9_percentile[[1]],], aes(.data$time, .data$Cplasma, linetype = '5th percentile'), color = 'black') +
  geom_line(data = acute_584_84_9[iteration == acute_584_84_9_percentile[[2]],], aes(.data$time, .data$Cplasma, linetype = '95th percentile'), color = 'black') +
  scale_linetype_manual(values = c('general person' = 'solid', 'average_person' = 'longdash', '5th percentile' = '1F', '95th percentile' = 'dotted'),
                        labels = c('general person' = 'General Person', 'average_person' = 'Average Person'), guide = 'none') +
  labs(title = 'Acute Exposure to 584-84-9', x = 'Time (d)', y = expression("Cplasma"~(mu~"M"))) + theme(axis.title = element_text(size = 14))


periodic_584_84_9 <- data.table::rbindlist(tar_read(periodic_scenario_figures)[[1]]$normal$periodic_norm_20$numeric, fill = TRUE)
periodic_584_84_9_new <- periodic_584_84_9[, .(AUC = max(AUC), Cplasma = max(Cplasma)), by = iteration]
periodic_584_84_9_new[, Cplasma_mod := Cplasma - Cplasma[501]]
periodic_584_84_9_new[, AUC_mod := AUC - AUC[501]]
periodic_584_84_9_percentile <- periodic_584_84_9_new[order(Cplasma_mod), iteration][c(25, 475)]

periodic_584_84_9_plot <- ggplot() + geom_line(data = merge.data.table(periodic_584_84_9, periodic_584_84_9_new[, .(iteration, Cplasma_mod)], by = 'iteration'),
                                              aes(.data$time, .data$Cplasma, color = as.factor(abs(.data$Cplasma_mod)), linetype = .data$person)) +
  scale_colour_grey(guide = 'none', start = 0.75, end = .85) + ggnewscale::new_scale_color() +
  geom_line(data = periodic_584_84_9[is.na(iteration),], aes(.data$time, .data$Cplasma, linetype = .data$person), color = 'black', linewidth = 0.5) +
  geom_line(data = periodic_584_84_9[iteration == acute_584_84_9_percentile[[1]],], aes(.data$time, .data$Cplasma, linetype = '5th percentile'), color = 'black') +
  geom_line(data = periodic_584_84_9[iteration == acute_584_84_9_percentile[[2]],], aes(.data$time, .data$Cplasma, linetype = '95th percentile'), color = 'black') +
  scale_linetype_manual(values = c('general person' = 'solid', 'average_person' = 'longdash', '5th percentile' = '1F', '95th percentile' = 'dotted'),
                        labels = c('general person' = 'General Person', 'average_person' = 'Average Person'), guide = 'none') +
  labs(title = 'Periodic Exposure to 584-84-9', x = 'Time (d)', y = expression("Cplasma"~(mu~"M"))) + theme(axis.title = element_text(size = 14))

constant_584_84_9 <- data.table::rbindlist(tar_read(constant_scenario_figures)[[1]]$normal$constant_norm_20$numeric, fill = TRUE)
constant_584_84_9_new <- periodic_584_84_9[, .(AUC = max(AUC), Cplasma = max(Cplasma)), by = iteration]
constant_584_84_9_new[, Cplasma_mod := Cplasma - Cplasma[501]]
constant_584_84_9_new[, AUC_mod := AUC - AUC[501]]
constant_584_84_9_percentile <- constant_584_84_9_new[order(Cplasma_mod), iteration][c(25, 475)]

constant_584_84_9_plot <- ggplot() + geom_line(data = merge.data.table(constant_584_84_9, constant_584_84_9_new[, .(iteration, Cplasma_mod)], by = 'iteration'),
                                              aes(.data$time, .data$Cplasma, color = as.factor(abs(.data$Cplasma_mod)), linetype = .data$person)) +
  scale_colour_grey(guide = 'none', start = 0.75, end = .85) + ggnewscale::new_scale_color() +
  geom_line(data = constant_584_84_9[is.na(iteration),], aes(.data$time, .data$Cplasma, linetype = .data$person), color = 'black', linewidth = 0.5) +
  geom_line(data = constant_584_84_9[iteration == acute_584_84_9_percentile[[1]],], aes(.data$time, .data$Cplasma, linetype = '5th percentile'), color = 'black') +
  geom_line(data = constant_584_84_9[iteration == acute_584_84_9_percentile[[2]],], aes(.data$time, .data$Cplasma, linetype = '95th percentile'), color = 'black') +
  scale_linetype_manual(values = c('general person' = 'solid', 'average_person' = 'longdash', '5th percentile' = '1F', '95th percentile' = 'dotted'),
                        labels = c('general person' = 'General Person', 'average_person' = 'Average Person')) +
  guides(linetype = guide_legend(title = 'Person type')) +
  labs(title = 'Constant Exposure to 584-84-9', x = 'Time (d)', y = expression("Cplasma"~(mu~"M"))) + theme(axis.title = element_text(size = 14))

row_584_84_9 <- cowplot::plot_grid(acute_584_84_9_plot, periodic_584_84_9_plot, constant_584_84_9_plot, align = 'h', nrow = 1)


# Save plots
# ggsave('./inst/NCPostdoc2026_exposure_normal_20s_plot_grayscale_percentiles_79_44_7.png',
#        plot = row_79_44_7, device = 'png', width = 7680, height = 4110, units = 'px')
# ggsave('./inst/NCPostdoc2026_exposure_normal_20s_plot_grayscale_percentiles_117_81_7.png',
#        plot = row_117_81_7, device = 'png', width = 7680, height = 4110, units = 'px')
# ggsave('./inst/NCPostdoc2026_exposure_normal_20s_plot_grayscale_percentiles_584_84_9.png',
#        plot = row_584_84_9, device = 'png', width = 7680, height = 4110, units = 'px')






# Dose-response parameter sweep plots
acute_79_44_7_average_person <- tar_read(acute_scenario_figures)[[2]]$normal$acute_norm_20$numeric$`Average Person`
acute_79_44_7_dr <- cowplot::plot_grid(ggplot(response_decay_exponential(plasma_data = acute_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-10, k_on = 10^-1),
                                              aes(time, response)) + geom_line() +ylim(0, 38) + 
                                              labs(title = expression(k['-']*'= 1E-10 '*s^-1*', '*k['+']*'=1E-1 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                              theme(text = element_text(size = 14)),
                                       ggplot(response_decay_exponential(plasma_data = acute_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-9, k_on = 10^0),
                                              aes(time, response)) + geom_line() +ylim(0, 38) + 
                                              labs(title = expression(k['-']*'= 1E-9 '*s^-1*', '*k['+']*'=1E0 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                              theme(text = element_text(size = 14)),
                                       ggplot(response_decay_exponential(plasma_data = acute_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-8, k_on = 10^1),
                                              aes(time, response)) + geom_line() +ylim(0, 38) + 
                                              labs(title = expression(k['-']*'= 1E-8 '*s^-1*', '*k['+']*'=1E1 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                              theme(text = element_text(size = 14)),
                                       ggplot(response_decay_exponential(plasma_data = acute_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-7, k_on = 10^2),
                                              aes(time, response)) + geom_line() +ylim(0, 38) + 
                                              labs(title = expression(k['-']*'= 1E-7 '*s^-1*', '*k['+']*'=1E2 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                              theme(text = element_text(size = 14)),
                                       ggplot(response_decay_exponential(plasma_data = acute_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-6, k_on = 10^3),
                                              aes(time, response)) + geom_line() +ylim(0, 38) + 
                                              labs(title = expression(k['-']*'= 1E-6 '*s^-1*', '*k['+']*'=1E3 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                              theme(text = element_text(size = 14)),
                                       ggplot(response_decay_exponential(plasma_data = acute_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-5, k_on = 10^4),
                                              aes(time, response)) + geom_line() +ylim(0, 38) + 
                                              labs(title = expression(k['-']*'= 1E-5 '*s^-1*', '*k['+']*'=1E4 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                              theme(text = element_text(size = 14)),
                                       nrow = 1, align = 'v')

periodic_79_44_7_average_person <- tar_read(periodic_scenario_figures)[[2]]$normal$periodic_norm_20$numeric$`Average Person`
periodic_79_44_7_dr <- cowplot::plot_grid(ggplot(response_decay_exponential(plasma_data = periodic_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-10, k_on = 10^-1),
                                                 aes(time, response)) + geom_line() +ylim(0, 38) + 
                                                 labs(title = expression(k['-']*'= 1E-10 '*s^-1*', '*k['+']*'=1E-1 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                                 theme(text = element_text(size = 14)),
                                          ggplot(response_decay_exponential(plasma_data = periodic_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-9, k_on = 10^0),
                                                 aes(time, response)) + geom_line() +ylim(0, 38) + 
                                                 labs(title = expression(k['-']*'= 1E-9 '*s^-1*', '*k['+']*'=1E0 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                                 theme(text = element_text(size = 14)),
                                          ggplot(response_decay_exponential(plasma_data = periodic_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-8, k_on = 10^1),
                                                 aes(time, response)) + geom_line() +ylim(0, 38) + 
                                                 labs(title = expression(k['-']*'= 1E-8 '*s^-1*', '*k['+']*'=1E1 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                                 theme(text = element_text(size = 14)),
                                          ggplot(response_decay_exponential(plasma_data = periodic_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-7, k_on = 10^2),
                                                 aes(time, response)) + geom_line() +ylim(0, 38) + 
                                                 labs(title = expression(k['-']*'= 1E-7 '*s^-1*', '*k['+']*'=1E2 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                                 theme(text = element_text(size = 14)),
                                          ggplot(response_decay_exponential(plasma_data = periodic_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-6, k_on = 10^3),
                                                 aes(time, response)) + geom_line() +ylim(0, 38) + 
                                                 labs(title = expression(k['-']*'= 1E-6 '*s^-1*', '*k['+']*'=1E3 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                                 theme(text = element_text(size = 14)),
                                          ggplot(response_decay_exponential(plasma_data = periodic_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-5, k_on = 10^4),
                                                 aes(time, response)) + geom_line() +ylim(0, 38) + 
                                                 labs(title = expression(k['-']*'= 1E-5 '*s^-1*', '*k['+']*'=1E4 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                                 theme(text = element_text(size = 14)),
                                          nrow = 1, align = 'v')

constant_79_44_7_average_person <- tar_read(constant_scenario_figures)[[2]]$normal$constant_norm_20$numeric$`Average Person`
constant_79_44_7_dr <- cowplot::plot_grid(ggplot(response_decay_exponential(plasma_data = constant_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-10, k_on = 10^-1),
                                                 aes(time, response)) + geom_line() +ylim(0, 38) + 
                                                 labs(title = expression(k['-']*'= 1E-10 '*s^-1*', '*k['+']*'=1E-1 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                                 theme(text = element_text(size = 14)),
                                          ggplot(response_decay_exponential(plasma_data = constant_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-9, k_on = 10^0),
                                                 aes(time, response)) + geom_line() +ylim(0, 38) + 
                                                 labs(title = expression(k['-']*'= 1E-9 '*s^-1*', '*k['+']*'=1E0 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                                 theme(text = element_text(size = 14)),
                                          ggplot(response_decay_exponential(plasma_data = constant_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-8, k_on = 10^1),
                                                 aes(time, response)) + geom_line() +ylim(0, 38) + 
                                                 labs(title = expression(k['-']*'= 1E-8 '*s^-1*', '*k['+']*'=1E1 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                                 theme(text = element_text(size = 14)),
                                          ggplot(response_decay_exponential(plasma_data = constant_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-7, k_on = 10^2),
                                                 aes(time, response)) + geom_line() +ylim(0, 38) + 
                                                 labs(title = expression(k['-']*'= 1E-7 '*s^-1*', '*k['+']*'=1E2 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                                 theme(text = element_text(size = 14)),
                                          ggplot(response_decay_exponential(plasma_data = constant_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-6, k_on = 10^3),
                                                 aes(time, response)) + geom_line() +ylim(0, 38) + 
                                                 labs(title = expression(k['-']*'= 1E-6 '*s^-1*', '*k['+']*'=1E3 '*M^-1*s^-1), x = 'Time (d)', y = 'Response') + 
                                                 theme(text = element_text(size = 14)),
                                          ggplot(response_decay_exponential(plasma_data = constant_79_44_7_average_person, max = 100, AC50 = 5, n = 1, k_off = 10^-5, k_on = 10^4),
                                                 aes(time, response)) + geom_line() +ylim(0, 38) + 
                                                 labs(title = expression(k['-']*'= 1E-5 '*s^-1*', '*k['+']*'=1E4 '*M^-1*s^-1), x = 'Time (d)', y = 'Response')+ 
                                                 theme(text = element_text(size = 14)),
                                          nrow = 1, align = 'v')

# Save plots
# ggsave('./inst/NCPostdoc2026_acute_dr_sweep_79_44_7.png',
#        plot = acute_79_44_7_dr, device = 'png', width = 7680, height = 4110, units = 'px')
# ggsave('./inst/NCPostdoc2026_periodic_dr_sweep_79_44_7.png',
#        plot = periodic_79_44_7_dr, device = 'png', width = 7680, height = 4110, units = 'px')
# ggsave('./inst/NCPostdoc2026_constant_dr_sweep_79_44_7.png',
#        plot = constant_79_44_7_dr, device = 'png', width = 7680, height = 4110, units = 'px')


#Hysteresis figures
# Use person 1, age 20, chemical 95-80-7, log_AC50 = 0, log_k_off = -6, log_k_on = 0
# These hash names may go out of date. Use tar_branches on acute_dr_sweep for given chemical (16th) and 
# parameter combination (87) [16 + 0:209*39,][87,] to capture correct acute_dr_sweep and acute_scenario branch
 
dr_data <- tar_read(acute_dr_sweep_7e77b371d11b1ea8)
plasma_data <- tar_read(acute_scenario_bb67352d80d8463c)

hysteresis_curve <- ggplot(cbind(data.table(dr_data$normal$acute_norm_20[[1]])[, .(time, response)],data.table(plasma_data$normal$acute_norm_20$numeric[[1]])[, .(Cplasma)]), aes(x = Cplasma, y = response)) + 
       geom_point() + labs(title = 'Counterclockwise hysteresis curve for 95-80-7', subtitle = expression(k['-']*' = 1E-6 '*M^-1*', '*k['+']*'=1 '*M^-1*s^-1*', '*AC['50']*' = 1'~mu~'M'), x = expression("Cplasma"~(mu~"M")), y = 'Percent of max response')
plasma_dr_curve <- ggplot(cbind(data.table(dr_data$normal$acute_norm_20[[1]])[, .(time, response)],data.table(plasma_data$normal$acute_norm_20$numeric[[1]])[, .(Cplasma)]), aes(x = time)) + 
       geom_point(aes(y = Cplasma), color = 'black') + geom_point(aes(y= response/10), color = '#fdae61') + 
       scale_y_continuous(name = expression('Cplasma'~(mu~'M')), 
                          sec.axis = sec_axis(trans = ~ . * 10, name = 'Percent of max response')) + 
       labs(title = 'Cplasma and Response curves for 95-80-7, acute exposure scenario', subtitle = expression(k['-']*' = 1E-6 '*s^-1*', '*k['+']*'=1 '*M^-1*s^-1*', '*AC['50']*' = 1'~mu~'M'), x = 'Time (d)')

#ggsave('./inst/NCPostdoc2026_hysteresis_95_80_7_horizontal.png', 
#       plot = cowplot::plot_grid(plasma_dr_curve, hysteresis_curve, nrow = 1),
#       device = 'png', width = 15360, height = 4110, units = 'px', limitsize = FALSE)
#ggsave('./inst/NCPostdoc2026_hysteresis_95_80_7_vertical.png', 
#       plot = cowplot::plot_grid(plasma_dr_curve, hysteresis_curve, nrow = 2),
#       device = 'png', width = 4110, height = 7680, units = 'px')
