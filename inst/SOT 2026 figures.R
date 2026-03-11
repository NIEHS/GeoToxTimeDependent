# This script generates figures and tables for SOT 2026 poster.
library(ggplot2)
library(Polychrome)
library(kableExtra)
library(ggnewscale)
library(cowplot)
library(data.table)


# Exposure scenario plots
# Generate exposure scenarios. Same three across all chemicals, populations
time <- seq(from = 0, to = 30, by = 1/24)
dose_a <- c(rep(0, 24), 100*exp(-time[1:697]))
dose_b <- 100/19.15*abs(sin(time))
dose_c <- rep(10/3, length(time))

acute_exposure_matrix <- matrix(c(time, dose_a), ncol = 2, dimnames = list(NULL, c('time', 'dose')))
periodic_exposure_matrix <- matrix(c(time, dose_b), ncol = 2, dimnames = list(NULL, c('time', 'dose')))
constant_exposure_matrix <- matrix(c(time, dose_c), ncol = 2, dimnames = list(NULL, c('time', 'dose')))

ggplot(data = acute_exposure_matrix, aes(x = time, y = dose, color = 'Acute')) +  geom_line() +
  geom_line(data = periodic_exposure_matrix, aes(color = 'Periodic')) +
  geom_line(data = constant_exposure_matrix, aes(color = 'Constant')) +
  scale_color_manual(values = c('#b2182b', '#fdae61', '#0571b0'), name = 'Exposure scenario') +
  labs(x = 'Time (d)', y = 'Dose (mg)') + theme(axis.title = element_text(size = 14))


# Generate tables
normal_cplasma_auc_table <- kableExtra::kbl(merge.data.table(normal_weight_table, normal_weight_css_table,
                                            by.x = c('individual.a', 'Age.a', 'casn'),
                                            by.y = c('individual', 'Age', 'casn'))[individual.a != 501,
                                                                                   .(Ratio_cplasma_max_ap = mean(Cplasma_max.a/Cplasma_max.p),
                                                                                     Ratio_cplasma_max_ac = mean(Cplasma_max.a/Cplasma_max.c),
                                                                                     Ratio_AUC_ap = mean(AUC.a/AUC.p),
                                                                                     Ratio_AUC_ac = mean(AUC.a/AUC.c),
                                                                                     AUC_a = mean(AUC.a),
                                                                                     Cplasma_max_a = mean(Cplasma_max.a),
                                                                                     css_day = mean(the.day)),
                                                                                   by = .(casn)][order(Ratio_cplasma_max_ap), ] |> setnames(old = 'casn', new = 'CASN'), caption = '`Normal` weight class blood plasma ratios and time to steady-state', digits = 2) |> kableExtra::kable_classic_2(latex_options = "striped", html_font = 'Calibri', full_width = FALSE)
#normal_cplasma_auc_table |> #kableExtra::save_kable(file = 'C:/Users/krusepm/Pictures/SOT 2026 Figures/normal_blood_plasma_table_expanded_classic_capitalized.png')
obese_cplasma_auc_table <- kableExtra::kbl(merge.data.table(obese_weight_table, obese_weight_css_table,
                                           by.x = c('individual.a', 'Age.a', 'casn'),
                                           by.y = c('individual', 'Age', 'casn'))[individual.a != 501,
                                                                                  .(Ratio_cplasma_max_ap = mean(Cplasma_max.a/Cplasma_max.p),
                                                                                    Ratio_cplasma_max_ac = mean(Cplasma_max.a/Cplasma_max.c),
                                                                                    Ratio_AUC_ap = mean(AUC.a/AUC.p),
                                                                                    Ratio_AUC_ac = mean(AUC.a/AUC.c),
                                                                                    AUC_a = mean(AUC.a),
                                                                                    Cplasma_max_a = mean(Cplasma_max.a),
                                                                                    css_day = mean(the.day)),
                                                                                  by = .(casn)][order(Ratio_cplasma_max_ap), ] |> setnames(old = 'casn', new = 'CASN'), caption = '`Obese` weight class blood plasma ratios and time to steady-state', digits = 2) |> kableExtra::kable_classic_2(latex_options = "striped", html_font = 'Calibri', full_width = FALSE)
#obese_cplasma_auc_table |> #kableExtra::save_kable(file = 'C:/Users/krusepm/Pictures/SOT 2026 Figures/obese_blood_plasma_table_expanded_classic_capitalized.png')

# Generate blood plasma stat plots
P40 <- Polychrome::createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
P40 <- Polychrome::sortByHue(P40)
P40 <- as.vector(t(matrix(P40, ncol=4)))
names(P40) <- NULL

cplasma_max_auc_plot_normal <- ggplot(merge.data.table(normal_weight_table, normal_weight_css_table,
                                                       by.x = c('individual.a', 'Age.a', 'casn'),
                                                       by.y = c('individual', 'Age', 'casn'))[individual.a != 501, .(Ratio_cplasma_max_ap = mean(Cplasma_max.a/Cplasma_max.p),
                                                                                      Ratio_AUC_ap = mean(AUC.a/AUC.p),
                                                                                      css_day = mean(the.day)), by = .(casn)],
       aes(x = Ratio_cplasma_max_ap, y = Ratio_AUC_ap, color = as.factor(casn))) +
  geom_point() + scale_color_manual(name = 'CASN', values = P40) +
  geom_function(fun = \(x) 1.0224/x + 0.9604, color = 'black') + theme(axis.title = element_text(size = 14))

cplasma_max_auc_plot_obese <- ggplot(merge.data.table(obese_weight_table, obese_weight_css_table,
                                                       by.x = c('individual.a', 'Age.a', 'casn'),
                                                       by.y = c('individual', 'Age', 'casn'))[individual.a != 501, .(Ratio_cplasma_max_ap = mean(Cplasma_max.a/Cplasma_max.p),
                                                                                                                     Ratio_AUC_ap = mean(AUC.a/AUC.p),
                                                                                                                     css_day = mean(the.day)), by = .(casn)],
                                      aes(x = Ratio_cplasma_max_ap, y = Ratio_AUC_ap, color = as.factor(casn))) +
  geom_point() + scale_color_manual(name = 'CASN', values = P40) +
  geom_function(fun = \(x) 1.0235/x + 0.9609, color = 'black') + theme(axis.title = element_text(size = 14))


acute_117_81_7 <- data.table::rbindlist(acute_exposure_117_81_7$normal$acute_norm_20$numeric, fill = TRUE)
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


periodic_117_81_7 <- data.table::rbindlist(periodic_exposure_117_81_7$normal$periodic_norm_20$numeric, fill = TRUE)
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

constant_117_81_7 <- data.table::rbindlist(constant_exposure_117_81_7$normal$constant_norm_20$numeric, fill = TRUE)
constant_117_81_7_new <- periodic_117_81_7[, .(AUC = max(AUC), Cplasma = max(Cplasma)), by = iteration]
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


acute_79_44_7 <- data.table::rbindlist(acute_exposure_79_44_7$normal$acute_norm_20$numeric, fill = TRUE)
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


periodic_79_44_7 <- data.table::rbindlist(periodic_exposure_79_44_7$normal$periodic_norm_20$numeric, fill = TRUE)
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

constant_79_44_7 <- data.table::rbindlist(constant_exposure_79_44_7$normal$constant_norm_20$numeric, fill = TRUE)
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

acute_584_84_9 <- data.table::rbindlist(acute_exposure_584_84_9$normal$acute_norm_20$numeric, fill = TRUE)
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


periodic_584_84_9 <- data.table::rbindlist(periodic_exposure_584_84_9$normal$periodic_norm_20$numeric, fill = TRUE)
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

constant_584_84_9 <- data.table::rbindlist(constant_exposure_584_84_9$normal$constant_norm_20$numeric, fill = TRUE)
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
# ggsave('./inst/SOT2026_exposure_normal_20s_plot_grayscale_percentiles_79_44_7.png',
#        plot = row_79_44_7, device = 'png', width = 2560, height = 1370, units = 'px')
# ggsave('./inst/SOT2026_exposure_normal_20s_plot_grayscale_percentiles_117_81_7.png',
#        plot = row_117_81_7, device = 'png', width = 2560, height = 1370, units = 'px')
# ggsave('./inst/SOT2026_exposure_normal_20s_plot_grayscale_percentiles_584_84_9.png',
#        plot = row_584_84_9, device = 'png', width = 2560, height = 1370, units = 'px')
