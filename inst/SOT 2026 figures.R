# This script generates figures and tables for SOT 2026 poster.
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(data.table)



acute_79_44_7 <- data.table::rbindlist(acute_exposure_79_44_7$normal$acute_norm_20$numeric, fill = TRUE)
acute_79_44_7_new <- acute_79_44_7[, .(AUC = max(AUC), Cplasma = max(Cplasma)), by = iteration]
acute_79_44_7_new[, Cplasma_mod := Cplasma - Cplasma[501]]
acute_79_44_7_new[, AUC_mod := AUC - AUC[501]]
acute_79_44_7_percentile <- acute_79_44_7_new[order(Cplasma_mod), iteration][c(25, 475)]

acute_79_44_7_plot <- ggplot() + geom_line(data = merge.data.table(acute_79_44_7, acute_79_44_7_new[, .(iteration, Cplasma_mod)], by = 'iteration'),
                                           aes(.data$time, .data$Cplasma, linetype = .data$person), color = '#E5E5E5') +
                      #scale_colour_grey(guide = 'none', start = 0.75, end = .85) +
                      ggnewscale::new_scale_color() +
                      geom_line(data = acute_79_44_7[is.na(iteration),], aes(.data$time, .data$Cplasma, linetype = .data$person), color = 'black', linewidth = 0.5) +
                      geom_line(data = acute_79_44_7[iteration == acute_79_44_7_percentile[[1]],], aes(.data$time, .data$Cplasma, linetype = '5th_percentile'), color = 'black') +
                      geom_line(data = acute_79_44_7[iteration == acute_79_44_7_percentile[[2]],], aes(.data$time, .data$Cplasma, linetype = '95th_percentile'), color = 'black') +
                      scale_linetype_manual(values = c('general person' = 'solid', 'average_person' = 'longdash', '5th_percentile' = '1F', '95th_percentile' = 'dotted'), guide = 'none') +
                      labs(title = 'Acute Exposure to 79-44-7')


periodic_79_44_7 <- data.table::rbindlist(periodic_exposure_79_44_7$normal$periodic_norm_20$numeric, fill = TRUE)
periodic_79_44_7_new <- periodic_79_44_7[, .(AUC = max(AUC), Cplasma = max(Cplasma)), by = iteration]
periodic_79_44_7_new[, Cplasma_mod := Cplasma - Cplasma[501]]
periodic_79_44_7_new[, AUC_mod := AUC - AUC[501]]
periodic_79_44_7_percentile <- periodic_79_44_7_new[order(Cplasma_mod), iteration][c(25, 475)]

periodic_79_44_7_plot <- ggplot() + geom_line(data = merge.data.table(periodic_79_44_7, periodic_79_44_7_new[, .(iteration, Cplasma_mod)], by = 'iteration'),
                                           aes(.data$time, .data$Cplasma, linetype = .data$person), color = '#E5E5E5') +
                        #scale_colour_grey(guide = 'none', start = 0.75, end = .85) +
                        ggnewscale::new_scale_color() +
                        geom_line(data = periodic_79_44_7[is.na(iteration),], aes(.data$time, .data$Cplasma, linetype = .data$person), color = 'black', linewidth = 0.5) +
                        geom_line(data = periodic_79_44_7[iteration == acute_79_44_7_percentile[[1]],], aes(.data$time, .data$Cplasma, linetype = '5th_percentile'), color = 'black') +
                        geom_line(data = periodic_79_44_7[iteration == acute_79_44_7_percentile[[2]],], aes(.data$time, .data$Cplasma, linetype = '95th_percentile'), color = 'black') +
                        scale_linetype_manual(values = c('general person' = 'solid', 'average_person' = 'longdash', '5th_percentile' = '1F', '95th_percentile' = 'dotted'), guide = 'none') +
                        labs(title = 'Periodic Exposure to 79-44-7')#+ guides(linetype = guide_legend(title = 'Person type'))

constant_79_44_7 <- data.table::rbindlist(constant_exposure_79_44_7$normal$constant_norm_20$numeric, fill = TRUE)
constant_79_44_7_new <- periodic_79_44_7[, .(AUC = max(AUC), Cplasma = max(Cplasma)), by = iteration]
constant_79_44_7_new[, Cplasma_mod := Cplasma - Cplasma[501]]
constant_79_44_7_new[, AUC_mod := AUC - AUC[501]]
constant_79_44_7_percentile <- constant_79_44_7_new[order(Cplasma_mod), iteration][c(25, 475)]

constant_79_44_7_plot <- ggplot() + geom_line(data = merge.data.table(constant_79_44_7, constant_79_44_7_new[, .(iteration, Cplasma_mod)], by = 'iteration'),
                                              aes(.data$time, .data$Cplasma, linetype = .data$person), color = '#E5E5E5') +
                        #scale_colour_grey(guide = 'none', start = 0.75, end = .85) +
                        ggnewscale::new_scale_color() +
                        geom_line(data = constant_79_44_7[is.na(iteration),], aes(.data$time, .data$Cplasma, linetype = .data$person), color = 'black', linewidth = 0.5) +
                        geom_line(data = constant_79_44_7[iteration == acute_79_44_7_percentile[[1]],], aes(.data$time, .data$Cplasma, linetype = '5th_percentile'), color = 'black') +
                        geom_line(data = constant_79_44_7[iteration == acute_79_44_7_percentile[[2]],], aes(.data$time, .data$Cplasma, linetype = '95th_percentile'), color = 'black') +
                        scale_linetype_manual(values = c('general person' = 'solid', 'average_person' = 'longdash', '5th_percentile' = '1F', '95th_percentile' = 'dotted')) +
                        guides(linetype = guide_legend(title = 'Person type')) +
                        labs(title = 'Constant Exposure to 79-44-7')

row_79_44_7 <- cowplot::plot_grid(acute_79_44_7_plot, periodic_79_44_7_plot, constant_79_44_7_plot, align = 'h', nrow = 1)
