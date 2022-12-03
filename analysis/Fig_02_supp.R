source(file.path('analysis', '00_config.R'))

(exp_name <- 'generic_PMC_RTSS_EIR_vaccSP_IIV')

fig2A = TRUE
fig2A_SUPP = FALSE
fig2B_SUPP = FALSE
fig2B_SUPP2 = FALSE

cases_df <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_year.csv')) %>%
  rename_with(~gsub('ipti', 'pmc', .x)) %>%
  filter(Annual_EIR < eir_max) %>%
  mutate(pmc_mode_fct = ifelse(pmc_rtss_cov == '0-0', 'None',
                               ifelse(pmc_rtss_cov == '0-0.8', 'RTS,S',
                                      ifelse(pmc_rtss_cov == '0.8-0', 'PMC-3', 'PMC-3 + RTS,S'))))

cases_df$pmc_mode_fct <- factor(cases_df$pmc_mode_fct,
                                levels = (c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')),
                                labels = (c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')))

table(cases_df$pmc_mode_fct, exclude = NULL)
table(cases_df$name, exclude = NULL)


if (fig2A_SUPP) {
  pdat <- cases_df %>%
    filter(name %in% c('clinical_cases', 'severe_cases') &
             pmc_coverage == 0 &
             rtss_coverage == 0) %>%
    ungroup()

  pplot <- pdat %>%
    mutate(name = ifelse(name == 'clinical_cases', 'Clinical cases', 'Severe cases')) %>%
    ggplot(aes(x = year, y = median_val, ymin = low_val, ymax = up_val)) +
    geom_ribbon(aes(fill = as.factor(Annual_EIR), group = Annual_EIR), alpha = 0.3) +
    geom_line(aes(col = as.factor(Annual_EIR), group = Annual_EIR), size = 0.8) +
    facet_wrap(~name, nrow = 1, scales = 'free') +
    scale_x_continuous(breaks = age_label_values, labels = age_labels) +
    scale_y_continuous(labels = comma) +
    scale_color_manual(values = getPalette) +
    scale_fill_manual(values = getPalette) +
    geom_hline(yintercept = 0, size = 0.5) +
    labs(subtitle = '',
         x = 'Age (years)',
         y = 'Annual cases\nper 1000 population',
         color = 'annual EIR',
         fill = 'annual EIR') +
    customTheme_nogrid +
    guides(colour = guide_legend(reverse = T),
           fill = guide_legend(reverse = T))

  print(pplot)

  f_save_plot(pplot, paste0('S1_Fig2A'), file.path(plot_dir), width = 8, height = 3.5, units = 'in',
              device_format = device_format)
  fwrite(pdat, file.path(plot_dir, 'csv', 'S1_Fig2A_dat.csv'))

  ## Clinical cases
  pdat %>%
    ungroup() %>%
    filter(name == 'clinical_cases') %>%
    mutate(Annual_EIR = paste0('eir ', Annual_EIR), age = paste0(' >', year),
           median = format_num(median_val),
           low = format_num(low_val),
           up = format_num(up_val),
           value = paste0(median, ' (', low, '-', up, ')')) %>%
    dplyr::select(age, Annual_EIR, value) %>%
    pivot_wider(names_from = Annual_EIR, values_from = value) %>%
    kbl() %>%
    kable_styling(bootstrap_options = tbl_opts, full_width = F, position = "left", fixed_thead = T)

  ### Severe cases
  pdat %>%
    ungroup() %>%
    filter(name == 'severe_cases') %>%
    mutate(Annual_EIR = paste0('eir ', Annual_EIR), age = paste0(' >', year),
           median = format_num(median_val),
           low = format_num(low_val),
           up = format_num(up_val),
           value = paste0(median, ' (', low, '-', up, ')')) %>%
    dplyr::select(age, Annual_EIR, value) %>%
    pivot_wider(names_from = Annual_EIR, values_from = value) %>%
    kbl() %>%
    kable_styling(bootstrap_options = tbl_opts, full_width = F, position = "left", fixed_thead = T)


  rm(pdat, pplot)

} #fig2A_SUPP

if (fig2B_SUPP) {

  pdat <- cases_df %>%
    filter(name %in% c('clinical_cases', 'severe_cases')) %>%
    mutate(yint = ifelse(name == 'clinical_cases', 5500, 70),
           yint_min = ifelse(name == 'clinical_cases', 0, 0))

  pdat_c <- pdat %>%
    filter(pmc_mode_fct == 'None') %>%
    select(-pmc_mode_fct)

  pplot <- pdat %>%
    filter(pmc_mode_fct %in% c('PMC-3', 'RTS,S')) %>%
    mutate(name = gsub('_', '', gsub('_', ' ', name))) %>%
    ggplot(aes(x = year, y = median_val, ymin = low_val, ymax = up_val)) +
    geom_hline(aes(yintercept = yint), size = 0.1, col = 'white', alpha = 0.96) +
    geom_hline(aes(yintercept = yint_min), size = 0.1, col = 'white', alpha = 0.96) +
    geom_ribbon(aes(fill = as.factor(Annual_EIR), group = Annual_EIR), alpha = 0.3) +  #,linetype="dashed"
    geom_line(aes(col = as.factor(Annual_EIR), group = Annual_EIR), size = 0.8) +  #,linetype="dashed"
    #geom_line(data =pdat_c, aes(col = as.factor(Annual_EIR), group = Annual_EIR), linetype = 'dashed', size = 0.8) +  #,linetype="dashed"
    scale_x_continuous(breaks = age_label_values, labels = age_labels) +
    scale_y_continuous(labels = comma) +
    scale_color_manual(values = getPalette) +
    scale_fill_manual(values = getPalette) +
    geom_hline(yintercept = 0, size = 0.5) +
    labs(subtitle = '',
         x = 'Age (years)',
         y = 'Annual cases \nper 1000 population',
         color = 'annual EIR',
         fill = 'annual EIR') +
    customTheme +
    guides(colour = guide_legend(reverse = T),
           fill = guide_legend(reverse = T),
           panel.spacing = unit(1.3, "lines")) +
    theme(panel.grid = element_blank()) +
    facet_wrap(name ~ pmc_mode_fct, nrow = 2, scales = 'free')

  print(pplot)

  f_save_plot(pplot, paste0('S1_Fig2A'),
              file.path(plot_dir), width = 8, height = 6, units = 'in', device_format = device_format)


} #fig2B_SUPP

if (fig2B_SUPP2) {

  pdat <- cases_df %>%
    filter(name %in% c('clinical_cases_averted', 'severe_cases_averted')) %>%
    mutate(yint = ifelse(name == 'clinical_cases', 6000, 110),
           yint_min = ifelse(name == 'clinical_cases', -500, -20))


  pplot <- pdat %>%
    filter(pmc_mode_fct %in% c('PMC-3', 'RTS,S')) %>%
    mutate(name = gsub('averted', '', gsub('_', ' ', name))) %>%
    ggplot(aes(x = year, y = median_val, ymin = low_val, ymax = up_val)) +
    geom_hline(aes(yintercept = yint), size = 0.1, col = 'white', alpha = 0.96) +
    geom_hline(aes(yintercept = yint_min), size = 0.1, col = 'white', alpha = 0.96) +
    geom_ribbon(aes(fill = as.factor(Annual_EIR), group = Annual_EIR), alpha = 0.3) +  #,linetype="dashed"
    geom_line(aes(col = as.factor(Annual_EIR), group = Annual_EIR), size = 0.8) +  #,linetype="dashed"
    scale_x_continuous(breaks = age_label_values, labels = age_labels) +
    scale_y_continuous(labels = comma) +
    scale_color_manual(values = getPalette) +
    scale_fill_manual(values = getPalette) +
    geom_hline(yintercept = 0, size = 0.5) +
    labs(subtitle = '',
         x = 'Age (years)',
         y = 'Annual cases averted\nper 1000 population',
         color = 'annual EIR',
         fill = 'annual EIR') +
    customTheme +
    guides(colour = guide_legend(reverse = T),
           fill = guide_legend(reverse = T),
           panel.spacing = unit(1.3, "lines")) +
    theme(panel.grid = element_blank()) +
    facet_wrap(name ~ pmc_mode_fct, nrow = 2, scales = 'free')

  print(pplot)

  f_save_plot(pplot, paste0('S1_fig_2B'),
              file.path(plot_dir), width = 8, height = 6, units = 'in', device_format = device_format)

  fwrite(pdat, file.path(plot_dir, 'csv', 'S1_fig_2B.csv'))


  dat_aggr %>%
    ungroup() %>%
    filter(name == 'clinical_cases' & pmc_mode_fct == 'None') %>%
    mutate(Annual_EIR = paste0('eir ', Annual_EIR), age = paste0(' >', year),
           median = format_num(median_val),
           low = format_num(low_val),
           up = format_num(up_val),
           value = paste0(median, ' (', low, '-', up, ')')) %>%
    dplyr::select(age, Annual_EIR, pmc_mode_fct, value) %>%
    pivot_wider(names_from = Annual_EIR, values_from = value) %>%
    kbl() %>%
    kable_styling(bootstrap_options = tbl_opts, full_width = F, position = "left", fixed_thead = T)


  dat_aggr %>%
    ungroup() %>%
    filter(name == 'clinical_cases_averted' & pmc_mode_fct %in% c('PMC-3', 'RTS,S')) %>%
    mutate(Annual_EIR = paste0('eir ', Annual_EIR), age = paste0(' >', year),
           median = format_num(median_val),
           low = format_num(low_val),
           up = format_num(up_val),
           value = paste0(median, ' (', low, '-', up, ')')) %>%
    dplyr::select(age, Annual_EIR, pmc_mode_fct, value) %>%
    pivot_wider(names_from = Annual_EIR, values_from = value) %>%
    kbl() %>%
    kable_styling(bootstrap_options = tbl_opts, full_width = F, position = "left", fixed_thead = T)


  dat_aggr %>%
    ungroup() %>%
    filter(name == 'severe_cases' & pmc_mode_fct == 'None') %>%
    mutate(Annual_EIR = paste0('eir ', Annual_EIR), age = paste0(' >', year),
           median = format_num(median_val),
           low = format_num(low_val),
           up = format_num(up_val),
           value = paste0(median, ' (', low, '-', up, ')')) %>%
    dplyr::select(age, Annual_EIR, pmc_mode_fct, value) %>%
    pivot_wider(names_from = Annual_EIR, values_from = value) %>%
    kbl() %>%
    kable_styling(bootstrap_options = tbl_opts, full_width = F, position = "left", fixed_thead = T)


  dat_aggr %>%
    ungroup() %>%
    filter(name == 'severe_cases_averted' & pmc_mode_fct != 'None') %>%
    mutate(Annual_EIR = paste0('eir ', Annual_EIR), age = paste0(' >', year),
           median = format_num(median_val),
           low = format_num(low_val),
           up = format_num(up_val),
           value = paste0(median, ' (', low, '-', up, ')')) %>%
    dplyr::select(age, Annual_EIR, pmc_mode_fct, value) %>%
    pivot_wider(names_from = Annual_EIR, values_from = value) %>%
    kbl() %>%
    kable_styling(bootstrap_options = tbl_opts, full_width = F, position = "left", fixed_thead = T)


} #fig2B_SUPP2
