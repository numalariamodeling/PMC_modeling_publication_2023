source(file.path('analysis', '00_config.R'))

(exp_name <- 'generic_PMC_RTSS_EIR_vaccSP_IIV')

Fig2A = TRUE
Fig2CD = TRUE

## Weekly or monthly age


if (Fig2CD) {
  cases_df_agegrp <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
    filter(age_group %in% c('U1', 'U2'))

  cases_df <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_year.csv')) %>%
    filter(year < 2) %>%
    mutate(age_group = paste0('>', year)) %>%
    dplyr::select(colnames(cases_df_agegrp)) %>%
    bind_rows(cases_df_agegrp) %>%
    filter(Annual_EIR < eir_max) %>%
    mutate(pmc_mode_fct = ifelse(pmc_rtss_cov == '0-0', 'None',
                                 ifelse(pmc_rtss_cov == '0-0.8', 'RTS,S',
                                        ifelse(pmc_rtss_cov == '0.8-0', 'PMC-3', 'PMC-3 + RTS,S'))))

  cases_df$pmc_mode_fct <- factor(cases_df$pmc_mode_fct,
                                  levels = (c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')),
                                  labels = (c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')))


  p1 <- cases_df %>%
    filter(age_group == 'U2' &
             pmc_mode_fct != 'None' &
             name %in% c('clinical_cases_averted', 'PE_clinical_incidence')) %>%
    ggplot() +
    geom_hline(yintercept = 0.6, alpha = 0) +
    geom_hline(yintercept = 0) +
    geom_line(aes(x = Annual_EIR, y = median_val, col = pmc_mode_fct)) +
    geom_ribbon(aes(x = Annual_EIR, ymin = low_val, ymax = up_val, fill = pmc_mode_fct), alpha = 0.3) +
    facet_wrap(~name, scales = 'free') +
    scale_x_log10(breaks = unique(dat_aggr$Annual_EIR)) +
    scale_color_manual(values = c(pmc_cols[c(1)], rtss_col, rtss_pmc_cols[2])) +
    scale_fill_manual(values = c(pmc_cols[c(1)], rtss_col, rtss_pmc_cols[2])) +
    customTheme_nogrid +
    labs(caption = '') +
    theme(legend.position = 'right')

  p2 <- cases_df %>%
    filter(age_group == 'U2' &
             pmc_mode_fct != 'None' &
             name %in% c('severe_cases_averted', 'PE_severe_incidence')) %>%
    mutate(low_val = ifelse(low_val < -0.1, -0.1, low_val)) %>%
    ggplot() +
    geom_hline(yintercept = 0.6, alpha = 0) +
    geom_hline(yintercept = 0) +
    geom_line(aes(x = Annual_EIR, y = median_val, col = pmc_mode_fct)) +
    geom_ribbon(aes(x = Annual_EIR, ymin = low_val, ymax = up_val, fill = pmc_mode_fct), alpha = 0.3) +
    facet_wrap(~name, scales = 'free') +
    scale_x_log10(breaks = unique(dat_aggr$Annual_EIR)) +
    scale_color_manual(values = c(pmc_cols[c(1)], rtss_col, rtss_pmc_cols[2])) +
    scale_fill_manual(values = c(pmc_cols[c(1)], rtss_col, rtss_pmc_cols[2])) +
    customTheme_nogrid +
    labs(caption = '') +
    theme(legend.position = 'right')


  pplot <- plot_combine(list(p1, p2), ncol = 2, labels = c('C', 'D'))
  print(pplot)

  f_save_plot(pplot, paste0('Fig2CD_U2'),
              file.path(plot_dir), width = 15, height = 3, units = 'in', device_format = device_format)


  #### Supp, U1
  p1 <- cases_df %>%
    filter(age_group == 'U1' &
             pmc_mode_fct != 'None' &
             name %in% c('clinical_cases_averted', 'PE_clinical_incidence')) %>%
    ggplot() +
    geom_hline(yintercept = 0.6, alpha = 0) +
    geom_hline(yintercept = 0) +
    geom_line(aes(x = Annual_EIR, y = median_val, col = pmc_mode_fct)) +
    geom_ribbon(aes(x = Annual_EIR, ymin = low_val, ymax = up_val, fill = pmc_mode_fct), alpha = 0.3) +
    facet_wrap(~name, scales = 'free') +
    scale_x_log10(breaks = unique(dat_aggr$Annual_EIR)) +
    scale_color_manual(values = c(pmc_cols[c(1)], rtss_col, rtss_pmc_cols[2])) +
    scale_fill_manual(values = c(pmc_cols[c(1)], rtss_col, rtss_pmc_cols[2])) +
    customTheme_nogrid +
    labs(caption = '') +
    theme(legend.position = 'right')

  p2 <- cases_df %>%
    filter(age_group == 'U1' &
             pmc_mode_fct != 'None' &
             name %in% c('severe_cases_averted', 'PE_severe_incidence')) %>%
    mutate(low_val = ifelse(low_val < -0.1, -0.1, low_val)) %>%
    ggplot() +
    geom_hline(yintercept = 0.6, alpha = 0) +
    geom_hline(yintercept = 0) +
    geom_line(aes(x = Annual_EIR, y = median_val, col = pmc_mode_fct)) +
    geom_ribbon(aes(x = Annual_EIR, ymin = low_val, ymax = up_val, fill = pmc_mode_fct), alpha = 0.3) +
    facet_wrap(~name, scales = 'free') +
    scale_x_log10(breaks = unique(dat_aggr$Annual_EIR)) +
    scale_color_manual(values = c(pmc_cols[c(1)], rtss_col, rtss_pmc_cols[2])) +
    scale_fill_manual(values = c(pmc_cols[c(1)], rtss_col, rtss_pmc_cols[2])) +
    customTheme_nogrid +
    labs(caption = '') +
    theme(legend.position = 'right')

  pplot <- plot_combine(list(p1, p2), ncol = 2, labels = c('C', 'D'))

  print(pplot)
  f_save_plot(pplot, paste0('Fig2CD_U1'),
              file.path(plot_dir), width = 15, height = 3, units = 'in', device_format = device_format)

  rm(pplot, cases_df, p1, p2)
}  #Fig2CD

##------------------------------------------
## Tables
##------------------------------------------
if (result_tables) {
  cases_df <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    filter(Annual_EIR < eir_max) %>%
    mutate(pmc_mode_fct = ifelse(pmc_rtss_cov == '0-0', 'None',
                                 ifelse(pmc_rtss_cov == '0-0.8', 'RTS,S',
                                        ifelse(pmc_rtss_cov == '0.8-0', 'PMC-3', 'PMC-3 + RTS,S'))))

  cases_df$pmc_mode_fct <- factor(cases_df$pmc_mode_fct,
                                  levels = rev(c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')),
                                  labels = rev(c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')))

  dat_aggr %>%
    ungroup() %>%
    filter(name == 'clinical_cases') %>%
    mutate(Annual_EIR = paste0('eir ', Annual_EIR),
           median = format_num(median_val),
           low = format_num(low_val),
           up = format_num(up_val),
           value = paste0(median, ' (', low, '-', up, ')')) %>%
    dplyr::select(age_group, Annual_EIR, pmc_mode_fct, value) %>%
    pivot_wider(names_from = Annual_EIR, values_from = value) %>%
    kbl() %>%
    kable_styling(bootstrap_options = tbl_opts, full_width = F, position = "left", fixed_thead = T)


  dat_aggr %>%
    ungroup() %>%
    filter(name == 'clinical_cases_averted') %>%
    mutate(Annual_EIR = paste0('eir ', Annual_EIR),
           median = format_num(median_val),
           low = format_num(low_val),
           up = format_num(up_val),
           value = paste0(median, ' (', low, '-', up, ')')) %>%
    dplyr::select(age_group, Annual_EIR, pmc_mode_fct, value) %>%
    pivot_wider(names_from = Annual_EIR, values_from = value) %>%
    kbl() %>%
    kable_styling(bootstrap_options = tbl_opts, full_width = F, position = "left", fixed_thead = T)


  dat_aggr %>%
    ungroup() %>%
    filter(name == 'severe_cases') %>%
    mutate(Annual_EIR = paste0('eir ', Annual_EIR),
           median = format_num(median_val),
           low = format_num(low_val),
           up = format_num(up_val),
           value = paste0(median, ' (', low, '-', up, ')')) %>%
    dplyr::select(age_group, Annual_EIR, pmc_mode_fct, value) %>%
    pivot_wider(names_from = Annual_EIR, values_from = value) %>%
    kbl() %>%
    kable_styling(bootstrap_options = tbl_opts, full_width = F, position = "left", fixed_thead = T)


  dat_aggr %>%
    ungroup() %>%
    filter(name == 'severe_cases_averted') %>%
    mutate(Annual_EIR = paste0('eir ', Annual_EIR),
           median = format_num(median_val),
           low = format_num(low_val),
           up = format_num(up_val),
           value = paste0(median, ' (', low, '-', up, ')')) %>%
    dplyr::select(age_group, Annual_EIR, pmc_mode_fct, value) %>%
    pivot_wider(names_from = Annual_EIR, values_from = value) %>%
    kbl() %>%
    kable_styling(bootstrap_options = tbl_opts, full_width = F, position = "left", fixed_thead = T)

}