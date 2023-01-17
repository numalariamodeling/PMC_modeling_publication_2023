##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## Fig_02.R
##---------------------

source(file.path('analysis', '_config.R'))

(exp_name <- 'generic_PMC_RTSS_EIR_vaccSP_IIV')

Fig2AB <- T
Fig2CD <- T
result_tables <- T

## Weekly or monthly age
if (Fig2AB) {
  dat_aggr <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_week.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    group_by(year) %>%
    filter(age > min(age)) %>%
    ungroup() %>%
    filter(Annual_EIR < eir_max)

  dat_aggr <- dat_aggr %>%
    group_by(pmc_rtss_cov,pmc_coverage,rtss_coverage, cm_coverage, name,  Annual_EIR) %>%
    arrange(age) %>%
    mutate(median_val_raw = median_val,
           median_val = rollmean(median_val, k = 3, fill = NA, align = 'left'),
           low_val = rollmean(low_val, k = 3, fill = NA, align = 'left'),
           up_val = rollmean(up_val, k = 3, fill = NA, align = 'left'))

  dat_aggr <- dat_aggr %>% mutate(pmc_mode_fct = ifelse(pmc_rtss_cov == '0-0', 'None',
                                                        ifelse(pmc_rtss_cov == '0-0.8', 'RTS,S',
                                                               ifelse(pmc_rtss_cov == '0.8-0', 'PMC-3', 'PMC-3 + RTS,S'))))

  dat_aggr$pmc_mode_fct <- factor(dat_aggr$pmc_mode_fct,
                                  levels = (c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')),
                                  labels = (c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')))


  pdat <- dat_aggr %>%
    filter(name %in% c('clinical_cases', 'severe_cases') &
             age < 365 * 2 &
             pmc_mode_fct != 'PMC-3 + RTS,S') %>%
    mutate(yint = ifelse(name == 'clinical_cases', 8000, 250),
           yint_min = ifelse(name == 'clinical_cases', 0, 0))

  pplot <- pdat %>%
    mutate(name = ifelse(name == 'clinical_cases', 'Clinical malaria', 'Severe malaria')) %>%
    ggplot(aes(x = age, y = median_val, ymin = low_val, ymax = up_val)) +
    geom_hline(aes(yintercept = yint), size = 0.1, col = 'white', alpha = 0.96) +
    geom_hline(aes(yintercept = yint_min), size = 0.1, col = 'white', alpha = 0.96) +
    geom_ribbon(aes(fill = as.factor(Annual_EIR), group = Annual_EIR), alpha = 0.3) +
    geom_line(aes(col = as.factor(Annual_EIR), group = Annual_EIR), size = 0.8) +
    facet_wrap(name ~ pmc_mode_fct, nrow = 2, scales = 'free') +
    scale_x_continuous(breaks = seq(0, 5 * 365, 30 * 6), labels = seq(0, 5 * 12, 6) / 12) +
    scale_y_continuous(labels = comma, expand = c(0, 0)) +
    scale_color_manual(values = getPalette) +
    scale_fill_manual(values = getPalette) +
    geom_hline(yintercept = 0, size = 0.5) +
    labs(subtitle = '',
         x = 'Age (years)',
         y = 'Cases\nper 1000 population per year',
         color = 'annual EIR',
         fill = 'annual EIR', caption='Values smoothed using 3 weeks rolling average') +
    customTheme_nogrid +
    guides(colour = guide_legend(reverse = T),
           fill = guide_legend(reverse = T))

  print(pplot)
  summary(pdat$age)

  f_save_plot(pplot, paste0('Fig2A_U2'),
              file.path(plot_dir), width = 10, height = 5, units = 'in', device_format = device_format)
  #fwrite(pdat, file.path(plot_dir, 'csv', 'Fig2A.csv'))
} # Fig2AB

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
    scale_x_log10(breaks = unique(cases_df$Annual_EIR)) +
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
    scale_x_log10(breaks = unique(cases_df$Annual_EIR)) +
    scale_color_manual(values = c(pmc_cols[c(1)], rtss_col, rtss_pmc_cols[2])) +
    scale_fill_manual(values = c(pmc_cols[c(1)], rtss_col, rtss_pmc_cols[2])) +
    customTheme_nogrid +
    labs(caption = '') +
    theme(legend.position = 'right')


  pplot <- plot_combine(list(p1, p2), ncol = 2, labels = c('C', 'D'))
  print(pplot)

  f_save_plot(pplot, paste0('Fig2CD_U2'),
              file.path(plot_dir), width = 15, height = 3, units = 'in', device_format = device_format)
  fwrite(cases_df, file.path(plot_dir, 'csv', 'Fig2CD_dat.csv'))

  rm(pplot, cases_df, p1, p2)
}  #Fig2CD

if (result_tables) {
  (exp_name <- 'generic_PMC_RTSS_EIR_vaccSP_IIV')

  ### Table for text
  cases_df <- fread(file.path(plot_dir, 'csv', 'Fig2CD_dat.csv'))

  table(cases_df$pmc_mode_fct)
  tdf <- cases_df %>%
    filter(age_group == 'U2' &
             pmc_mode_fct != 'None' &
             name %in% c('clinical_cases_averted', 'severe_cases_averted')) %>%
    dplyr::select(Annual_EIR, age_group, name, pmc_mode_fct, mean_val) %>%
    pivot_wider(names_from = pmc_mode_fct, values_from = mean_val) %>%
    mutate(combo_to_rtss = `PMC-3 + RTS,S` / `RTS,S`,
           combo_to_pmc = `PMC-3 + RTS,S` / `PMC-3`)

  tapply(tdf$combo_to_rtss, tdf$name, summary)
  tapply(tdf$combo_to_pmc, tdf$name, summary)


  cases_df <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    filter(Annual_EIR < eir_max) %>%
    mutate(pmc_mode_fct = ifelse(pmc_rtss_cov == '0-0', 'None',
                                 ifelse(pmc_rtss_cov == '0-0.8', 'RTS,S',
                                        ifelse(pmc_rtss_cov == '0.8-0', 'PMC-3', 'PMC-3 + RTS,S'))))

  cases_df$pmc_mode_fct <- factor(cases_df$pmc_mode_fct,
                                  levels = rev(c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')),
                                  labels = rev(c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')))

  cases_df %>%
    ungroup() %>%
    filter(name == 'clinical_cases') %>%
    mutate(Annual_EIR = paste0('eir ', Annual_EIR),
           median = format_num(median_val),
           low = format_num(low_val),
           up = format_num(up_val),
           value = paste0(median, ' (', low, '-', up, ')')) %>%
    dplyr::select(age_group, Annual_EIR, pmc_mode_fct, value) %>%
    pivot_wider(names_from = Annual_EIR, values_from = value)


  cases_df %>%
    ungroup() %>%
    filter(name == 'clinical_cases_averted') %>%
    mutate(Annual_EIR = paste0('eir ', Annual_EIR),
           median = format_num(median_val),
           low = format_num(low_val),
           up = format_num(up_val),
           value = paste0(median, ' (', low, '-', up, ')')) %>%
    dplyr::select(age_group, Annual_EIR, pmc_mode_fct, value) %>%
    pivot_wider(names_from = Annual_EIR, values_from = value)


  cases_df %>%
    ungroup() %>%
    filter(name == 'severe_cases') %>%
    mutate(Annual_EIR = paste0('eir ', Annual_EIR),
           median = format_num(median_val),
           low = format_num(low_val),
           up = format_num(up_val),
           value = paste0(median, ' (', low, '-', up, ')')) %>%
    dplyr::select(age_group, Annual_EIR, pmc_mode_fct, value) %>%
    pivot_wider(names_from = Annual_EIR, values_from = value)

  cases_df %>%
    ungroup() %>%
    filter(name == 'severe_cases_averted') %>%
    mutate(Annual_EIR = paste0('eir ', Annual_EIR),
           median = format_num(median_val),
           low = format_num(low_val),
           up = format_num(up_val),
           value = paste0(median, ' (', low, '-', up, ')')) %>%
    dplyr::select(age_group, Annual_EIR, pmc_mode_fct, value) %>%
    pivot_wider(names_from = Annual_EIR, values_from = value)
}

