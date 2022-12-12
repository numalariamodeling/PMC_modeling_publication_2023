source(file.path('analysis', '_config.R'))

(exp_name <- 'generic_PMCmode_RTSS_vaccSP_IIV')
scenario_labels <- c('None', 'PMC-3', 'PMC-4', 'PMC-5', 'PMC-6', 'PMC-7', 'RTS,S', 'PMC-3 + RTS,S')

outcome_cols <- c('clinical_cases_averted', 'severe_cases_averted')

cases_df_agegrp <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>% filter(age_group %in% c('U1', 'U2'))
cases_df <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_year.csv')) %>%
  mutate(year = ifelse(year >= 2020, year - 2020, year)) %>%
  filter(year < 2) %>%
  mutate(age_group = paste0('>', year)) %>%
  dplyr::select(colnames(cases_df_agegrp)) %>%
  bind_rows(cases_df_agegrp)

pdat <- cases_df %>%
  mutate(pmc_mode = ifelse(pmc_rtss_cov == '0-0.8', 'RTS,S', pmc_mode),
         pmc_mode = ifelse(pmc_rtss_cov == '0-0', 'None', pmc_mode)) %>%
  gen_pmc() %>%
  filter((pmc_rtss_cov == '0.8-0.8' & pmc_mode == '3tp') | (pmc_rtss_cov %in% c('0-0.8', '0.8-0'))) %>%
  filter(!(pmc_mode %in% c('4tp', '5tp'))) %>%
  mutate(yint = ifelse(name == 'clinical_cases', 110, 1.5), yint_min = 0,
         pmc = ifelse(pmc_rtss_cov == '0.8-0.8', 'PMC-3+RTS,S', as.character(pmc)))


pdat <- cases_df %>%
  filter((rtss_coverage == 0.8 & pmc_coverage == 0 | (rtss_coverage == 0 & pmc_coverage != 0))) %>%
  mutate(pmc_mode = ifelse(rtss_coverage == 0.8, 'RTS,S', pmc_mode))

table(pdat$pmc_coverage, pdat$rtss_coverage)
table(pdat$pmc_mode, pdat$rtss_coverage)
table(pdat$name, pdat$age_group)
dat_hline <- pdat %>% filter(rtss_coverage == 0.8 &
                               pmc_coverage == 0 &
                               name %in% outcome_cols &
                               age_group %in% c('>0', '>1'))

pplot <- pdat %>%
  filter(rtss_coverage == 0 &
           name %in% outcome_cols &
           age_group %in% c('>0', '>1')) %>%
  ggplot() +
  geom_hline(yintercept = 0, size = 0.4) +
  geom_hline(data = dat_hline, aes(yintercept = median_val), col = rtss_col) +
  geom_errorbar(aes(x = pmc_coverage, group = pmc_mode, y = median_val, ymin = low_val, ymax = up_val), width = 0, alpha = 0.4) +
  geom_line(aes(x = pmc_coverage, y = median_val, linetype = 'dashed', group = interaction(pmc_mode, rtss_coverage)), col = 'grey') +
  geom_point(aes(x = pmc_coverage, fill = pmc_mode, y = median_val, group = interaction(pmc_mode, rtss_coverage)),
             position = position_dodge(width = 0), shape = 21, col = 'black', size = 3) +
  facet_wrap(name ~ age_group, nrow = 2, scales = 'free_y') +
  scale_x_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2)) +
  scale_color_manual(values = c(pmc_cols)) +
  scale_fill_manual(values = c(pmc_cols)) +
  customTheme_nogrid

print(pplot)

f_save_plot(pplot, paste0('S1_fig_x'),
            file.path(plot_dir), width = 10, height = 8, units = 'in', device_format = device_format)


pdat %>%
  filter(age_group != 'U1' &
           name %in% c('clinical_cases', 'clinical_cases_averted', 'PE_clinical_incidence') &
           age_group != 'U1') %>%
  mutate(median = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(median_val * 100, 1), format_num(median_val)),
         low = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(low_val * 100, 1), format_num(low_val)),
         up = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(up_val * 100, 1), format_num(up_val)),
         value = paste0(median, ' (', low, '-', up, ')')) %>%
  dplyr::select(age_group, pmc_mode, pmc_coverage, name, value) %>%
  mutate(name = gsub('clinical_', '', name)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  pivot_wider(names_from = age_group, values_from = c('cases', 'cases_averted', 'PE_incidence')) %>%
  kbl() %>%
  kable_styling(bootstrap_options = tbl_opts, full_width = F, position = "left", fixed_thead = T)


pdat %>%
  filter(age_group == 'U2' & name %in% c('severe_cases', 'severe_cases_averted', 'PE_severe_incidence')) %>%
  mutate(median = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(median_val * 100, 1), format_num(median_val)),
         low = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(low_val * 100, 1), format_num(low_val)),
         up = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(up_val * 100, 1), format_num(up_val)),
         value = paste0(median, ' (', low, '-', up, ')')) %>%
  dplyr::select(age_group, pmc_mode, pmc_coverage, name, value) %>%
  mutate(name = gsub('severe_', '', name)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  pivot_wider(names_from = age_group, values_from = c('cases', 'cases_averted', 'PE_incidence')) %>%
  kbl() %>%
  kable_styling(bootstrap_options = tbl_opts, full_width = F, position = "left", fixed_thead = T)


###
if (cov_heatmap) {


  exp_name_rtss_cov_scen = 'generic_single_RTSS_vaccSP_IIV'
  dat_agrgrp_rtss_cov <- fread(file.path(simout_dir, exp_name_rtss_cov_scen, 'simdat_aggr_agegroup.csv')) %>% filter(pmc_coverage == 0 & age_group %in% c('U1', 'U2'))
  pdat_rtss <- fread(file.path(simout_dir, exp_name_rtss_cov_scen, 'simdat_aggr_year.csv')) %>%
    filter(pmc_coverage == 0 & rtss_coverage > 0.1) %>%
    mutate(year = ifelse(year >= 2020, year - 2020, year)) %>%
    filter(year < 2) %>%
    mutate(age_group = paste0('>', year)) %>%
    dplyr::select(colnames(dat_agrgrp_rtss_cov)) %>%
    bind_rows(dat_agrgrp_rtss_cov) %>%
    mutate(pmc_mode = 'RTS,S',
           pmc_coverage = rtss_coverage)
  pdat_rtss <- gen_pmc(pdat_rtss)


  pdat_pmc <- pdat %>%
    filter(pmc_mode != 'RTS,S') %>%
    dplyr::select(-rtss_coverage)
  pdat_comb <- pdat_rtss %>%
    bind_rows(pdat_pmc) %>%
    mutate(pmc_coverage = as.numeric(pmc_coverage))

  pplot <- pdat_comb %>%
    filter((rtss_coverage > 0 | pmc_mode != 'RTS,S'),
           name %in% c('PE_clinical_incidence', 'PE_severe_incidence') &
             age_group %in% c('U2')) %>%
    mutate(median_val = median_val * 100,
           pmc_coverage = pmc_coverage * 100,
           name = ifelse(name == 'PE_clinical_incidence', 'Clinical malaria U2', 'Severe malaria U2')) %>%
    ggplot(aes(x = pmc_mode, fill = median_val, y = pmc_coverage, label = round(median_val, 0)),
           position = position_dodge(width = 0)) +
    geom_tile() +
    facet_wrap(~name, nrow = 1) +
    scale_fill_gradientn(colours = rev(RdYlBupalette), limits = c(0, 50)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 0.2), labels = seq(0, 100, 0.2)) +
    scale_x_discrete(expand = c(0, 0), labels = c('3', '4a', '4b', '5a', '5b', '6', '7', 'RTS,S')) +
    labs(y = 'Coverage', x = 'PMC/RTS,S scenario', fill = '% reduction') +
    customTheme_nogrid

  pplot_wtxt <- pplot + geom_text()
  print(pplot_wtxt)

  f_save_plot(pplot, paste0('S1_fig_3C'),
              file.path(plot_dir), width = 7, height = 3, units = 'in', device_format = device_format)

  f_save_plot(pplot_wtxt, paste0('S1_fig_3C_wtxt'),
              file.path(plot_dir), width = 7, height = 3, units = 'in', device_format = device_format)


} #


tdat_a <- pdat %>%
  gen_pmc() %>%
  filter(pmc_coverage == c(0.8), age_group == 'U2' &
    pmc != 'RTS,S' &
    name %in% c('clinical_cases_averted', 'severe_cases_averted') &
    pmc2 != 'PMC-4a' &
    pmc2 != 'PMC-5a') %>%
  arrange(age_group, pmc_rtss_cov, name, pmc2) %>%
  group_by(age_group, pmc_rtss_cov, name) %>%
  mutate(median_diff = median_val - lag(median_val),
         pmc_4_5 = 'a')

exp_key <- c('age_group', 'pmc_rtss_cov', 'name')
tdat_a <- data.table(tdat_a, key = exp_key)
tdat_a[, median_diff_PMC3 := median_val - median_val[pmc == 'PMC-3'], by = exp_key]

tdat_b <- pdat %>%
  gen_pmc() %>%
  filter(pmc_coverage == c(0.8), age_group == 'U2' &
    pmc != 'RTS,S' &
    name %in% c('clinical_cases_averted', 'severe_cases_averted') &
    pmc2 != 'PMC-4b' &
    pmc2 != 'PMC-5b') %>%
  arrange(age_group, pmc_rtss_cov, name, pmc2) %>%
  group_by(age_group, pmc_rtss_cov, name) %>%
  mutate(median_diff = median_val - lag(median_val),
         pmc_4_5 = 'b')

exp_key <- c('age_group', 'pmc_rtss_cov', 'name')
tdat_b <- data.table(tdat_b, key = exp_key)

tdat_b[, median_diff_PMC3 := median_val - median_val[pmc == 'PMC-3'], by = exp_key]

tdat <- tdat_a %>%
  bind_rows(tdat_b) %>%
  mutate(name = gsub('averted', 'averted U2', gsub('_', ' ', name)))

pplot <- ggplot(data = tdat) +
  geom_hline(yintercept = 0) +
  geom_col(aes(x = pmc, y = median_diff, fill = pmc, group = pmc_4_5, linetype = pmc_4_5),
           alpha = 0.7, col = 'black', size = 0.4, width = 0.4, position = position_dodge(w = 0.4), show.legend = F) +
  geom_line(aes(x = pmc, y = median_diff_PMC3, linetype = pmc_4_5, group = pmc_4_5)) +
  facet_wrap(~name, scales = 'free', ncol = 1) +
  labs(y = 'diiference to PMC-3', x = '', linetype = 'Age PMC-4 and 5\ndoses given') +
  scale_linetype_manual(values = c('solid', 'dashed')) +
  scale_color_manual(values = pmc_cols[c(1, 2, 4, 6, 7)]) +
  scale_fill_manual(values = pmc_cols[c(1, 2, 4, 6, 7)]) +
  customTheme_nogrid


print(pplot)

if (SAVE) {
  f_save_plot(pplot, paste0('fig_SI_X'),
              file.path(plot_dir), width = 8, height = 5, units = 'in', device_format = device_format)

}


rtss_fine_age = F
if (rtss_fine_age) {
  (exp_name <- 'generic_PMCmode_RTSS_vaccSP_IIV')


  dat_aggr <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_week.csv'))
  dat_aggr$pmc_mode_fct <- factor(dat_aggr$pmc_mode,
                                  levels = c('3tp', '5tp', '7tp2ndyr'),
                                  labels = c('PMC-3', 'PMC-5', 'PMC-7'))

  dat_aggr %>%
    group_by(Annual_EIR, pmc_rtss_cov) %>%
    unique() %>%
    tally() %>%
    pivot_wider(names_from = pmc_rtss_cov, values_from = n)

  y <- 2
  pdat1 <- subset(dat_aggr, age <= y * 365 &
    name == 'clinical_cases' &
    pmc_coverage == 0.8 &
    pmc_mode %in% c('3tp', '5tp', '7tp2ndyr'))

  pdat2 <- dat_aggr %>%
    filter(rtss_coverage == 0 &
             pmc_coverage == 0 &
             age <= y * 365 &
             name == 'clinical_cases') %>%
    dplyr::select(-pmc_mode, -pmc_mode_fct, -rtss_coverage, -pmc_coverage)

  pdat3 <- subset(dat_aggr, age <= y * 365 &
    name == 'clinical_cases' &
    rtss_coverage != 0 &
    pmc_coverage == 0) %>%
    dplyr::select(-pmc_mode, -pmc_mode_fct)


  pplot <- ggplot(data = pdat1) +
    geom_vline(xintercept = c(365)) +
    geom_line(data = pdat2,
              aes(x = age, y = median_val), col = 'darkgrey') +
    geom_ribbon(data = pdat2,
                aes(x = age, ymin = low_val, ymax = up_val),
                alpha = 0.4, fill = 'darkgrey') +
    geom_ribbon(aes(x = age, ymin = low_val, ymax = up_val,
                    group = interaction(pmc_mode_fct, rtss_coverage),
                    fill = pmc_mode_fct), alpha = 0.4) +
    geom_ribbon(data = pdat3,
                aes(x = age, ymin = low_val, ymax = up_val), alpha = 0.3, fill = rtss_col) +
    geom_line(aes(x = age, y = median_val, group = interaction(pmc_mode_fct, rtss_coverage),
                  col = pmc_mode_fct, linetype = as.factor(rtss_coverage)), size = 1) +
    geom_line(data = pdat3,
              aes(x = age, y = median_val),
              col = rtss_col, size = 1) +
    facet_wrap(~pmc_mode_fct, nrow = 1, scales = 'free_x') +
    scale_linetype_manual(values = c('solid', 'dashed')) +
    scale_color_manual(values = c(pmc_cols[c(1, 2, 5)], rtss_col)) +
    scale_fill_manual(values = c(pmc_cols[c(1, 2, 5)], rtss_col)) +
    scale_x_continuous(lim = c(0, y * 365),
                       breaks = seq(0, y * 365, 30 * 6),
                       labels = seq(0, y * 12, 6) / 12) +
    labs(col = '', fill = '', x = 'Age (years)',
         y = 'Weekly cases \nper 1000 population') +
    theme(panel.grid.major = element_blank(),
          legend.position = 'None') +
    customTheme_nogrid

  print(pplot)
  f_save_plot(pplot, paste0('S1_fig_x'),
              file.path(plot_dir), width = 8, height = 3.5, units = 'in', device_format = device_format)


}


## WIth RTSS
if (with_rtss_heatmap2) {
  exp_name <- 'generic_PMCmode_RTSS_cov_vaccSP_IIV'

  dat_agrgrp <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
    filter(age_group %in% c('U1', 'U2'))
  dat_aggr <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_year.csv')) %>%
    filter(year < 2) %>%
    mutate(age_group = paste0('>', year)) %>%
    dplyr::select(colnames(dat_agrgrp)) %>%
    bind_rows(dat_agrgrp)

  dat_aggr$pmc_mode_fct <- factor(dat_aggr$pmc_mode,
                                  levels = c('3tp', '5tp', '7tp2ndyr'),
                                  labels = c('PMC-3', 'PMC-5', 'PMC-7'))

  table(dat_aggr$pmc_mode_fct, dat_aggr$pmc_rtss_cov, exclude = NULL)

  #dat_table(dat_agrgrp, qc(age_group, pmc_coverage, rtss_coverage, name))

  pcols <- c('clinical_cases_averted',
             'severe_cases_averted',
             'PE_clinical_incidence',
             'PE_severe_incidence')

  pdat <- dat_aggr %>% filter(name %in% pcols &
                                age_group == 'U2' &
                                pmc_coverage != 0)
  pdat_c <- dat_aggr %>%
    filter(name %in% pcols &
             age_group == 'U2' &
             pmc_coverage == 0) %>%
    select(-pmc_mode, -pmc_mode_fct) %>%
    left_join(pdat[, c('age_group', 'name', 'pmc_mode', 'pmc_mode_fct')])

  pdat <- pdat %>% bind_rows(pdat_c)

  pplot1 <- ggplot(data = subset(pdat, name == 'PE_clinical_incidence')) +
    geom_tile(aes(x = as.factor(pmc_coverage * 100), y = as.factor(rtss_coverage * 100), fill = median_val * 100)) +
    scale_fill_gradientn(colours = rev(RdYlBupalette), limits = c(0, 80)) +
    facet_wrap(name ~ pmc_mode_fct) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    customTheme_nogrid +
    labs(x = 'PMC coverage (%)', y = 'RTS,S coverage (%)', fill = '% reduction')

  pplot2 <- ggplot(data = subset(pdat, name == 'PE_severe_incidence')) +
    geom_tile(aes(x = as.factor(pmc_coverage * 100), y = as.factor(rtss_coverage * 100), fill = median_val * 100)) +
    scale_fill_gradientn(colours = rev(RdYlBupalette), limits = c(0, 80)) +
    facet_wrap(name ~ pmc_mode_fct) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    customTheme_nogrid +
    labs(x = 'PMC coverage (%)', y = 'RTS,S coverage (%)', fill = '% reduction')

  pplot <- plot_combine(list(pplot1, pplot2), ncol = 1)

  print(pplot)


  f_save_plot(pplot, paste0('S1_fig_heatmap'),
              file.path(plot_dir), width = 8, height = 6, units = 'in', device_format = device_format)

  #### tables

  pdat <- pdat %>%
    group_by(ipti_coverage, rtss_coverage, name, age_group) %>%
    summarize(median_val_low = min(median_val),
              median_val_mean = mean(median_val),
              median_val_up = max(median_val))

  pdat %>%
    ggplot() +
    geom_hline(data = pdat_c, aes(yintercept = median_val), col = rtss_col) +
    geom_line(aes(x = as.factor(ipti_coverage), y = median_val_mean, group = rtss_coverage), col = 'grey', size = 1, linetype = 'dashed') +
    #geom_errorbar(aes(x = as.factor(ipti_coverage), ymin = median_val_low, ymax = median_val_up, group = as.factor(rtss_coverage)), width = 0.01) +
    geom_point(aes(x = as.factor(ipti_coverage), y = median_val_mean, fill = as.factor(rtss_coverage)), size = 3, shape = 21) +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    facet_wrap(~name, scales = 'free_y') +
    customTheme_nogrid


}

