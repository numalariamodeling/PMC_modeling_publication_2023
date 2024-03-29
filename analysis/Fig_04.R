##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## Fig_04.R
##---------------------

source(file.path('analysis', '_config.R'))
(exp_name <- 'generic_PMCmode_RTSS_vaccSP_IIV')

fig4A <- T
fig4B <- T
fig4B_cum <- T
fig4B_PE <- T
fig4C <- T
fig4C_cum <- T
result_tables <- T

y <- 2.0

cases_df <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_week.csv')) %>%
  rename_with(~gsub('ipti', 'pmc', .x))
table(cases_df$pmc_coverage, cases_df$rtss_coverage, exclude = NULL)

if (fig4A) {
  pdat1 <- cases_df %>%
    mutate(pmc_mode = ifelse(pmc_rtss_cov == '0-0.8', 'RTS,S', pmc_mode),
           pmc_mode = ifelse(pmc_rtss_cov == '0-0', 'None', pmc_mode)) %>%
    gen_pmc_mode_fct()

  pdat1 <- pdat1 %>%
    filter(age <= y * 365 & name %in% c('clinical_cases', 'severe_cases', 'PE_clinical_incidence', 'PE_severe_incidence')) %>%
    filter((pmc_rtss_cov == '0.8-0.8' & pmc_mode == '3tp') | (pmc_rtss_cov %in% c('0-0.8', '0.8-0'))) %>%
    filter(!(pmc_mode %in% c('4tp', '5tp'))) %>%
    mutate(yint = ifelse(name == 'clinical_cases', 110, 1.5), yint_min = 0,
           pmc_mode_fct = ifelse(pmc_rtss_cov == '0.8-0.8', 'PMC-3+RTS,S', as.character(pmc_mode_fct)))

  scenario_labels <- c('PMC-3', 'PMC-4', 'PMC-4', 'PMC-5', 'PMC-5', 'PMC-6', 'PMC-7', 'RTS,S', 'PMC-3+RTS,S')
  pdat1$pmc_mode_fct <- factor(pdat1$pmc_mode_fct, levels = scenario_labels, labels = scenario_labels)

  table(pdat1$pmc_mode, pdat1$pmc_mode_fct, exclude = NULL)
  table(pdat1$rtss_coverage, pdat1$pmc_mode_fct, exclude = NULL)

  pdat2 <- cases_df %>%
    dplyr::filter(rtss_coverage == 0 & pmc_coverage == 0) %>%
    dplyr::select(-pmc_mode, -rtss_coverage, -pmc_coverage)
  pdat2 <- subset(pdat2, age <= y * 365 & name %in% c('clinical_cases', 'severe_cases',
                                                      'PE_clinical_incidence', 'PE_severe_incidence'))

  out_vars <- c('clinical_cases', 'severe_cases')
  pplot <- ggplot(data = subset(pdat1, name %in% out_vars)) +
    geom_vline(xintercept = c(365)) +
    geom_hline(aes(yintercept = yint), size = 0.1, col = 'white', alpha = 0.96) +
    geom_hline(aes(yintercept = yint_min), size = 0.1, col = 'white', alpha = 0.96) +
    geom_line(data = subset(pdat2, name %in% out_vars),
              aes(x = age, y = median_val), col = 'darkgrey') +
    geom_ribbon(data = subset(pdat2, name %in% out_vars),
                aes(x = age, ymin = low_val, ymax = up_val, fill = pmc_mode_fct),
                alpha = 0.4, fill = 'darkgrey') +
    geom_ribbon(data = subset(pdat1, name %in% out_vars),
                aes(x = age, ymin = low_val, ymax = up_val, fill = pmc_mode_fct, group = pmc_mode_fct),
                alpha = 0.4) +
    geom_line(data = subset(pdat1, name %in% out_vars),
              aes(x = age, y = median_val, col = pmc_mode_fct)) +
    facet_wrap(pmc_mode_fct ~ name, ncol = 2, scales = 'free_y', strip.position = "right") +
    scale_color_manual(values = custom_cols) +
    scale_fill_manual(values = custom_cols) +
    scale_x_continuous(lim = c(0, y * 365),
                       breaks = seq(0, y * 365, 30 * 6),
                       labels = seq(0, y * 12, 6) / 12) +
    customTheme +
    labs(col = '', fill = '', x = 'Age (years)',
         y = 'Cases \nper 1000 population') +
    theme(panel.grid.major = element_blank(),
          legend.position = 'None')

  print(pplot)

  f_save_plot(pplot, paste0('Fig4A'),
              file.path(plot_dir), width = 8, height = 10, units = 'in', device_format = device_format)
  fwrite(pdat1, file.path(plot_dir, 'csv', 'Fig4A_dat.csv'))

  #### PE
  out_vars <- c('PE_clinical_incidence', 'PE_severe_incidence')
  pplot <- ggplot(data = subset(pdat1, name %in% out_vars)) +
    geom_vline(xintercept = c(365)) +
    geom_hline(aes(yintercept = yint), size = 0.1, col = 'white', alpha = 0.96) +
    geom_hline(aes(yintercept = yint_min), size = 0.1, col = 'white', alpha = 0.96) +
    geom_line(data = subset(pdat2, name %in% out_vars),
              aes(x = age, y = median_val), col = 'darkgrey') +
    geom_ribbon(data = subset(pdat2, name %in% out_vars),
                aes(x = age, ymin = low_val, ymax = up_val, fill = pmc_mode_fct),
                alpha = 0.4, fill = 'darkgrey') +
    geom_ribbon(data = subset(pdat1, name %in% out_vars),
                aes(x = age, ymin = low_val, ymax = up_val, fill = pmc_mode_fct, group = pmc_mode_fct),
                alpha = 0.4) +
    geom_line(data = subset(pdat1, name %in% out_vars),
              aes(x = age, y = median_val, col = pmc_mode_fct)) +
    facet_wrap(pmc_mode_fct ~ name, ncol = 2, scales = 'free_y', strip.position = "right") +
    scale_color_manual(values = custom_cols) +
    scale_fill_manual(values = custom_cols) +
    scale_x_continuous(lim = c(0, y * 365),
                       breaks = seq(0, y * 365, 30 * 6),
                       labels = seq(0, y * 12, 6) / 12) +
    customTheme +
    labs(col = '', fill = '', x = 'Age (years)',
         y = '% reduction in cases') +
    theme(panel.grid.major = element_blank(),
          legend.position = 'None')


  print(pplot)
  f_save_plot(pplot, paste0('Fig4A_PE'),
              file.path(plot_dir), width = 8, height = 10, units = 'in', device_format = device_format)
  fwrite(pdat1, file.path(plot_dir, 'csv', 'Fig4A_PE_dat.csv'))

} #fig4A

if (fig4B) {

  outcome_cols <- c('clinical_cases_averted', 'severe_cases_averted')

  cases_df_agegrp <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    filter(age_group %in% c('U1', 'U2'))

  cases_df <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_year.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    mutate(year = ifelse(year >= 2020, year - 2020, year)) %>%
    filter(year < 2) %>%
    mutate(age_group = paste0('>', year)) %>%
    dplyr::select(colnames(cases_df_agegrp)) %>%
    bind_rows(cases_df_agegrp) %>%
    as.data.table()

  pdat <- cases_df %>%
    mutate(pmc_mode = ifelse(pmc_rtss_cov == '0-0.8', 'RTS,S', pmc_mode),
           pmc_mode = ifelse(pmc_rtss_cov == '0-0', 'None', pmc_mode)) %>%
    gen_pmc_mode_fct() %>%
    dplyr::filter((pmc_rtss_cov == '0.8-0.8' & pmc_mode == '3tp') | (pmc_rtss_cov %in% c('0-0.8', '0.8-0'))) %>%
    dplyr::filter(!(pmc_mode %in% c('4tp', '5tp'))) %>%
    dplyr::mutate(yint = ifelse(name == 'clinical_cases', 110, 1.5), yint_min = 0,
                  pmc_mode_fct = ifelse(pmc_rtss_cov == '0.8-0.8', 'PMC-3 + RTS,S', as.character(pmc_mode_fct)))

  pdat$pmc_mode_fct <- factor(pdat$pmc_mode_fct,
                              levels = scenario_labels,
                              labels = gsub('[+]', '\n+', scenario_labels))

  pdat <- pdat %>%
    mutate(yint = ifelse(name == 'clinical_cases_averted', 1500, 25),
           yint_min = ifelse(name == 'clinical_cases_averted', -50, -5))
  table(pdat$age_group)
  pdat <- subset(pdat, age_group != 'U1' &
    pmc_mode_fct != 'None' &
    name %in% outcome_cols)
  pdat$age_group <- factor(pdat$age_group, levels = c('>0', '>1', 'U2'), labels = c('0-1', '1-2', '0-2'))
  pdat$name <- factor(pdat$name,
                      levels = c('clinical_cases_averted', 'severe_cases_averted'),
                      labels = c('Clinical malaria in children U2 (0-2 years)', 'Severe malaria in children U2 (0-2 years)'))


  pplot <- ggplot(data = pdat) +
    geom_hline(yintercept = 0, size = 0.2) +
    geom_hline(aes(yintercept = yint), alpha = 0) +
    geom_hline(aes(yintercept = yint_min), alpha = 0) +
    geom_col(aes(x = pmc_mode_fct, y = median_val, alpha = age_group, group = interaction(age_group, pmc_mode_fct), fill = pmc_mode_fct),
             position = position_dodge(width = 0.75), size = 0.3, width = 0.7, col = 'black') +
    geom_errorbar(aes(x = pmc_mode_fct, y = median_val, ymin = low_val, ymax = up_val, group = interaction(age_group, pmc_mode_fct)),
                  position = position_dodge(width = 0.75), width = 0.01) +
    facet_wrap(~name, nrow = 2, scales = 'free_y') +
    scale_y_continuous(labels = comma, expand = c(0, 0)) +
    theme(panel.spacing = unit(1, "lines"),
          axis.ticks.x = element_blank()) +
    labs(x = '',
         y = 'Annual cases averted\n per 1000 population',
         alpha = 'Age group\n(years)',
         fill = 'Scenario',
         col = 'Scenario') +
    #scale_x_discrete(labels = c(2:7)) +
    scale_alpha_manual(values = c(0.5, 0.75, 1)) +
    scale_fill_manual(values = custom_cols) +
    customTheme_nogrid

  print(pplot)
  f_save_plot(pplot, paste0('Fig4B'),
              file.path(plot_dir), width = 7, height = 6, units = 'in', device_format = device_format)
  fwrite(pdat, file.path(plot_dir, 'csv', 'Fig4B_dat.csv'))

  rm(pplot, cases_df, cases_df_agegrp)
} #fig4B

if (fig4B_PE) {

  outcome_cols <- c('PE_clinical_incidence', 'PE_severe_incidence')

  cases_df_agegrp <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    filter(age_group %in% c('U1', 'U2'))

  cases_df <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_year.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    mutate(year = ifelse(year >= 2020, year - 2020, year)) %>%
    filter(year < 2) %>%
    mutate(age_group = paste0('>', year)) %>%
    dplyr::select(colnames(cases_df_agegrp)) %>%
    bind_rows(cases_df_agegrp) %>%
    as.data.table()

  pdat <- cases_df %>%
    mutate(pmc_mode = ifelse(pmc_rtss_cov == '0-0.8', 'RTS,S', pmc_mode),
           pmc_mode = ifelse(pmc_rtss_cov == '0-0', 'None', pmc_mode)) %>%
    gen_pmc_mode_fct() %>%
    dplyr::filter((pmc_rtss_cov == '0.8-0.8' & pmc_mode == '3tp') | (pmc_rtss_cov %in% c('0-0.8', '0.8-0'))) %>%
    dplyr::filter(!(pmc_mode %in% c('4tp', '5tp'))) %>%
    dplyr::mutate(pmc_mode_fct = ifelse(pmc_rtss_cov == '0.8-0.8', 'PMC-3 + RTS,S', as.character(pmc_mode_fct)))

  pdat$pmc_mode_fct <- factor(pdat$pmc_mode_fct,
                              levels = scenario_labels,
                              labels = gsub('[+]', '\n+', scenario_labels))

  table(pdat$age_group)
  pdat <- subset(pdat, age_group != 'U1' &
    pmc_mode_fct != 'None' &
    name %in% outcome_cols)
  pdat$age_group <- factor(pdat$age_group, levels = c('>0', '>1', 'U2'), labels = c('0-1', '1-2', '0-2'))
  pdat$name <- factor(pdat$name,
                      levels = c('PE_clinical_incidence', 'PE_severe_incidence'),
                      labels = c('Clinical malaria by age in children U2 (0-2 years)',
                                 'Severe malaria in children U2 (0-2 years)'))


  pplot <- ggplot(data = pdat) +
    geom_hline(yintercept = c(0.45), alpha = 0) +
    geom_hline(yintercept = c(0), size = 0.2) +
    geom_col(aes(x = pmc_mode_fct, y = median_val, alpha = age_group,
                 group = interaction(age_group, pmc_mode_fct), fill = pmc_mode_fct),
             position = position_dodge(width = 0.75), size = 0.3, width = 0.7, col = 'black') +
    geom_errorbar(aes(x = pmc_mode_fct, y = median_val, ymin = low_val, ymax = up_val,
                      group = interaction(age_group, pmc_mode_fct)),
                  position = position_dodge(width = 0.75), width = 0.01) +
    facet_wrap(~name, nrow = 2, scales = 'free_y') +
    scale_y_continuous(expand = c(0, 0),
                       breaks = seq(-0.1, 0.5, 0.1),
                       labels = seq(-0.1, 0.5, 0.1) * 100) +
    theme(panel.spacing = unit(1, "lines"),
          axis.ticks.x = element_blank()) +
    labs(x = '',
         y = '% reduction in cases',
         alpha = 'Age group\n(years)',
         fill = 'Scenario',
         col = 'Scenario') +
    #scale_x_discrete(labels = c(2:7)) +
    scale_alpha_manual(values = c(0.5, 0.75, 1)) +
    scale_fill_manual(values = custom_cols) +
    customTheme_nogrid

  print(pplot)
  f_save_plot(pplot, paste0('Fig4B_PE'),
              file.path(plot_dir), width = 7, height = 6, units = 'in', device_format = device_format)
  fwrite(pdat, file.path(plot_dir, 'csv', 'Fig4B_PE_dat.csv'))

  rm(pplot, cases_df, cases_df_agegrp)
} #fig4B_PE

if (fig4B_cum) {

  outcome_cols <- c('clinical_cases_averted', 'severe_cases_averted')

  cases_df_agegrp <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    filter(age_group %in% c('U1', 'U2')) %>%
    filter(name %in% paste0(outcome_cols, '_cum')) %>%
    mutate(name = gsub('_cum', '', name))

  cases_df <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_year.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    mutate(year = ifelse(year >= 2020, year - 2020, year)) %>%
    filter(year < 2) %>%
    mutate(age_group = paste0('>', year)) %>%
    dplyr::select(colnames(cases_df_agegrp)) %>%
    bind_rows(cases_df_agegrp) %>%
    as.data.table()

  pdat <- cases_df %>%
    mutate(pmc_mode = ifelse(pmc_rtss_cov == '0-0.8', 'RTS,S', pmc_mode),
           pmc_mode = ifelse(pmc_rtss_cov == '0-0', 'None', pmc_mode)) %>%
    gen_pmc_mode_fct() %>%
    dplyr::filter((pmc_rtss_cov == '0.8-0.8' & pmc_mode == '3tp') | (pmc_rtss_cov %in% c('0-0.8', '0.8-0'))) %>%
    dplyr::filter(!(pmc_mode %in% c('4tp', '5tp'))) %>%
    dplyr::mutate(yint = ifelse(name == 'clinical_cases', 110, 1.5), yint_min = 0,
                  pmc_mode_fct = ifelse(pmc_rtss_cov == '0.8-0.8', 'PMC-3 + RTS,S', as.character(pmc_mode_fct)))

  pdat$pmc_mode_fct <- factor(pdat$pmc_mode_fct,
                              levels = scenario_labels,
                              labels = scenario_labels)

  pdat <- pdat %>%
    mutate(yint = ifelse(name == 'clinical_cases_averted', 2250, 35), yint_min = 0)

  pdatU1 <- pdat %>%
    filter(age_group == 'U1' & name %in% outcome_cols)

  w <- 0.7
  pplot <- ggplot(data = subset(pdat, age_group %in% c('>0', '>1') &
    pmc_mode_fct != 'None' &
    name %in% outcome_cols)) +
    geom_hline(yintercept = 0, size = 0.2) +
    geom_hline(aes(yintercept = yint), size = 0.1, col = 'white', alpha = 0.96) +
    geom_hline(aes(yintercept = yint_min), size = 0.1, col = 'white', alpha = 0.96) +
    geom_col(data = subset(pdat, age_group == 'U2' &
      pmc_mode_fct != 'None' &
      name %in% outcome_cols),
             aes(x = pmc_mode_fct, y = median_val, fill = pmc_mode_fct),
             position = position_dodge(), size = 0.3, width = w - 0.1) +
    geom_col(aes(x = pmc_mode_fct, y = median_val, linetype = age_group,
                 group = interaction(age_group, pmc_mode_fct)),
             fill = 'white', col = 'black', position = position_dodge(), alpha = 0.7, size = 0.3, width = w) +
    geom_errorbar(aes(x = pmc_mode_fct, y = median_val, ymin = low_val, ymax = up_val,
                      group = interaction(age_group, pmc_mode_fct)),
                  position = position_dodge(width = w), width = 0.01) +
    geom_errorbar(data = subset(pdat, age_group == 'U2' &
      pmc_mode_fct != 'None' &
      name %in% outcome_cols),
                  aes(x = pmc_mode_fct, y = median_val, ymin = low_val, ymax = up_val, group = pmc_mode_fct),
                  position = position_dodge(), width = 0.01) +
    facet_wrap(~name, nrow = 2, scales = 'free') +
    scale_y_continuous(labels = comma, expand = c(0, 0)) +
    theme(panel.spacing = unit(1, "lines"),
          axis.ticks.x = element_blank()) +
    labs(x = 'Scenario',
         y = 'Cumulative cases averted\n per 1000 population',
         fill = 'Scenario', col = 'Scenario') +
    scale_fill_manual(values = custom_cols) +
    customTheme_nogrid


  print(pplot)
  f_save_plot(pplot, paste0('Fig4B_cum'),
              file.path(plot_dir), width = 7, height = 6, units = 'in', device_format = device_format)
  fwrite(pdat, file.path(plot_dir, 'csv', 'Fig4B_cum_dat.csv'))

  rm(pplot, cases_df, cases_df_agegrp)
} #fig4B_cum

if (fig4B_PE) {
  outcome_cols <- c('PE_clinical_incidence', 'PE_severe_incidence')

  cases_df_agegrp <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    filter(age_group %in% c('U1', 'U2'))

  cases_df <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_year.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    mutate(year = ifelse(year >= 2020, year - 2020, year)) %>%
    filter(year < 2) %>%
    mutate(age_group = paste0('>', year)) %>%
    dplyr::select(colnames(cases_df_agegrp)) %>%
    bind_rows(cases_df_agegrp) %>%
    as.data.table()

  pdat <- cases_df %>%
    mutate(pmc_mode = ifelse(pmc_rtss_cov == '0-0.8', 'RTS,S', pmc_mode),
           pmc_mode = ifelse(pmc_rtss_cov == '0-0', 'None', pmc_mode)) %>%
    gen_pmc_mode_fct() %>%
    dplyr::filter((pmc_rtss_cov == '0.8-0.8' & pmc_mode == '3tp') | (pmc_rtss_cov %in% c('0-0.8', '0.8-0'))) %>%
    dplyr::filter(!(pmc_mode %in% c('4tp', '5tp'))) %>%
    dplyr::mutate(yint = ifelse(name == 'clinical_cases', 110, 1.5), yint_min = 0,
                  pmc_mode_fct = ifelse(pmc_rtss_cov == '0.8-0.8', 'PMC-3 + RTS,S', as.character(pmc_mode_fct)))


  pdat$pmc_mode_fct <- factor(pdat$pmc_mode_fct,
                              levels = scenario_labels,
                              labels = scenario_labels)

  pdat <- pdat %>%
    mutate(yint = 1, yint_min = -0.1)

  pdatU1 <- pdat %>%
    filter(age_group == 'U1' & name %in% outcome_cols)

  pplot <- ggplot(data = subset(pdat, age_group %in% c('>0', '>1') &
    pmc_mode_fct != 'None' &
    name %in% outcome_cols)) +
    geom_hline(yintercept = 0, size = 0.2) +
    geom_hline(aes(yintercept = yint), size = 0.1, col = 'white', alpha = 0.96) +
    geom_hline(aes(yintercept = yint_min), size = 0.1, col = 'white', alpha = 0.96) +
    geom_col(data = subset(pdat, age_group == 'U2' &
      pmc_mode_fct != 'None' &
      name %in% outcome_cols),
             aes(x = pmc_mode_fct, y = median_val, fill = pmc_mode_fct),
             position = position_dodge(), size = 0.3, width = w - 0.1) +
    geom_col(aes(x = pmc_mode_fct, y = median_val, linetype = age_group,
                 group = interaction(age_group, pmc_mode_fct)),
             fill = 'white', col = 'black', position = position_dodge(), alpha = 0.7, size = 0.3, width = w) +
    geom_errorbar(aes(x = pmc_mode_fct, y = median_val, ymin = low_val, ymax = up_val,
                      group = interaction(age_group, pmc_mode_fct)),
                  position = position_dodge(width = w), width = 0.01) +
    geom_errorbar(data = subset(pdat, age_group == 'U2' &
      pmc_mode_fct != 'None' &
      name %in% outcome_cols),
                  aes(x = pmc_mode_fct, y = median_val, ymin = low_val, ymax = up_val, group = pmc_mode_fct),
                  position = position_dodge(), width = 0.01) +
    facet_wrap(~name, nrow = 2, scales = 'free') +
    scale_y_continuous(labels = comma, expand = c(0, 0)) +
    theme(panel.spacing = unit(1, "lines"),
          axis.ticks.x = element_blank()) +
    labs(x = 'Scenario',
         y = '% reduction',
         fill = 'Scenario', col = 'Scenario') +
    scale_fill_manual(values = custom_cols) +
    customTheme_nogrid

  print(pplot)

  f_save_plot(pplot, paste0('Fig4B_PE'),
              file.path(plot_dir), width = 7, height = 6, units = 'in', device_format = device_format)
  fwrite(pdat, file.path(plot_dir, 'csv', 'Fig4B_PE_dat.csv'))

} #fig4B_PE

if (fig4C) {

  outcome_cols <- c('clinical_cases_averted', 'severe_cases_averted')
  cases_df_agegrp <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    filter(age_group %in% c('U2')) %>%
    filter(name %in% outcome_cols)

  pdat <- cases_df_agegrp %>%
    filter(pmc_rtss_cov %in% c('0-0', '0.8-0', '0-0.8', '0.8-0.8')) %>%
    filter(pmc_mode %in% c('3tp', '5tp', '7tp2ndyr'))

  table(pdat$name, exclude = NULL)
  table(pdat$pmc_rtss_cov, pdat$pmc_mode, exclude = NULL)
  pdat <- gen_pmc_mode_fct(pdat)

  pdat <- pdat %>% filter(pmc_rtss_cov != '0-0')

  pdat_rtss <- subset(pdat, rtss_coverage != 0 & pmc_coverage == 0) %>%
    dplyr::select(-pmc_mode, -pmc_mode_fct) %>%
    left_join(unique(pdat[, c('age_group', 'pmc_mode', 'pmc_mode_fct')]))

  rtss_alone <- pdat %>%
    filter(pmc_coverage == 0 &
             rtss_coverage != 0) %>%
    dplyr::select(-pmc_mode, -pmc_mode_fct, -pmc_rtss_cov)

  rtss_alone <- as.data.frame(unique(pdat[, c('pmc_mode_fct', 'pmc_mode')])) %>%
    left_join(rtss_alone, by = character())

  pdat <- pdat %>%
    mutate(yint = ifelse(name == 'clinical_cases_averted', 1700, 20), yint_min = 0)

  pdat1 <- subset(pdat, pmc_coverage != 0 & rtss_coverage != 0)
  #### RTSS in addition to PMC
  pplot1 <- ggplot(data = pdat1) +
    geom_hline(aes(yintercept = yint), alpha = 0) +
    geom_col(aes(x = pmc_mode_fct, y = median_val, col = pmc_mode_fct),
             fill = 'NA', position = position_dodge(width = 0.9), linetype = 'dashed',
             size = 1, width = 0.8, show.legend = F, alpha = 0.6) +
    geom_col(data = subset(pdat, ((pmc_coverage != 0 & rtss_coverage == 0))),
             aes(x = pmc_mode_fct, y = median_val, fill = pmc_mode_fct),
             position = position_dodge(width = 0.9), size = 0.3, width = 0.8) +
    geom_errorbar(aes(x = pmc_mode_fct, y = median_val, ymin = low_val, ymax = up_val), width = 0.01) +
    facet_wrap(~name, nrow = 2, scales = 'free_y') +
    scale_y_continuous(labels = comma, expand = c(0, 0)) +
    labs(subtitle = 'RTS,S in addition to PMC\n',
         x = 'Age (months)',
         y = 'Cumulative cases averted per\n1000 population in children U2',
         color = '', fill = '') +
    scale_color_manual(values = c(pmc_cols[c(1, 2, 5)], rtss_col)) +
    scale_fill_manual(values = c(pmc_cols[c(1, 2, 5)], rtss_col)) +
    customTheme_nogrid +
    theme(legend.position = 'None')


  rtss_alone <- rtss_alone %>%
    dplyr::select(pmc_mode_fct, Annual_EIR, age_group, cm_coverage, name, median_val, low_val, up_val) %>%
    rename(median_val_rtss = median_val,
           low_val_rtss = low_val,
           up_val_rtss = up_val)

  pdat2 <- subset(pdat, pmc_coverage == 0.8 & rtss_coverage == 0.8) %>%
    left_join(rtss_alone) %>%
    mutate(median_val_subrtss = median_val - median_val_rtss,
           low_val_subrtss = low_val - low_val_rtss,
           up_val_subrtss = up_val - up_val_rtss)

  pplot2 <- ggplot(data = pdat2) +
    geom_hline(aes(yintercept = yint), alpha = 0) +
    geom_col(aes(x = pmc_mode_fct, y = median_val, col = pmc_mode_fct), position = position_dodge(width = 0.9),
             fill = 'NA', linetype = 'dashed', size = 1, width = 0.8, show.legend = F, alpha = 0.6) +
    geom_col(aes(x = pmc_mode_fct, y = median_val_subrtss, fill = pmc_mode_fct),
             position = position_dodge(width = 0.9), size = 1, width = 0.8, show.legend = F) +
    geom_errorbar(aes(x = pmc_mode_fct, y = median_val, ymin = low_val, ymax = up_val), width = 0.01) +
    facet_wrap(~name, nrow = 2, scales = 'free_y') +
    scale_y_continuous(labels = comma, expand = c(0, 0)) +
    labs(subtitle = 'PMC in addition to RTS,S\n',
         x = 'Age (months)',
         y = 'Cumulative cases averted per\n1000 population in children U2',
         color = '', fill = '') +
    scale_color_manual(values = c(pmc_cols[c(1, 2, 5)], rtss_col)) +
    scale_fill_manual(values = c(pmc_cols[c(1, 2, 5)], rtss_col)) +
    customTheme_nogrid +
    theme(legend.position = 'None')

  pplot <- plot_combine(list(pplot1, pplot2), ncol = 2)
  print(pplot)

  f_save_plot(pplot, paste0('Fig4C'),
              file.path(plot_dir), width = 10, height = 5, units = 'in', device_format = device_format)
  fwrite(pdat, file.path(plot_dir, 'csv', 'Fig4C_dat_.csv'))
  fwrite(pdat1, file.path(plot_dir, 'csv', 'Fig4CA_dat_.csv'))
  fwrite(pdat2, file.path(plot_dir, 'csv', 'Fig4CB_dat_.csv'))

  pdat %>%
    filter(name == 'clinical_cases_averted') %>%
    dplyr::select(age_group, pmc_mode_fct, pmc_rtss_cov, median_val)
  rtss_alone %>%
    filter(name == 'clinical_cases_averted') %>%
    dplyr::select(age_group, pmc_mode_fct, median_val_rtss)


} #fig4C

if (fig4C_cum) {

  outcome_cols <- c('clinical_cases_averted', 'severe_cases_averted')
  cases_df_agegrp <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    filter(age_group %in% c('U2')) %>%
    filter(name %in% paste0(outcome_cols, '_cum')) %>%
    mutate(name = gsub('_cum', '', name))

  pdat <- cases_df_agegrp %>%
    filter(pmc_rtss_cov %in% c('0-0', '0.8-0', '0-0.8', '0.8-0.8')) %>%
    filter(pmc_mode %in% c('3tp', '5tp', '7tp2ndyr'))

  table(pdat$name, exclude = NULL)
  table(pdat$pmc_rtss_cov, pdat$pmc_mode, exclude = NULL)
  pdat <- gen_pmc_mode_fct(pdat)

  pdat <- pdat %>% filter(pmc_rtss_cov != '0-0')

  pdat_rtss <- subset(pdat, rtss_coverage != 0 & pmc_coverage == 0) %>%
    dplyr::select(-pmc_mode, -pmc_mode_fct) %>%
    left_join(unique(pdat[, c('age_group', 'pmc_mode', 'pmc_mode_fct')]))

  rtss_alone <- pdat %>%
    filter(pmc_coverage == 0 &
             rtss_coverage != 0) %>%
    dplyr::select(-pmc_mode, -pmc_mode_fct, -pmc_rtss_cov)

  rtss_alone <- as.data.frame(unique(pdat[, c('pmc_mode_fct', 'pmc_mode')])) %>%
    left_join(rtss_alone, by = character())

  pdat <- pdat %>%
    mutate(yint = ifelse(name == 'clinical_cases_averted', 3000, 40), yint_min = 0)

  #### RTSS in addition to PMC
  pplot1 <- ggplot(data = subset(pdat, pmc_coverage != 0 & rtss_coverage != 0)) +
    geom_hline(aes(yintercept = yint), alpha = 0) +
    geom_col(aes(x = pmc_mode_fct, y = median_val, col = pmc_mode_fct),
             fill = 'NA', position = position_dodge(width = 0.9), linetype = 'dashed',
             size = 1, width = 0.8, show.legend = F, alpha = 0.6) +
    geom_col(data = subset(pdat, ((pmc_coverage != 0 & rtss_coverage == 0))),
             aes(x = pmc_mode_fct, y = median_val, fill = pmc_mode_fct),
             position = position_dodge(width = 0.9), size = 0.3, width = 0.8) +
    geom_errorbar(aes(x = pmc_mode_fct, y = median_val, ymin = low_val, ymax = up_val), width = 0.01) +
    facet_wrap(~name, nrow = 2, scales = 'free_y') +
    scale_y_continuous(labels = comma, expand = c(0, 0)) +
    labs(subtitle = 'RTS,S in addition to PMC\n',
         x = 'Age (months)',
         y = 'Cumulative cases averted per\n1000 population in children U2',
         color = '', fill = '') +
    scale_color_manual(values = c(pmc_cols[c(1, 2, 5)], rtss_col)) +
    scale_fill_manual(values = c(pmc_cols[c(1, 2, 5)], rtss_col)) +
    customTheme_nogrid +
    theme(legend.position = 'None')


  pplot2 <- ggplot(data = subset(pdat, pmc_coverage == 0.8 & rtss_coverage == 0.8)) +
    geom_hline(aes(yintercept = yint), alpha = 0) +
    geom_col(aes(x = pmc_mode_fct, y = median_val, col = pmc_mode_fct),
             position = position_dodge(width = 0.9), fill = 'NA', linetype = 'dashed',
             size = 1, width = 0.8, show.legend = F, alpha = 0.6) +
    geom_col(aes(x = pmc_mode_fct, y = median_val, fill = pmc_mode_fct),
             position = position_dodge(width = 0.9), size = 1, width = 0.8, show.legend = F) +
    geom_col(data = pdat_rtss,
             aes(x = pmc_mode_fct, y = median_val), fill = 'white', position = position_dodge(width = 0.9),
             size = 0.3, width = 0.8) +
    geom_errorbar(aes(x = pmc_mode_fct, y = median_val, ymin = low_val, ymax = up_val), width = 0.01) +
    facet_wrap(~name, nrow = 2, scales = 'free_y') +
    scale_y_continuous(labels = comma, expand = c(0, 0)) +
    labs(subtitle = 'PMC in addition to RTS,S\n',
         x = 'Age (months)',
         y = 'Cumulative cases averted per\n1000 population in children U2',
         color = '', fill = '') +
    scale_color_manual(values = c(pmc_cols[c(1, 2, 5)], rtss_col)) +
    scale_fill_manual(values = c(pmc_cols[c(1, 2, 5)], rtss_col)) +
    customTheme_nogrid +
    theme(legend.position = 'None')


  pplot <- plot_combine(list(pplot1, pplot2), ncol = 2)
  print(pplot)

  f_save_plot(pplot, paste0('Fig4C_cum'),
              file.path(plot_dir), width = 10, height = 5, units = 'in', device_format = device_format)
  fwrite(pdat, file.path(plot_dir, 'csv', 'Fig4C_dat_cum.csv'))


} #fig4C_cum

if (fig4_alternative) {

  scenario_labels <- c('None', 'PMC-3', 'PMC-4', 'PMC-5', 'PMC-6', 'PMC-7', 'RTS,S', 'PMC-3 + RTS,S')
  agegrp <- c('U2')
  outcome_cols <- c('clinical_cases_averted', 'severe_cases_averted')
  pmc_modes <- c("3tp", "4tp2ndyr", "5tp2ndyr", "6tp2ndyr", "7tp2ndyr")  # "4tp"  "5tp"
  pmc_labels <- c('PMC-3', 'PMC-4', 'PMC-5', 'PMC-6', 'PMC-7')

  model1 <- T
  if (model1) {
    exp_name <- 'generic_PMCmode_RTSS_vaccSP_IIV'
    exp_name_rtss_cov_scen <- 'generic_single_RTSS_vaccSP_IIV'
    exp_name_pmc3rtss <- 'generic_PMCmode_RTSS_cov_vaccSP_IIV'

    cases_df_pmc <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
      filter(age_group %in% agegrp) %>%
      filter(pmc_mode %in% pmc_modes,
             pmc_coverage > 0 & rtss_coverage == 0) %>%
      mutate(coverage = pmc_coverage) %>%
      dplyr::select(-pmc_coverage, -pmc_rtss_cov)

    cases_df <- fread(file.path(simout_dir, exp_name_rtss_cov_scen, 'simdat_aggr_agegroup.csv')) %>%
      filter(age_group %in% agegrp,
             pmc_coverage == 0 & rtss_coverage > 0) %>%
      mutate(coverage = rtss_coverage,
             pmc_mode = 'rtss') %>%
      dplyr::select(-pmc_coverage, -pmc_rtss_cov) %>%
      bind_rows(cases_df_pmc)

    cases_df_0 <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
      filter(age_group %in% agegrp) %>%
      filter(pmc_mode %in% c('3tp'),
             pmc_coverage == 0 & rtss_coverage == 0) %>%
      mutate(coverage = 0) %>%
      dplyr::select(-pmc_coverage, -pmc_rtss_cov, -pmc_mode) %>%
      left_join(cases_df[, c('name', 'Annual_EIR', 'age_group', 'pmc_mode')]) ## repeat by pmc_mode

    cases_df_pmctss <- fread(file.path(simout_dir, exp_name_pmc3rtss, 'simdat_aggr_agegroup.csv')) %>%
      filter(age_group %in% agegrp) %>%
      filter(pmc_mode %in% c('3tp'),
             pmc_coverage == rtss_coverage) %>%
      mutate(coverage = rtss_coverage,
             pmc_mode = 'rtss_pmc3') %>%
      dplyr::select(-pmc_coverage, -pmc_rtss_cov)

    cases_df <- cases_df %>%
      bind_rows(cases_df_0) %>%
      bind_rows(cases_df_pmctss)

    cases_df$scen <- factor(cases_df$pmc_mode,
                            levels = c(pmc_modes, 'rtss', 'rtss_pmc3'),
                            labels = c(pmc_labels, 'RTS,S', 'PMC-3 + RTS,S'))
  } #model1


  plot_df_averted <- cases_df %>%
    filter(name %in% c('clinical_cases_averted', 'severe_cases_averted')) %>%
    mutate(name = ifelse(name == 'clinical_cases_averted', 'clinical', 'severe'))

  plot_df_PE <- cases_df %>%
    filter(name %in% c('PE_clinical_incidence', 'PE_severe_incidence')) %>%
    rename_with(~sub('_val', '_valPE', .x)) %>%
    mutate(name = ifelse(name == 'PE_clinical_incidence', 'clinical', 'severe'))

  plot_df <- plot_df_averted %>%
    left_join(plot_df_PE) %>%
    mutate(yconversion = up_val / up_valPE,
           yconversion = ifelse(is.na(yconversion), 1, yconversion))

  plot_df %>%
    filter(name == 'clinical') %>%
    summarize(yconversion = max(yconversion))
  plot_df %>%
    filter(name == 'severe') %>%
    summarize(yconversion = max(yconversion))
  #
  ggplot(data = subset(plot_df, name == 'clinical')) +
    geom_line(aes(x = coverage, y = mean_valPE, col = scen), show.legend = F) +
    geom_ribbon(aes(x = coverage, ymin = low_valPE, ymax = up_valPE, fill = scen), show.legend = F, alpha = 0.3) +
    geom_point(aes(x = coverage, y = mean_valPE, col = scen, fill = scen, shape = 'generic')) +
    scale_shape_manual(values = c(21, 23, 25)) +
    labs(shape = 'Model type',
         y = 'Cases averted',
         col = 'Scenario', fill = 'Scenario') +
    facet_wrap(~name, scales = 'free') +
    f_getCustomTheme()

  p1 <- ggplot(data = subset(plot_df, name == 'clinical')) +
    geom_line(aes(x = coverage, y = mean_val, col = scen), show.legend = F) +
    geom_ribbon(aes(x = coverage, ymin = low_val, ymax = up_val, fill = scen), show.legend = F, alpha = 0.3) +
    geom_point(aes(x = coverage, y = mean_val, col = scen, fill = scen), shape=21) +
    scale_y_continuous(lim=c(NA, 1500),"Annual cases averted\nper 1000 population", sec.axis = sec_axis(~. / 3100*100, name = "PE (% reduction)")) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25) * 100) +
    scale_color_manual(values = custom_cols) +
    scale_fill_manual(values = custom_cols) +
    labs(y = 'Cases averted', x = 'Coverage (%)',
         col = 'Scenario', fill = 'Scenario') +
    facet_wrap(~name, scales = 'free') +
    f_getCustomTheme()


  p2 <- ggplot(data = subset(plot_df, name == 'severe')) +
    geom_line(aes(x = coverage, y = mean_val, col = scen), show.legend = F) +
    geom_ribbon(aes(x = coverage, ymin = low_val, ymax = up_val, fill = scen), show.legend = F, alpha = 0.3) +
    geom_point(aes(x = coverage, y = mean_val, col = scen, fill = scen), shape=21) +
    scale_y_continuous(lim=c(NA, 25), "Annual cases averted\nper 1000 population", sec.axis = sec_axis(~. / 30*100, name = "PE (% reduction)")) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25) * 100) +
    scale_color_manual(values = custom_cols) +
    scale_fill_manual(values = custom_cols) +
    labs(y = 'Cases averted', x = 'Coverage (%)',
         col = 'Scenario', fill = 'Scenario') +
    facet_wrap(~name, scales = 'free') +
    f_getCustomTheme()


  pplot <- plot_combine(list(p1, p2), labels = c('a', 'b'), ncol = 1)
  pplot

  f_save_plot(pplot, paste0('Fig_04C_alternative'), file.path(plot_dir),
              width = 6, height = 8, units = 'in', device_format = device_format)


}


### Tables
if (result_tables) {
  cases_df_agegrp <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    filter(age_group %in% c('U1', 'U2'))

  cases_df <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_year.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    mutate(year = ifelse(year >= 2020, year - 2020, year)) %>%
    filter(year < 2) %>%
    mutate(age_group = paste0('>', year)) %>%
    dplyr::select(colnames(cases_df_agegrp)) %>%
    bind_rows(cases_df_agegrp)

  pdat <- cases_df %>%
    mutate(pmc_mode = ifelse(pmc_rtss_cov == '0-0.8', 'RTS,S', pmc_mode),
           pmc_mode = ifelse(pmc_rtss_cov == '0-0', 'None', pmc_mode)) %>%
    gen_pmc_mode_fct() %>%
    filter((pmc_rtss_cov == '0.8-0.8' & pmc_mode == '3tp') | (pmc_rtss_cov %in% c('0-0.8', '0.8-0'))) %>%
    filter(!(pmc_mode %in% c('4tp', '5tp'))) %>%
    mutate(yint = ifelse(name == 'clinical_cases', 110, 1.5), yint_min = 0,
           pmc_mode_fct = ifelse(pmc_rtss_cov == '0.8-0.8', 'PMC-3+RTS,S', as.character(pmc_mode_fct)))

  pdat$pmc_mode_fct <- factor(pdat$pmc_mode_fct, levels = scenario_labels, labels = scenario_labels)
  #dat_table(pdat, qc(age_group, Annual_EIR, name))

  pdat %>%
    filter(name %in% c('clinical_cases', 'clinical_cases_averted', 'PE_clinical_incidence') & age_group != 'U1') %>%
    mutate(median = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(median_val * 100, 1), format_num(median_val)),
           low = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(low_val * 100, 1), format_num(low_val)),
           up = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(up_val * 100, 1), format_num(up_val)),
           value = paste0(median, ' (', low, '-', up, ')')) %>%
    select(age_group, pmc_mode_fct, name, value) %>%
    pivot_wider(names_from = name, values_from = value)

  pdat %>%
    filter(name %in% c('severe_cases', 'severe_cases_averted', 'PE_severe_incidence') & age_group != 'U1') %>%
    mutate(median = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(median_val * 100, 1), format_num(median_val)),
           low = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(low_val * 100, 1), format_num(low_val)),
           up = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(up_val * 100, 1), format_num(up_val)),
           value = paste0(median, ' (', low, '-', up, ')')) %>%
    select(age_group, pmc_mode_fct, name, value) %>%
    pivot_wider(names_from = name, values_from = value)

}