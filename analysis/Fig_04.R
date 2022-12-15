source(file.path('analysis', '_config.R'))
exp_name <- 'generic_PMCmode_RTSS_vaccSP_IIV'

fig4A <- TRUE
fig4B <- TRUE
fig4B_cum <- TRUE
fig4C <- TRUE
fig4C_cum <- TRUE
y <- 2.0
scenario_labels <- c('None', 'PMC-3', 'PMC-4', 'PMC-5', 'PMC-6', 'PMC-7', 'RTS,S', 'PMC-3 + RTS,S')

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
        geom_hline(yintercept = c( 0.45), alpha=0) +
    geom_hline(yintercept = c(0), size = 0.2) +
    geom_col(aes(x = pmc_mode_fct, y = median_val, alpha = age_group, group = interaction(age_group, pmc_mode_fct), fill = pmc_mode_fct),
             position = position_dodge(width = 0.75), size = 0.3, width = 0.7, col = 'black') +
    geom_errorbar(aes(x = pmc_mode_fct, y = median_val, ymin = low_val, ymax = up_val, group = interaction(age_group, pmc_mode_fct)),
                  position = position_dodge(width = 0.75), width = 0.01) +
    facet_wrap(~name, nrow = 2, scales = 'free_y') +
    scale_y_continuous( expand = c(0, 0), breaks=seq(-0.1, 0.5, 0.1), labels=seq(-0.1, 0.5, 0.1)*100) +
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
      pmc_mode_fct2 != 'None' &
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
             aes(x = pmc_mode_fct, y = median_val, fill = pmc_mode_fct), position = position_dodge(), size = 0.3, width = w - 0.1) +
    geom_col(aes(x = pmc_mode_fct, y = median_val, linetype = age_group, group = interaction(age_group, pmc_mode_fct)),
             fill = 'white', col = 'black', position = position_dodge(), alpha = 0.7, size = 0.3, width = w) +
    geom_errorbar(aes(x = pmc_mode_fct, y = median_val, ymin = low_val, ymax = up_val, group = interaction(age_group, pmc_mode_fct)),
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
    geom_col(aes(x = pmc_mode_fct, y = median_val, col = pmc_mode_fct),
             position = position_dodge(width = 0.9), fill = 'NA', linetype = 'dashed', size = 1, width = 0.8, show.legend = F, alpha = 0.6) +
    geom_col(aes(x = pmc_mode_fct, y = median_val_subrtss, fill = pmc_mode_fct),
             position = position_dodge(width = 0.9), size = 1, width = 0.8, show.legend = F) +
    #geom_col(data = pdat_rtss,
    #         aes(x = pmc_mode_fct, y = median_val), fill = 'white', position = position_dodge(width = 0.9), size = 0.3, width = 0.8) +
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


  #### PMC in addition to RTSS
  # pplot2 <- ggplot(data = subset(pdat, pmc_coverage == 0.8 & rtss_coverage == 0.8)) +
  #   geom_hline(aes(yintercept = yint), alpha = 0) +
  #   geom_col(aes(x = pmc_mode_fct, y = median_val, col = pmc_mode_fct),
  #            fill = 'NA', position = position_dodge(width = 0.9), linetype = 'dashed', size = 1, width = 0.8, show.legend = F, alpha = 0.6) +
  #   geom_col(data = pdat_rtss,
  #            aes(x = pmc_mode_fct, y = median_val, fill = pmc_mode_fct), position = position_dodge(width = 0.9), size = 0.3, width = 0.8) +
  #   geom_errorbar(aes(x = pmc_mode_fct, y = median_val, ymin = low_val, ymax = up_val), width = 0.01) +
  #   facet_wrap(~name, nrow = 2, scales = 'free_y') +
  #   scale_y_continuous(labels = comma, expand = c(0, 0)) +
  #   labs(subtitle = 'PMC in addition to RTS,S\n',
  #        x = 'Age (months)',
  #        y = 'Cumulative cases averted per\n1000 population in children U2',
  #        color = '', fill = '') +
  #   scale_color_manual(values = c(pmc_cols[c(1, 2, 5)], rtss_col)) +
  #   scale_fill_manual(values = c(pmc_cols[c(1, 2, 5)], rtss_col)) +
  #   customTheme_nogrid +
  #   theme(legend.position = 'None')


  pplot2 <- ggplot(data = subset(pdat, pmc_coverage == 0.8 & rtss_coverage == 0.8)) +
    geom_hline(aes(yintercept = yint), alpha = 0) +
    geom_col(aes(x = pmc_mode_fct, y = median_val, col = pmc_mode_fct),
             position = position_dodge(width = 0.9), fill = 'NA', linetype = 'dashed', size = 1, width = 0.8, show.legend = F, alpha = 0.6) +
    geom_col(aes(x = pmc_mode_fct, y = median_val, fill = pmc_mode_fct),
             position = position_dodge(width = 0.9), size = 1, width = 0.8, show.legend = F) +
    geom_col(data = pdat_rtss,
             aes(x = pmc_mode_fct, y = median_val), fill = 'white', position = position_dodge(width = 0.9), size = 0.3, width = 0.8) +
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

### Tables
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


dat_table(pdat, qc(age_group, Annual_EIR, name))

pdat %>%
  filter(name %in% c('clinical_cases', 'clinical_cases_averted', 'PE_clinical_incidence') & age_group != 'U1') %>%
  mutate(median = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(median_val * 100, 1), format_num(median_val)),
         low = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(low_val * 100, 1), format_num(low_val)),
         up = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(up_val * 100, 1), format_num(up_val)),
         value = paste0(median, ' (', low, '-', up, ')')) %>%
  select(age_group, pmc_mode_fct2, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  kbl() %>%
  kable_styling(bootstrap_options = tbl_opts, full_width = F, position = "left", fixed_thead = T)

pdat %>%
  filter(name %in% c('severe_cases', 'severe_cases_averted', 'PE_severe_incidence') & age_group != 'U1') %>%
  mutate(median = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(median_val * 100, 1), format_num(median_val)),
         low = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(low_val * 100, 1), format_num(low_val)),
         up = ifelse(name %in% c('PE_clinical_incidence', 'PE_severe_incidence'), format_num(up_val * 100, 1), format_num(up_val)),
         value = paste0(median, ' (', low, '-', up, ')')) %>%
  select(age_group, pmc_mode_fct2, name, value) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  kbl() %>%
  kable_styling(bootstrap_options = tbl_opts, full_width = F, position = "left", fixed_thead = T)
