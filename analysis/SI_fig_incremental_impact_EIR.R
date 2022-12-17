source(file.path('analysis', '_config.R'))

exp_name <- 'generic_PMCmode_EIR_vaccSP_IIV'

dat_aggr <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
  rename_with(~gsub('ipti', 'pmc', .x)) %>%
  filter(age_group %in% c('U2') &
           pmc_mode != '4tp' &
           pmc_mode != '5tp')

table(dat_aggr$name)


dat_aggr$pmc_mode <- factor(dat_aggr$pmc_mode,
                            levels = c("3tp", "4tp2ndyr", "5tp2ndyr", "6tp2ndyr", "7tp2ndyr"),
                            labels = c("3 (2.5, 3.5,    9            )",
                                       "4 (2.5, 3.5,    9, 12        )",
                                       "5 (2.5, 3.5,    9, 12, 15    )",
                                       "6 (2.5, 3.5, 6, 9, 12, 15    )",
                                       "7 (2.5, 3.5, 6, 9, 12, 15, 18)"))

dat1 <- dat_aggr %>%
  filter(name == 'clinical_cases_averted' & pmc_coverage == 0.8) %>%
  dplyr::select(Annual_EIR, pmc_mode, median_val) %>%
  group_by(Annual_EIR) %>%
  mutate(diff = median_val - lag(median_val)) %>%
  filter(pmc_mode != '3 (2.5, 3.5,    9            )')

dat2 <- dat_aggr %>%
  filter(name == 'PE_clinical_incidence' & pmc_coverage == 0.8) %>%
  dplyr::select(Annual_EIR, pmc_mode, median_val) %>%
  group_by(Annual_EIR) %>%
  mutate(diff = (median_val - lag(median_val)) * 100) %>%
  filter(pmc_mode != '3 (2.5, 3.5,    9            )')

dat1s <- dat_aggr %>%
  filter(name == 'severe_cases_averted' & pmc_coverage == 0.8) %>%
  dplyr::select(Annual_EIR, pmc_mode, median_val) %>%
  group_by(Annual_EIR) %>%
  mutate(diff = median_val - lag(median_val)) %>%
  filter(pmc_mode != '3 (2.5, 3.5,    9            )')

dat2s <- dat_aggr %>%
  filter(name == 'PE_severe_incidence' & pmc_coverage == 0.8) %>%
  dplyr::select(Annual_EIR, pmc_mode, median_val) %>%
  group_by(Annual_EIR) %>%
  mutate(diff = (median_val - lag(median_val)) * 100) %>%
  filter(pmc_mode != '3 (2.5, 3.5,    9            )')

dat1 %>%
  dplyr::select(-median_val) %>%
  pivot_wider(names_from = pmc_mode, values_from = c(diff))

dat2 %>%
  dplyr::select(-median_val) %>%
  pivot_wider(names_from = pmc_mode, values_from = c(diff))


dat1s %>%
  dplyr::select(-median_val) %>%
  pivot_wider(names_from = pmc_mode, values_from = c(diff))

dat2s %>%
  dplyr::select(-median_val) %>%
  pivot_wider(names_from = pmc_mode, values_from = c(diff))


p1 <- ggplot(data = dat1) +
  geom_hline(yintercept = 0) +
  geom_col(aes(x = as.factor(Annual_EIR), y = diff, fill = pmc_mode,
               group = interaction(pmc_mode, Annual_EIR)), col = 'black',
           position = position_dodge(width = 0.75), width = 0.7) +
  #facet_wrap(~Annual_EIR, scales = 'free_y', nrow = 1, labeller = labeller(Annual_EIR = label_both)) +
  scale_fill_manual(values = pmc_cols[2:6]) +
  labs(x = 'Annual EIR', y = 'Difference to previous dose\n(cases averted per 1000 population)',
       fill = 'PMC doses\n(age in months)') +
  customTheme_nogrid


p2 <- ggplot(data = dat2) +
  geom_hline(yintercept = 0) +
  geom_col(aes(x = as.factor(Annual_EIR), y = diff, fill = pmc_mode,
               group = interaction(pmc_mode, Annual_EIR)), col = 'black',
           position = position_dodge(width = 0.75), width = 0.7) +
  #facet_wrap(~Annual_EIR, scales = 'free_y', nrow = 1, labeller = labeller(Annual_EIR = label_both)) +
  scale_fill_manual(values = pmc_cols[2:6]) +
  labs(x = 'Annual EIR', y = 'Relative difference to previous dose\n(in percentage points)', fill = 'PMC doses\n(age in months)') +
  customTheme_nogrid

pplot <- plot_combine(list(p1, p2), ncol = 1, labels = c('A', 'B'))
pplot

f_save_plot(pplot, paste0('SI_fig_incremental'), file.path(plot_dir), width = 12, height = 6, units = 'in', device_format = device_format)
#fwrite(plotdat, file.path(plot_dir, 'csv', 'SI_fig_incremental_dat.csv'))

alternative = T
if (alternative) {

  dat_aggr <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    filter(age_group %in% c('U2') &
             pmc_mode != '6tp2ndyr' &
             pmc_mode != '7tp2ndyr')

  table(dat_aggr$name)

  dat_aggr$infant_ext <- factor(dat_aggr$pmc_mode,
                                levels = c("3tp", "4tp", "4tp2ndyr", "5tp", "5tp2ndyr"),
                                labels = c("<=12 months",
                                           "6 months dose",
                                           "12 months or later",
                                           "6 months dose",
                                           "12 months or later"))

  dat_aggr$pmc_mode <- factor(dat_aggr$pmc_mode,
                              levels = c("3tp", "4tp", "4tp2ndyr", "5tp", "5tp2ndyr"),
                              labels = c("3 (2.5, 3.5,    9            )",
                                         "4 (+ 6)",
                                         "4 (+ 12)",
                                         "5 (+6, 12)",
                                         "5 (+12, 15)"))


  dat1 <- dat_aggr %>%
    filter(name == 'clinical_cases_averted' & pmc_coverage == 0.8) %>%
    dplyr::select(Annual_EIR, pmc_mode, infant_ext, median_val) %>%
    group_by(Annual_EIR) %>%
    mutate(diff = median_val - median_val[pmc_mode == '3 (2.5, 3.5,    9            )']) %>%
    filter(pmc_mode != '3 (2.5, 3.5,    9            )')

  dat2 <- dat_aggr %>%
    filter(name == 'PE_clinical_incidence' & pmc_coverage == 0.8) %>%
    dplyr::select(Annual_EIR, pmc_mode, infant_ext, median_val) %>%
    group_by(Annual_EIR) %>%
    mutate(diff = median_val - median_val[pmc_mode == '3 (2.5, 3.5,    9            )']) %>%
    filter(pmc_mode != '3 (2.5, 3.5,    9            )')

  dat1s <- dat_aggr %>%
    filter(name == 'severe_cases_averted' & pmc_coverage == 0.8) %>%
    dplyr::select(Annual_EIR, pmc_mode, infant_ext, median_val) %>%
    group_by(Annual_EIR) %>%
    mutate(diff = median_val - median_val[pmc_mode == '3 (2.5, 3.5,    9            )']) %>%
    filter(pmc_mode != '3 (2.5, 3.5,    9            )')

  dat2s <- dat_aggr %>%
    filter(name == 'PE_severe_incidence' & pmc_coverage == 0.8) %>%
    dplyr::select(Annual_EIR, pmc_mode, infant_ext, median_val) %>%
    group_by(Annual_EIR) %>%
    mutate(diff = median_val - median_val[pmc_mode == '3 (2.5, 3.5,    9            )']) %>%
    filter(pmc_mode != '3 (2.5, 3.5,    9            )')

  dat1 %>%
    dplyr::select(-median_val) %>%
    pivot_wider(names_from = pmc_mode, values_from = c(diff))

  dat2 %>%
    dplyr::select(-median_val) %>%
    pivot_wider(names_from = pmc_mode, values_from = c(diff))


  dat1s %>%
    dplyr::select(-median_val) %>%
    pivot_wider(names_from = pmc_mode, values_from = c(diff))

  dat2s %>%
    dplyr::select(-median_val) %>%
    pivot_wider(names_from = pmc_mode, values_from = c(diff))


  p1 <- ggplot(data = dat1) +
    geom_hline(yintercept = 0) +
    geom_col(aes(x = as.factor(Annual_EIR), y = diff, alpha = infant_ext, fill = pmc_mode,
                 group = interaction(pmc_mode, Annual_EIR)), col = 'black',
             position = position_dodge(width = 0.75), width = 0.7) +
    #facet_wrap(~Annual_EIR, scales = 'free_y', nrow = 1, labeller = labeller(Annual_EIR = label_both)) +
    scale_fill_manual(values = c('#fc8d59', '#fc8d59', '#99d594', '#99d594')) +
    scale_alpha_manual(values = c(0.5, 1)) +
    labs(title = 'Clinical malaria U2\n',
         x = 'Annual EIR', alpha = '', y = 'Difference to previous dose\n(cases averted per 1000 population)',
         fill = 'PMC doses\n(age in months)') +
    customTheme_nogrid


  p2 <- ggplot(data = dat2) +
    geom_hline(yintercept = 0) +
    geom_col(aes(x = as.factor(Annual_EIR), y = diff, alpha = infant_ext, fill = pmc_mode,
                 group = interaction(pmc_mode, Annual_EIR)), col = 'black',
             position = position_dodge(width = 0.75), width = 0.7) +
    #facet_wrap(~Annual_EIR, scales = 'free_y', nrow = 1, labeller = labeller(Annual_EIR = label_both)) +
    scale_fill_manual(values = c('#fc8d59', '#fc8d59', '#99d594', '#99d594')) +
    scale_alpha_manual(values = c(0.5, 1)) +
    labs(title = 'Clinical malaria U2\n',
         x = 'Annual EIR', alpha = '', y = 'Relative difference to previous dose\n(in percentage points)', fill = 'PMC doses\n(age in months)') +
    customTheme_nogrid

  pplot <- plot_combine(list(p1, p2), ncol = 1, labels = c('A', 'B'))
  pplot

  p1s <- ggplot(data = dat1s) +
    geom_hline(yintercept = 0) +
    geom_col(aes(x = as.factor(Annual_EIR), y = diff, alpha = infant_ext, fill = pmc_mode,
                 group = interaction(pmc_mode, Annual_EIR)), col = 'black',
             position = position_dodge(width = 0.75), width = 0.7) +
    #facet_wrap(~Annual_EIR, scales = 'free_y', nrow = 1, labeller = labeller(Annual_EIR = label_both)) +
    scale_fill_manual(values = c('#fc8d59', '#fc8d59', '#99d594', '#99d594')) +
    scale_alpha_manual(values = c(0.5, 1)) +
    labs(title = 'Severe malaria U2\n', x = 'Annual EIR', alpha = '', y = 'Difference to previous dose\n(cases averted per 1000 population)',
         fill = 'PMC doses\n(age in months)') +
    customTheme_nogrid


  p2s <- ggplot(data = dat2s) +
    geom_hline(yintercept = 0) +
    geom_col(aes(x = as.factor(Annual_EIR), y = diff, alpha = infant_ext, fill = pmc_mode,
                 group = interaction(pmc_mode, Annual_EIR)), col = 'black',
             position = position_dodge(width = 0.75), width = 0.7) +
    #facet_wrap(~Annual_EIR, scales = 'free_y', nrow = 1, labeller = labeller(Annual_EIR = label_both)) +
    scale_fill_manual(values = c('#fc8d59', '#fc8d59', '#99d594', '#99d594')) +
    scale_alpha_manual(values = c(0.5, 1)) +
    labs(title = 'Severe malaria U2\n', x = 'Annual EIR', alpha = '',
         y = 'Relative difference to previous dose\n(in percentage points)',
         fill = 'PMC doses\n(age in months)') +
    customTheme_nogrid

  pplot2 <- plot_combine(list(p1, p1s, p2, p2s), ncol = 2, labels = c('A', 'B', 'C', 'D'))
  pplot2


  f_save_plot(pplot2, paste0('SI_fig_incremental_v2'), file.path(plot_dir), width = 12, height = 8, units = 'in', device_format = device_format)
  #fwrite(plotdat, file.path(plot_dir, 'csv', 'SI_fig_incremental_dat.csv'))


}