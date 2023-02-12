##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## Fig_03.R
##---------------------
source(file.path('analysis', '_config.R'))

(exp_name <- 'generic_PMC3mode_redistribution_vaccSP_IIV')

dat_aggr <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv'))
dat_aggr_mth <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_month.csv'))

table(dat_aggr$pmc_mode, dat_aggr$pmc_coverage)
table(dat_aggr$age_group, dat_aggr$Annual_EIR)
dat_aggr$age_group <- factor(dat_aggr$age_group, levels = age_groups, labels = age_groups)

yr <- function(x) {
  return(as.character(round(x / (365 / 12), 1)))
}

yr(c(76, 106, 274))
yr(c(106, 274))
yr(c(106, 182, 274))
yr(c(106, 274, 365))
yr(c(106, 274, 456))


dat_aggr_mth$pmc_mode_fct <- factor(dat_aggr_mth$pmc_mode,
                                    levels = c('2tp', '3tp', '3atp', '3btp', '3ctp'),
                                    labels = c('        3.5,      9', '2.5, 3.5,      9', '        3.5, 6,  9',
                                               '        3.5,      9, 12', '        3.5,      9,     15'))

dat_aggr$pmc_mode_fct <- factor(dat_aggr$pmc_mode,
                                levels = c('2tp', '3tp', '3atp', '3btp', '3ctp'),
                                labels = c('        3.5,      9', '2.5, 3.5,      9', '        3.5, 6,  9',
                                           '        3.5,      9, 12', '        3.5,      9,     15'))

pplot <- dat_aggr %>%
  filter(age_group %in% c('U2') &
           pmc_coverage != 0 &
           name %in% c('clinical_cases_averted', 'severe_cases_averted')) %>%
  mutate(name = ifelse(name == 'clinical_cases_averted', 'Clinical malaria U2', 'Severe malaria U2')) %>%
  ggplot() +
  geom_col(aes(x = as.factor(Annual_EIR), y = mean_val, fill = pmc_mode_fct, col = pmc_mode_fct),
           width = 0.5, position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(x = as.factor(Annual_EIR), ymin = low_val, ymax = up_val, group = pmc_mode_fct),
                width = 0, position = position_dodge(width = 0.9)) +
  facet_wrap(~name, nrow = 2, scales = 'free') +
  geom_hline(yintercept = c(0)) +
  labs(x = 'annual EIR', y = 'Cases averted\nper 1000 population U2',
       fill = 'Age at PMC dose\n(months)', col = 'Age at PMC dose\n(months)') +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = c('white', pmc_cols)) +
  scale_color_manual(values = c('deepskyblue3', pmc_cols)) +
  customTheme

pplot_PE <- dat_aggr %>%
  filter(age_group %in% c('U2') &
           pmc_coverage != 0 &
           name %in% c('PE_clinical_incidence', 'PE_severe_incidence')) %>%
  mutate(name = ifelse(name == 'PE_clinical_incidence', 'Clinical malaria U2', 'Severe malaria U2')) %>%
  group_by(Annual_EIR) %>%
  mutate(low_val = ifelse(low_val <= -0.1, -0.1, low_val)) %>%
  ggplot() +
  geom_col(aes(x = as.factor(Annual_EIR), y = mean_val * 100, fill = pmc_mode_fct, col = pmc_mode_fct),
           width = 0.5, position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(x = as.factor(Annual_EIR), ymin = low_val * 100, ymax = up_val * 100, group = pmc_mode_fct),
                width = 0, position = position_dodge(width = 0.9)) +
  facet_wrap(~name, nrow = 2, scales = 'free') +
  geom_hline(yintercept = c(0)) +
  labs(x = 'annual EIR', y = '% reduction in cases',
       fill = 'Age at PMC dose\n(months)', col = 'Age at PMC dose\n(months)') +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = c('white', pmc_cols)) +
  scale_color_manual(values = c('deepskyblue3', pmc_cols)) +
  customTheme


f_save_plot(pplot, paste0('Fig03'),
            file.path(plot_dir), width = 12, height = 6, units = 'in', device_format = device_format)


f_save_plot(pplot_PE, paste0('Fig03_PE'),
            file.path(plot_dir), width = 12, height = 6, units = 'in', device_format = device_format)

