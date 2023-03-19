##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## Fig_A1.3.3_clinicalTreatment_agevariation.R
##---------------------

source(file.path('analysis', '_config.R'))

exp_name <- 'generic_PMC_CM_vaccSP_IIV'
dat_aggr <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv'))
table(dat_aggr$pmc_rtss_cov, dat_aggr$cm_coverage)
summary(dat_aggr)

dat_aggr <- dat_aggr %>%
  filter(Annual_EIR < eir_max) %>%
  mutate(cm_coverage = as.factor(cm_coverage * 100))

##--------------------------------
pA1 <- dat_aggr %>%
  filter(age_group == 'U1' &
           name == 'clinical_cases_averted' &
           Annual_EIR == 32) %>%
  mutate(name = 'Clinical cases U1') %>%
  ggplot(aes(x = pmc_coverage, y = median_val, ymin = low_val, ymax = up_val,
             fill = cm_coverage)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(alpha = 0.2) +
  geom_line(aes(col = cm_coverage)) +
  facet_wrap(~name, nrow = 2, scales = 'free') +
  scale_color_viridis_d(option = 'A', end = 0.7, direction = -1) +
  scale_fill_viridis_d(option = 'A', end = 0.7, direction = -1) +
  labs(col = '% cases treated', fill = '% cases treated',
       y = 'Cases averted\nper 1000 population',
       x = 'PMC coverage') +
  customTheme

pA2 <- dat_aggr %>%
  filter(age_group == 'U1' &
           name == 'severe_cases_averted' &
           Annual_EIR == 32) %>%
  mutate(name = 'Severe cases U1') %>%
  ggplot(aes(x = pmc_coverage, y = median_val, ymin = low_val, ymax = up_val,
             fill = cm_coverage)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(alpha = 0.2) +
  geom_line(aes(col = cm_coverage)) +
  facet_wrap(~name, nrow = 2, scales = 'free') +
  scale_color_viridis_d(option = 'A', end = 0.7, direction = -1) +
  scale_fill_viridis_d(option = 'A', end = 0.7, direction = -1) +
  labs(col = '% cases treated', fill = '% cases treated',
       y = 'Cases averted\nper 1000 population',
       x = 'PMC coverage') +
  customTheme


pB1 <- dat_aggr %>%
  filter(age_group == 'U1' &
           name == 'PE_clinical_incidence' &
           Annual_EIR == 32) %>%
  mutate(name = 'Clinical cases U1') %>%
  ggplot(aes(x = pmc_coverage, y = median_val * 100, ymin = low_val * 100, ymax = up_val * 100,
             fill = cm_coverage)) +
  geom_hline(yintercept = c(-10, 40), alpha = 0) +
  geom_hline(yintercept = 0) +
  geom_ribbon(alpha = 0.2) +
  geom_line(aes(col = cm_coverage)) +
  facet_wrap(~name, nrow = 2, scales = 'free') +
  scale_color_viridis_d(option = 'A', end = 0.7, direction = -1) +
  scale_fill_viridis_d(option = 'A', end = 0.7, direction = -1) +
  labs(col = '% cases treated', fill = '% cases treated', y = '% reduction in cases', x = 'PMC coverage') +
  customTheme

pB2 <- dat_aggr %>%
  filter(age_group == 'U1' &
           name == 'PE_severe_incidence' &
           Annual_EIR == 32) %>%
  mutate(name = 'Severe cases U1') %>%
  ggplot(aes(x = pmc_coverage, y = median_val * 100, ymin = low_val * 100, ymax = up_val * 100,
             fill = cm_coverage)) +
  geom_hline(yintercept = c(-10, 40), alpha = 0) +
  geom_hline(yintercept = 0) +
  geom_ribbon(alpha = 0.2) +
  geom_line(aes(col = cm_coverage)) +
  facet_wrap(~name, nrow = 2, scales = 'free') +
  scale_color_viridis_d(option = 'A', end = 0.7, direction = -1) +
  scale_fill_viridis_d(option = 'A', end = 0.7, direction = -1) +
  labs(col = '% cases treated', fill = '% cases treated', y = '% reduction in cases', x = 'PMC coverage') +
  customTheme

pplot <- plot_combine(list(pA1, pB1, pA2, pB2), ncol = 2, labels = c('A', 'B', '', ''))
pplot
f_save_plot(pplot, paste0('Fig A1.3.2'),
            file.path(plot_dir), width = 8, height = 5, units = 'in', device_format = device_format)


##--------------------------------
### Fig A1.3.3
##--------------------------------
exp_name <- 'generic_custom_CM_PMC_vaccSP_IIV'
dat_aggr <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv'))
table(dat_aggr$pmc_rtss_cov, dat_aggr$cm_coverage)
table(dat_aggr$name, dat_aggr$Annual_EIR)
summary(dat_aggr)

dat_aggr <- dat_aggr %>%
  filter(Annual_EIR < eir_max)

##----

pA1 <- dat_aggr %>%
  filter(pmc_coverage != 0 &
           age_group == 'U5' &
           name == 'clinical_cases_averted') %>%
  mutate(name = 'Clinical cases U5') %>%
  ggplot(aes(x = as.factor(Annual_EIR), y = median_val, ymin = low_val, ymax = up_val, fill = cm_coverage)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(alpha = 0.2, width = 0, position = position_dodge(width = 0.8)) +
  geom_col(aes(col = cm_coverage), width = 0.7, position = position_dodge(width = 0.8)) +
  facet_wrap(~name, nrow = 2, scales = 'free') +
  scale_color_viridis_d(option = 'A', end = 0.7, direction = -1) +
  scale_fill_viridis_d(option = 'A', end = 0.7, direction = -1) +
  labs(col = '% cases treated', fill = '% cases treated', y = 'Cases averted\nper 1000 population', x = 'Annual EIR') +
  customTheme

pA2 <- dat_aggr %>%
  filter(pmc_coverage != 0 &
           age_group == 'U5' &
           name == 'severe_cases_averted') %>%
  mutate(name = 'Severe cases U5') %>%
  ggplot(aes(x = as.factor(Annual_EIR), y = median_val, ymin = low_val, ymax = up_val, fill = cm_coverage)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(alpha = 0.2, width = 0, position = position_dodge(width = 0.8)) +
  geom_col(aes(col = cm_coverage), width = 0.7, position = position_dodge(width = 0.8)) +
  facet_wrap(~name, nrow = 2, scales = 'free') +
  scale_color_viridis_d(option = 'A', end = 0.7, direction = -1) +
  scale_fill_viridis_d(option = 'A', end = 0.7, direction = -1) +
  labs(col = '% cases treated', fill = '% cases treated',
       y = 'Cases averted\nper 1000 population',
       x = 'Annual EIR') +
  customTheme


pB1 <- dat_aggr %>%
  filter(pmc_coverage != 0 &
           age_group == 'U5' &
           name == 'PE_clinical_incidence') %>%
  mutate(name = 'Clinical cases U5') %>%
  ggplot(aes(x = as.factor(Annual_EIR), y = median_val * 100,
             ymin = low_val * 100,
             ymax = up_val * 100, fill = cm_coverage)) +
  geom_hline(yintercept = c(0, 10), alpha = 0) +
  geom_hline(yintercept = 0) +
  geom_errorbar(alpha = 0.2, width = 0, position = position_dodge(width = 0.8)) +
  geom_col(aes(col = cm_coverage), width = 0.7, position = position_dodge(width = 0.8)) +
  facet_wrap(~name, nrow = 2, scales = 'free') +
  scale_color_viridis_d(option = 'A', end = 0.7, direction = -1) +
  scale_fill_viridis_d(option = 'A', end = 0.7, direction = -1) +
  labs(col = '% cases treated', fill = '% cases treated', y = '% reduction in cases', x = 'Annual EIR') +
  customTheme

pB2 <- dat_aggr %>%
  filter(pmc_coverage != 0 &
           age_group == 'U5' &
           name == 'PE_severe_incidence') %>%
  mutate(name = 'Severe cases U5') %>%
  ggplot(aes(x = as.factor(Annual_EIR), y = median_val * 100,
             ymin = low_val * 100,
             ymax = up_val * 100, fill = cm_coverage)) +
  geom_hline(yintercept = c(-10, 10), alpha = 0) +
  geom_hline(yintercept = 0) +
  geom_errorbar(alpha = 0.2, width = 0, position = position_dodge(width = 0.8)) +
  geom_col(aes(col = cm_coverage), width = 0.7, position = position_dodge(width = 0.8)) +
  facet_wrap(~name, nrow = 2, scales = 'free') +
  scale_color_viridis_d(option = 'A', end = 0.7, direction = -1) +
  scale_fill_viridis_d(option = 'A', end = 0.7, direction = -1) +
  labs(col = '% cases treated', fill = '% cases treated', y = '% reduction in cases', x = 'Annual EIR') +
  customTheme

pplot <- plot_combine(list(pA1, pB1, pA2, pB2), ncol = 2, labels = c('A', 'B', '', ''))
pplot
f_save_plot(pplot, paste0('Fig A1.3.3'),
            file.path(plot_dir), width = 8, height = 5, units = 'in', device_format = device_format)

## calculate treatment coverage
dat_aggr_mth <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_month.csv'))

dat_aggr_mth_treatments <- dat_aggr_mth %>%
  filter(name %in% c('New_Clinical_Cases', 'New_Severe_Cases',
                     'Received_Treatment', 'Received_Severe_Treatment')) %>%
  select(age, year, Annual_EIR, pmc_coverage, cm_coverage, name, mean_val) %>%
  pivot_wider(names_from = name, values_from = mean_val) %>%
  group_by(age, year, Annual_EIR, pmc_coverage, cm_coverage) %>%
  mutate(frac_clinical_treated = Received_Treatment / New_Clinical_Cases,
         frac_severe_treated = Received_Severe_Treatment / New_Severe_Cases) %>%
  pivot_longer(cols = -c(age, year, Annual_EIR, pmc_coverage, cm_coverage))


dat_aggr_mth %>%
  filter(pmc_coverage == 0 &
           Annual_EIR == 32 &
           name %in% c('Received_Treatment', 'Received_Severe_Treatment')) %>%
  ggplot(aes(x = age, y = median_val, ymin = low_val, ymax = up_val, fill = cm_coverage)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(col = cm_coverage), width = 0.7, position = position_dodge(width = 0.8)) +
  facet_wrap(~name, nrow = 2, scales = 'free') +
  scale_color_viridis_d(option = 'A', end = 0.7, direction = -1) +
  scale_fill_viridis_d(option = 'A', end = 0.7, direction = -1) +
  labs(col = '% cases treated', fill = '% cases treated', y = '% reduction in cases', x = 'Age') +
  customTheme

pplot <- dat_aggr_mth_treatments %>%
  filter(pmc_coverage == 0 &
           Annual_EIR == 32 &
           name %in% c('frac_clinical_treated', 'frac_severe_treated')) %>%
  ggplot(aes(x = age / (365), y = value, fill = cm_coverage)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(col = cm_coverage)) +
  facet_wrap(~name, nrow = 2, scales = 'free') +
  scale_color_viridis_d(option = 'A', end = 0.7, direction = -1) +
  scale_fill_viridis_d(option = 'A', end = 0.7, direction = -1) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.1), labels = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = c(0:6), labels = c(0:6)) +
  labs(col = '% cases treated', fill = '% cases treated', y = '% reduction in cases', x = 'Age (years)') +
  customTheme
pplot

f_save_plot(pplot, paste0('Fig A1.3.3_smallfig'),
            file.path(plot_dir), width = 8, height = 5, units = 'in', device_format = device_format)
