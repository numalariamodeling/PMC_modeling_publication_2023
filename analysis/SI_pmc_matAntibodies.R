#setwd(file.path(getwd(),'manuscript'))
model <- 'generic'
source(file.path('00_config.R'))


exp_name <- 'generic_PMC_matAb_vaccSP_IIV'
dat_aggr <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv'))
dat_aggr_mth <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_month.csv'))

plotdat <- dat_aggr %>%
  filter(age_group %in% c('U1') &
           ipti_coverage != 0 &
           name %in% c('PE_clinical_incidence', 'PE_severe_incidence')) %>%
  mutate(name = ifelse(name == 'PE_clinical_incidence', 'Clinical malaria U1', 'Severe malaria U1'))

plotdat_aggr <- dat_aggr %>%
  filter(age_group %in% c('U1') &
           ipti_coverage != 0 &
           name %in% c('PE_clinical_incidence', 'PE_severe_incidence')) %>%
  mutate(name = ifelse(name == 'PE_clinical_incidence', 'Clinical malaria U1', 'Severe malaria U1')) %>%
  ungroup() %>%
  group_by(name, age_group, Annual_EIR) %>%
  summarize(min_val = min(mean_val),
            max_val = max(mean_val))


pplot <- ggplot(data = plotdat) +
  geom_hline(yintercept = c(0, 0.5), alpha = 0) +
  geom_errorbar(data = plotdat_aggr, aes(x = Annual_EIR, ymin = min_val, ymax = max_val), width = 0, size = 0.8) +
  geom_point(data = subset(plotdat, matAbs == 0.1327), aes(x = Annual_EIR, y = mean_val), col = 'deepskyblue3', size = 1.6) +
  facet_wrap(~name, nrow = 1) +
  labs(x = 'annual EIR (log scale)') +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.5, 0.05), labels = seq(0, 0.5, 0.05)) +
  scale_x_log10(breaks = unique(plotdat$Annual_EIR), labels = unique(plotdat$Annual_EIR)) +
  customTheme
pplot
f_save_plot(pplot, paste0('SI_fig_PMC_maternalantibodies'),
            file.path(plot_dir), width =8, height =6, units = 'in', device_format = device_format)




pplot <- dat_aggr_mth %>%
  filter(age <= 365 &
           ipti_coverage == 0 &
           name %in% c('PfHRP2_Prevalence', 'clinical_cases', 'severe_cases')) %>%
  mutate(EIR = Annual_EIR,
         name = case_when(name == 'PfHRP2_Prevalence' ~ ' PfPR',
                          name == 'clinical_cases' ~ 'Clinical malaria',
                          name == 'severe_cases' ~ 'Severe malaria')) %>%
  ggplot() +
  geom_hline(yintercept = c(0, 0.5), alpha = 0) +
  geom_ribbon(aes(x = age/(365/12), ymin = low_val, ymax = up_val, fill = as.factor(matAbs)), alpha = 0.4) +
  geom_line(aes(x = age/(365/12), y = mean_val, col = as.factor(matAbs))) +
  #geom_smooth(aes(x = age, y = mean_val , col=as.factor(matAbs)), se=FALSE, span=0.5) +
  facet_wrap(name ~ EIR, nrow = 3, scales = 'free', labeller = labeller(EIR = label_both, name = label_value)) +
    scale_color_viridis_d(option = 'B', end = 0.8, direction = 1) +
  scale_fill_viridis_d(option = 'B', end = 0.8, direction = 1) +
  scale_x_continuous(breaks=c(0:12), labels=c(0:12)) +
  labs(x = 'Age (months)', y = 'Projected mean value',
       fill='maternal antibody\nprotection\nscaling',
       col='maternal antibody\nprotection\nscaling') +
  customTheme
pplot
f_save_plot(pplot, paste0('SI_fig_cases_maternalantibodies'),
            file.path(plot_dir), width = 8, height = 6, units = 'in', device_format = device_format)



