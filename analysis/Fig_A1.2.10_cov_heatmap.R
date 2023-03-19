##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## Fig_A1.2.10_cov_heatmap.R
##---------------------
source(file.path('analysis', '_config.R'))


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
  dplyr::select(-pmc_mode, -pmc_mode_fct) %>%
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

f_save_plot(pplot, paste0('Fig A1.2.10'),
            file.path(plot_dir), width = 8, height = 6, units = 'in', device_format = device_format)

#### tables
pdat %>%
  group_by(pmc_coverage, rtss_coverage, name, age_group) %>%
  summarize(median_val_low = min(median_val),
            median_val_mean = mean(median_val),
            median_val_up = max(median_val))




