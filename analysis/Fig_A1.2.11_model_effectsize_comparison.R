##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## Fig_A1.2.11_model_effectsize_comparison
##---------------------

source(file.path('analysis', '_config.R'))

scenario_labels <- c('None', 'PMC-3', 'PMC-4', 'PMC-5', 'PMC-6', 'PMC-7', 'RTS,S', 'PMC-3 + RTS,S')
agegrp <- c('U1', 'U2')
outcome_cols <- c('clinical_cases_averted', 'severe_cases_averted')

model1 <- T
model2 <- T

if (model1) {
  exp_name <- 'generic_PMCmode_RTSS_vaccSP_IIV'
  exp_name_rtss_cov_scen = 'generic_single_RTSS_vaccSP_IIV'

  cases_df_pmc <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
    filter(age_group %in% agegrp) %>%
    filter(pmc_mode %in% c('3tp', '5tp2ndyr', '7tp2ndyr'),
           rtss_coverage == 0) %>%
    mutate(coverage = pmc_coverage) %>%
    dplyr::select(-pmc_coverage, -pmc_rtss_cov)

  cases_df <- fread(file.path(simout_dir, exp_name_rtss_cov_scen, 'simdat_aggr_agegroup.csv')) %>%
    filter(age_group %in% agegrp,
           pmc_coverage == 0) %>%
    mutate(coverage = rtss_coverage,
           pmc_mode = 'rtss') %>%
    dplyr::select(-pmc_coverage, -pmc_rtss_cov) %>%
    bind_rows(cases_df_pmc)

  cases_df$scen <- factor(cases_df$pmc_mode,
                          levels = c('3tp', '5tp2ndyr', '7tp2ndyr', 'rtss'),
                          labels = c('PMC-3', 'PMC-5', 'PMC-7', 'RTS,S'))
} #model1

###
if (model2) {
  cases_df_nga <- fread(file.path(simout_dir, 'NGA_simdat_aggr_agegroup.csv')) %>%
    filter(age_group %in% agegrp &
             scen != 'counterfactual' &
             scen != 'PMC-3 + RTS,S' &
             statistic == 'mean_val') %>%
    rename(State = seasonality)


  ### add actual coverage values
  cov_df <- fread(file.path(data_path, 'dhs_epi_dat_extended.csv')) %>%
    group_by(State) %>%
    summarize(vacc_cov_adj = mean(vacc_cov_adj))

  cases_df_nga <- cases_df_nga %>%
    left_join(cov_df) %>%
    mutate(coverage = ifelse(coveragemode == 'target', 0.8, vacc_cov_adj))


} #model2

yvals = c(0.62, 0.76, 0.77, 0.88, 0.77, 0.8, 1.07, 0.93)
yvals = (1 - yvals)
xvals = c(0.945, 0.88, 0.81, 0.51, 0.88, 0.82, 0.82, 0.865)
xydat = as.data.frame(cbind(xvals, yvals))

dat = as.data.frame(cbind(c(0, 0, 0, 1, 1, 1), c(1, 1, 1, 0.79, 0.74, 0.85), c('mean', 'lo', 'up', 'mean', 'lo', 'up')))
colnames(dat) <- c('coverage', 'effect', 'measure')
dat$effect <- (1 - as.numeric(dat$effect))
dat$coverage <- as.numeric(dat$coverage)
dat_wide = dat %>% pivot_wider(names_from = measure, values_from = effect)

p1 <- ggplot(data = subset(cases_df, coverage != 0 &
  scen == 'PMC-3' &
  age_group == 'U1' &
  name == 'PE_clinical_incidence')) +
  geom_ribbon(data = dat_wide, aes(y = mean, x = coverage, ymin = lo, ymax = up), fill = 'grey', alpha = 0.3) +
  geom_line(data = dat_wide, aes(y = mean, x = coverage), col = 'grey') +
  geom_line(aes(x = coverage, y = mean_val, col = scen), show.legend = F) +
  geom_ribbon(aes(x = coverage, ymin = low_val, ymax = up_val, fill = scen), show.legend = F, alpha = 0.3) +
  geom_point(aes(x = coverage, y = mean_val, col = scen, fill = scen, shape = 'generic')) +
  geom_point(data = subset(cases_df_nga, scen == 'PMC-3' & age_group == 'U1'),
             aes(x = coverage, y = PE_clinical_cases, fill = scen, shape = coveragemode),
             col = 'black', size = 2) +
  geom_point(data = xydat, aes(x = xvals, y = yvals), shape = 21, fill = 'black', size = 2) +
  scale_y_continuous(lim = c(0, 0.5)) +
  scale_shape_manual(values = c(21, 23, 25)) +
  labs(shape = 'Model type', col = 'Scenario', fill = 'Scenario') +
  facet_wrap(~age_group) +
  f_getCustomTheme()


temp_df <- subset(cases_df, (scen != 'PMC-5' | age_group == 'U2') &
  coverage != 0 &
  name == 'PE_clinical_incidence')
temp_df_nga <- subset(cases_df_nga, (scen != 'PMC-5' | age_group == 'U2'))

p2 <- ggplot(data = temp_df) +
  geom_line(aes(x = coverage, y = mean_val, col = scen), show.legend = F) +
  geom_ribbon(aes(x = coverage, ymin = low_val, ymax = up_val, fill = scen), show.legend = F, alpha = 0.3) +
  geom_point(aes(x = coverage, y = mean_val, col = scen, fill = scen, shape = 'generic')) +
  geom_point(data = temp_df_nga,
             aes(x = coverage, y = PE_clinical_cases, fill = scen, shape = coveragemode),
             col = 'black', size = 2) +
  scale_y_continuous(lim = c(0, 0.5)) +
  scale_shape_manual(values = c(21, 23, 25)) +
  labs(shape = 'Model type', col = 'Scenario', fill = 'Scenario') +
  facet_wrap(~age_group) +
  f_getCustomTheme()

pplot <- plot_grid(p1, p2, labels = c('a', 'b'), ncol = 2, rel_widths = c(0.7, 1))
pplot

f_save_plot(pplot, paste0('Fig A1.2.11'), file.path(plot_dir),
            width = 12, height = 4, units = 'in', device_format = device_format)
