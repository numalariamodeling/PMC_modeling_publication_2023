##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## Fig_A1.2.9_PMCmode_cov_EIR.R
##---------------------

source(file.path('analysis', '_config.R'))

scenario_labels <- c('None', 'PMC-3', 'PMC-4', 'PMC-5', 'PMC-6', 'PMC-7', 'RTS,S', 'PMC-3 + RTS,S')
agegrp <- c('U2')
outcome_cols <- c('clinical_cases_averted', 'severe_cases_averted')
pmc_modes <- c("3tp", "4tp2ndyr", "5tp2ndyr", "6tp2ndyr", "7tp2ndyr")  # "4tp"  "5tp"
pmc_labels <- c('PMC-3', 'PMC-4', 'PMC-5', 'PMC-6', 'PMC-7')

eir_levels <- c(4, 8, 16, 32, 64, 128)
exp_name_eir <- c(paste0('generic_PMCmode_RTSS_cov_EIR', eir_levels, '_vaccSP_IIV'))


loadDat = T
if (loadDat) {

  cases_df_list <- list()
  for (exp_name in exp_name_eir) {

    cases_df_pmc <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
      filter(age_group %in% agegrp) %>%
      filter(pmc_mode %in% pmc_modes,
             pmc_coverage > 0 & rtss_coverage == 0) %>%
      mutate(coverage = pmc_coverage) %>%
      dplyr::select(-pmc_coverage, -pmc_rtss_cov)

    cases_df_rtss <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
      filter(age_group %in% agegrp,
             pmc_coverage == 0 & rtss_coverage > 0) %>%
      mutate(coverage = rtss_coverage,
             pmc_mode = 'rtss') %>%
      dplyr::select(-pmc_coverage, -pmc_rtss_cov)

    cases_df_0 <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
      filter(age_group %in% agegrp) %>%
      filter(pmc_mode %in% c('3tp'),
             pmc_coverage == 0 & rtss_coverage == 0) %>%
      mutate(coverage = 0) %>%
      dplyr::select(-pmc_coverage, -pmc_rtss_cov, -pmc_mode) %>%
      left_join(unique(cases_df_pmc[, c('name', 'Annual_EIR', 'age_group', 'pmc_mode')])) ## repeat by pmc_mode

    cases_df_pmctss <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
      filter(age_group %in% agegrp) %>%
      filter(pmc_mode %in% c('3tp'),
             pmc_coverage == rtss_coverage) %>%
      mutate(coverage = rtss_coverage,
             pmc_mode = 'rtss_pmc3') %>%
      dplyr::select(-pmc_coverage, -pmc_rtss_cov)

    cases_df_list[[length(cases_df_list) + 1]] <- cases_df_pmc %>%
      bind_rows(cases_df_0) %>%
      bind_rows(cases_df_pmctss)

    if (!nrow(cases_df_rtss) == 0) cases_df_list[[length(cases_df_list) + 1]] <- cases_df_rtss

  }
  cases_df <- cases_df_list %>% bind_rows()

}


cases_df$scen <- factor(cases_df$pmc_mode,
                        levels = c(pmc_modes, 'rtss', 'rtss_pmc3'),
                        labels = c(pmc_labels, 'RTS,S', 'PMC-3 + RTS,S'))

table(cases_df$scen, cases_df$Annual_EIR)
#table(cases_df$scen , cases_df$coverage)

plot_df_pppa <- cases_df %>%
  filter(name %in% c('clinical_cases_pppa', 'severe_cases_pppa')) %>%
  mutate(name = ifelse(name == 'clinical_cases_pppa', 'clinical', 'severe'))


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

plot_df_c <- plot_df %>%
  filter(name == 'clinical' & !(pmc_mode %in% c('rtss', 'rtss_pmc3'))) %>%
  group_by(Annual_EIR) %>%
  mutate(yconversion = max(yconversion))
plot_df_s <- plot_df %>%
  filter(name == 'severe' & !(pmc_mode %in% c('rtss', 'rtss_pmc3'))) %>%
  group_by(Annual_EIR) %>%
  mutate(yconversion = max(yconversion))

ggplot(data = plot_df_c) +
  geom_line(aes(x = coverage, y = mean_valPE, col = as.factor(Annual_EIR)), show.legend = F) +
  geom_ribbon(aes(x = coverage, ymin = low_valPE, ymax = up_valPE, fill = as.factor(Annual_EIR)), show.legend = F, alpha = 0.3) +
  geom_point(aes(x = coverage, y = mean_valPE, col = as.factor(Annual_EIR), fill = as.factor(Annual_EIR), shape = 'generic')) +
  scale_shape_manual(values = c(21, 23, 25)) +
  facet_wrap(~scen) +
  labs(y = 'Cases averted',
       col = 'Annual_EIR', fill = 'Annual_EIR') +
  customTheme_nogrid

p1 <- ggplot(data = plot_df_c) +
  geom_line(aes(x = coverage, y = mean_val, col = as.factor(Annual_EIR)), show.legend = F) +
  geom_ribbon(aes(x = coverage, ymin = low_val, ymax = up_val, fill = as.factor(Annual_EIR)), show.legend = F, alpha = 0.3) +
  # geom_point(aes(x = coverage, y = mean_val, col = as.factor(Annual_EIR), fill = as.factor(Annual_EIR)), shape = 21) +
  facet_wrap(~scen, nrow = 1) +
  scale_y_continuous(lim = c(NA, 1500), "Annual cases averted\nper 1000 population", sec.axis = sec_axis(~. / 3100 * 100, name = "PE (% reduction)")) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25) * 100) +
  scale_color_manual(values = getPalette) +
  scale_fill_manual(values = getPalette) +
  labs(y = 'Cases averted', x = 'Coverage (%)',
       col = 'Annual_EIR', fill = 'Annual_EIR') +
  customTheme_nogrid +
  guides(colour = guide_legend(reverse = T),
         fill = guide_legend(reverse = T))

p2 <- ggplot(data = plot_df_s) +
  geom_line(aes(x = coverage, y = mean_val, col = as.factor(Annual_EIR)), show.legend = F) +
  geom_ribbon(aes(x = coverage, ymin = low_val, ymax = up_val, fill = as.factor(Annual_EIR)), show.legend = F, alpha = 0.3) +
  #geom_point(aes(x = coverage, y = mean_val, col = as.factor(Annual_EIR), fill = as.factor(Annual_EIR)), shape = 21) +
  scale_y_continuous(lim = c(NA, 30), "Annual cases averted\nper 1000 population", sec.axis = sec_axis(~. / 30 * 100, name = "PE (% reduction)")) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25) * 100) +
  facet_wrap(~scen, nrow = 1) +
  scale_color_manual(values = getPalette) +
  scale_fill_manual(values = getPalette) +
  labs(y = 'Cases averted', x = 'Coverage (%)',
       col = 'Annual_EIR', fill = 'Annual_EIR') +
  customTheme_nogrid +
  guides(colour = guide_legend(reverse = T),
         fill = guide_legend(reverse = T))


pplot <- plot_combine(list(p1, p2), labels = c('a', 'b'), ncol = 1)
pplot

f_save_plot(pplot, paste0('Fig A1.2.9'), file.path(plot_dir),
            width = 10, height = 6, units = 'in', device_format = device_format)



