##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## Fig_A1.02_supp_age_incidence.R
##---------------------
source(file.path('analysis', '_config.R'))

exp_name <- 'generic_PMC_RTSS_EIR_vaccSP_IIV'


cases_df <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_week.csv')) %>%
  rename_with(~gsub('ipti', 'pmc', .x)) %>%
  filter(Annual_EIR < eir_max) %>%
  mutate(pmc_mode_fct = ifelse(pmc_rtss_cov == '0-0', 'None',
                               ifelse(pmc_rtss_cov == '0-0.8', 'RTS,S',
                                      ifelse(pmc_rtss_cov == '0.8-0', 'PMC-3', 'PMC-3 + RTS,S'))))

cases_df$pmc_mode_fct <- factor(cases_df$pmc_mode_fct,
                                levels = (c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')),
                                labels = (c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')))



