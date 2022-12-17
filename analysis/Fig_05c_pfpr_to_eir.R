## Fig_05c_pfpr_to_eir.R

(exp_name <- 'generic_PMC3_RTSScov_EIR_constant_vaccSP_IIV')

baseline_pfpr <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
  filter(age_group %in% c('U5', 'U2') &
           pmc_coverage == 0 &
           rtss_coverage == 0) %>%
  filter(name %in% c('PfHRP2_Prevalence', 'clinical_cases', 'severe_cases')) %>%
  mutate(name = case_when(name == 'PfHRP2_Prevalence' ~ 'pfpr',
                          name == 'clinical_cases' ~ 'clinicalinc',
                          name == 'severe_cases' ~ 'severeinc')) %>%
  pivot_longer(cols = c(mean_val, median_val, low_val, up_val, min_val, max_val), names_to = 'statistic') %>%
  pivot_wider(names_from = age_group, values_from = value,
              names_glue = "{age_group}_{.value}") %>%
  pivot_wider(names_from = name,
              values_from = c(U2_value, U5_value),
              names_glue = "{name}_{.value}") %>%
  dplyr::select(-pmc_coverage, -rtss_coverage, -pmc_mode, -pmc_rtss_cov) %>%
  rename_with(~gsub("_value", "", .x))

#fwrite(baseline_pfpr, file.path(simout_dir, exp_name, 'baseline_pfpr.csv'))
baseline_pfpr <- fread(file.path(simout_dir, exp_name, 'baseline_pfpr.csv'))


## add country DHS data
pfpr_df <- fread(file = file.path("data_files", "dhs_pfpr_df.csv")) %>%
  dplyr::select(State, RDT_2018, microscopy_2018)

baseline_pfpr_sub <- baseline_pfpr %>%
  filter(cm_coverage == 0.4)

y <- log(baseline_pfpr_sub$Annual_EIR)
x <- baseline_pfpr_sub$pfpr_U5

model <- loess(y ~ x)
pred_df = data.frame('RDT_2018' = sort(pfpr_df$RDT_2018))
pred_df$predEIR <- exp(predict(model, pred_df$RDT_2018, type = 'response'))
pfpr_df <- pfpr_df %>% left_join(pred_df)

# ggplot(data = subset(baseline_pfpr, statistic == 'mean_val')) +
#   geom_point(aes(x = Annual_EIR, y = pfpr_U5, group = cm_coverage)) +
#   geom_line(aes(x = Annual_EIR, y = pfpr_U5, group = cm_coverage)) +
#   scale_x_log10() +
#   geom_point(data = pfpr_df, aes(x = predEIR, y = RDT_2018, group = State), col = 'deepskyblue3')

fwrite(pfpr_df, file.path("data_files", "dhs_pfpr_eir_df.csv"))

