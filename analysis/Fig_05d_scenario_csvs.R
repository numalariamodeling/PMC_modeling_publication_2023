##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## Fig_05c_scenario_csvs.R
##---------------------

source(file.path('analysis', '_config.R'))
source(file.path('analysis', '_fig05_helper_functions.R'))

epi_df <- fread(file = file.path("data_files", "dhs_epi_dat_extended.csv")) %>%
  dplyr::select(State, vacc_age, round_new, vacc_cov_adj) %>%
  rename(round = round_new, pmc_coverage = vacc_cov_adj)

pfpr_df <- fread(file = file.path("data_files", "dhs_pfpr_eir_df.csv")) %>%
  dplyr::select(State, RDT_2018, microscopy_2018, predEIR) %>%
  rename(annual_EIR = predEIR)

# cm_df <- fread(file = file.path("data_files", "dhs_act_df.csv")) %>%
#   dplyr::select(State, U5_coverage) %>%
#   mutate(State = gsub('Akwa lbom', 'Akwa Ibom', State)) %>%
#   rename(cm_coverage = U5_coverage)

fread(file = file.path("data_files", "dhs_act_df.csv")) %>%
  dplyr::select(State, U5_coverage) %>%
  mutate(State = gsub('Akwa lbom', 'Akwa Ibom', State)) %>%
  mutate(adult_coverage = 0.625,
         severe_cases = 0.8,
         simday = 0,
         duration = -1,
         run_col = 'run') # %>%
#  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/CM/NGA_CM.csv'))

pfpr_df %>%
  left_join(cm_df) %>%
  mutate(pmc_coverage = 'custom',
         pmc_mode = '3tp',
         rtss_coverage = 0,
         CM_filename = 'CM/NGA_CM.csv',
         RTSS_filename = 'RTSS/RTSS_constant_0coverage.csv',
         PMC_filename = 'PMC/NGA_pmc_3tp.csv') %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/NGA_pmc_3tp.csv'))

pfpr_df %>%
  left_join(cm_df) %>%
  mutate(pmc_coverage = 'custom',
         pmc_mode = '7tp2ndyr',
         rtss_coverage = 0,
         CM_filename = 'CM/NGA_CM.csv',
         RTSS_filename = 'RTSS/RTSS_constant_0coverage.csv',
         PMC_filename = 'PMC/NGA_pmc_7tp2ndyr.csv') %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/NGA_pmc_7tp2ndyr.csv'))

pfpr_df %>%
  left_join(cm_df) %>%
  mutate(pmc_coverage = 'custom',
         pmc_mode = '5tp2ndyr',
         rtss_coverage = 0,
         CM_filename = 'CM/NGA_CM.csv',
         RTSS_filename = 'RTSS/RTSS_constant_0coverage.csv',
         PMC_filename = 'PMC/NGA_pmc_5tp2ndyr.csv') %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/NGA_pmc_5tp2ndyr.csv'))

### ith RTSS
pfpr_df %>%
  left_join(cm_df) %>%
  mutate(pmc_coverage = 'custom',
         pmc_mode = '3tp',
         rtss_coverage = 'custom',
         CM_filename = 'CM/NGA_CM.csv',
         RTSS_filename = 'RTSS/RTSS_constant_0coverage.csv',
         PMC_filename = 'PMC/NGA_pmc_3tp.csv') %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/NGA_rtss_pmc_3tp.csv'))

pfpr_df %>%
  left_join(cm_df) %>%
  mutate(pmc_coverage = 0,
         pmc_mode = '3tp',
         rtss_coverage = 'custom',
         CM_filename = 'CM/NGA_CM.csv',
         RTSS_filename = 'RTSS/NGA_rtss.csv',
         PMC_filename = 'PMC/PMC_3tp_0coverage.csv') %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/NGA_rtss.csv'))


#### sub dir files
pmc_df <- pfpr_df %>%
  dplyr::select(-RDT_2018, -microscopy_2018, -annual_EIR) %>%
  mutate(repetitions = 1,
         tsteps_btwn_repetitions = -1,
         agemin = 0, agemax = 2,
         deploy_type = 'EPI_cohort',
         run_col = 'run')

epi_df1 <- epi_df %>%
  filter(vacc_age %in% touchpoints[['3tp']]) %>%
  mutate(pmc_mode = '3tp') %>%
  group_by(vacc_age) %>%
  mutate(round = cur_group_id(),
         pmc_touchpoints = vacc_age,
         PMC_day = vacc_age) %>%
  left_join(pmc_df) %>%
  rename(coverage_levels = pmc_coverage) %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/PMC/NGA_pmc_3tp.csv'))


epi_df2 <- epi_df %>%
  filter(vacc_age %in% touchpoints[['5tp2ndyr']]) %>%
  mutate(pmc_mode = '5tp2ndyr') %>%
  group_by(vacc_age) %>%
  mutate(round = cur_group_id(),
         pmc_touchpoints = vacc_age,
         PMC_day = vacc_age) %>%
  left_join(pmc_df) %>%
  rename(coverage_levels = pmc_coverage) %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/PMC/NGA_pmc_5tp2ndyr.csv'))


epi_df3 <- epi_df %>%
  filter(vacc_age %in% touchpoints[['7tp']]) %>%
  mutate(pmc_mode = '7tp2ndyr') %>%
  group_by(vacc_age) %>%
  mutate(round = cur_group_id(),
         pmc_touchpoints = vacc_age,
         PMC_day = vacc_age) %>%
  left_join(pmc_df) %>%
  rename(coverage_levels = pmc_coverage) %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/PMC/NGA_pmc_7tp2ndyr.csv'))


## rtss
epi_df3 <- epi_df %>%
  filter(vacc_age %in% c(274, 548)) %>%
  mutate(rtss_types = ifelse(vacc_age == 274, 'simple', 'booster'),
         initial_killing = ifelse(vacc_age == 274, 0.8, 0.4),
         decay_time_constant = 592.4066512,
         RTSS_day = ifelse(vacc_age == 274, 274, 730)) %>%
  mutate(coverage_levels = ifelse(vacc_age == 274, pmc_coverage, 0.5)) %>%
  mutate(repetitions = 1,
         tsteps_btwn_repetitions = -1,
         agemin = 0, agemax = 5,
         distribution_name = 'CONSTANT_DISTRIBUTION',
         run_col = 'run') %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/RTSS/NGA_rtss.csv'))


##### OTHER

fread(file.path(pmc_path, 'SouthernNGA_extractedMonthlyEIR.csv')) %>%
  rename(month = Month) %>%
  group_by(State, month) %>%
  summarize(Monthly_EIR = mean(Monthly_EIR)) %>%
  group_by(State) %>%
  mutate(Monthly_EIR = Monthly_EIR / sum(Monthly_EIR)) %>%
  ungroup() %>%
  mutate(State = gsub('Akwa lbom', 'Akwa Ibom', State)) %>%
  pivot_wider(names_from = State, values_from = Monthly_EIR) %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/Seasonality/SouthernNGA_extractedMonthlyEIR.csv'))


season_df <- fread(file.path(pmc_path, 'SouthernNGA_extractedMonthlyEIR.csv')) %>%
  group_by(State, Month) %>%
  summarize(Monthly_EIR = mean(Monthly_EIR)) %>%
  group_by(State) %>%
  mutate(Monthly_EIR = Monthly_EIR / sum(Monthly_EIR)) %>%
  ungroup() %>%
  mutate(State = gsub('Akwa lbom', 'Akwa Ibom', State),
         Month = ifelse(Month < 10, paste0('0', Month), Month),
         Month = paste0('eir_', Month)) %>%
  pivot_wider(names_from = Month, values_from = Monthly_EIR)

country_df <- pfpr_df %>%
  left_join(cm_df) %>%
  left_join(epi_df) %>%
  left_join(season_df)

#fwrite(country_df, file.path('nga_input_df.csv'))


##### TARGET COVERAGE
epi_df <- fread(file = file.path("dhs_epi_dat_extended.csv")) %>%
  dplyr::select(State, vacc_age, round_new, vacc_cov_adj) %>%
  rename(round = round_new, pmc_coverage = vacc_cov_adj)

pfpr_df %>%
  left_join(cm_df) %>%
  mutate(pmc_coverage = 'custom',
         pmc_mode = '3tp',
         rtss_coverage = 0,
         CM_filename = 'CM/NGA_CM.csv',
         RTSS_filename = 'RTSS/RTSS_constant_0coverage.csv',
         PMC_filename = 'PMC/NGA_pmc_80_3tp.csv') %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/NGA_pmc_3tp_targetcov.csv'))

pfpr_df %>%
  left_join(cm_df) %>%
  mutate(pmc_coverage = 'custom',
         pmc_mode = '7tp2ndyr',
         rtss_coverage = 0,
         CM_filename = 'CM/NGA_CM.csv',
         RTSS_filename = 'RTSS/RTSS_constant_0coverage.csv',
         PMC_filename = 'PMC/NGA_pmc_80_7tp2ndyr.csv') %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/NGA_pmc_7tp2ndyr_targetcov.csv'))

pfpr_df %>%
  left_join(cm_df) %>%
  mutate(pmc_coverage = 'custom',
         pmc_mode = '5tp2ndyr',
         rtss_coverage = 0,
         CM_filename = 'CM/NGA_CM.csv',
         RTSS_filename = 'RTSS/RTSS_constant_0coverage.csv',
         PMC_filename = 'PMC/NGA_pmc_80_5tp2ndyr.csv') %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/NGA_pmc_5tp2ndyr_targetcov.csv'))

### ith RTSS
pfpr_df %>%
  left_join(cm_df) %>%
  mutate(pmc_coverage = 'custom',
         pmc_mode = '3tp',
         rtss_coverage = 'custom',
         CM_filename = 'CM/NGA_CM.csv',
         RTSS_filename = 'RTSS/RTSS_constant_0coverage.csv',
         PMC_filename = 'PMC/NGA_pmc_80_3tp.csv') %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/NGA_rtss_pmc_3tp_targetcov.csv'))

pfpr_df %>%
  left_join(cm_df) %>%
  mutate(pmc_coverage = 0,
         pmc_mode = '3tp',
         rtss_coverage = 'custom',
         CM_filename = 'CM/NGA_CM.csv',
         RTSS_filename = 'RTSS/NGA_rtss_80.csv',
         PMC_filename = 'PMC/PMC_3tp_0coverage.csv') %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/NGA_rtss_targetcov.csv'))


#### sub dir files
pmc_df <- pfpr_df %>%
  dplyr::select(-RDT_2018, -microscopy_2018, -annual_EIR) %>%
  mutate(repetitions = 1,
         tsteps_btwn_repetitions = -1,
         agemin = 0, agemax = 2,
         deploy_type = 'EPI_cohort',
         run_col = 'run')

epi_df1 <- epi_df %>%
  filter(vacc_age %in% touchpoints[['3tp']]) %>%
  mutate(pmc_mode = '3tp') %>%
  group_by(vacc_age) %>%
  mutate(round = cur_group_id(),
         pmc_touchpoints = vacc_age,
         PMC_day = vacc_age) %>%
  left_join(pmc_df) %>%
  mutate(pmc_coverage = 0.8) %>%
  rename(coverage_levels = pmc_coverage) %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/PMC/NGA_pmc_80_3tp.csv'))


epi_df2 <- epi_df %>%
  filter(vacc_age %in% touchpoints[['5tp2ndyr']]) %>%
  mutate(pmc_mode = '5tp2ndyr') %>%
  group_by(vacc_age) %>%
  mutate(round = cur_group_id(),
         pmc_touchpoints = vacc_age,
         PMC_day = vacc_age) %>%
  left_join(pmc_df) %>%
  mutate(pmc_coverage = 0.8) %>%
  rename(coverage_levels = pmc_coverage) %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/PMC/NGA_pmc_80_5tp2ndyr.csv'))


epi_df3 <- epi_df %>%
  filter(vacc_age %in% touchpoints[['7tp']]) %>%
  mutate(pmc_mode = '7tp2ndyr') %>%
  group_by(vacc_age) %>%
  mutate(round = cur_group_id(),
         pmc_touchpoints = vacc_age,
         PMC_day = vacc_age) %>%
  left_join(pmc_df) %>%
  mutate(pmc_coverage = 0.8) %>%
  rename(coverage_levels = pmc_coverage) %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/PMC/NGA_pmc_80_7tp2ndyr.csv'))


## rtss
epi_df3 <- epi_df %>%
  filter(vacc_age %in% c(274, 548)) %>%
  mutate(rtss_types = ifelse(vacc_age == 274, 'simple', 'booster'),
         initial_killing = ifelse(vacc_age == 274, 0.8, 0.4),
         decay_time_constant = 592.4066512,
         RTSS_day = ifelse(vacc_age == 274, 274, 730)) %>%
  mutate(coverage_levels = ifelse(vacc_age == 274, 0.8, 0.8)) %>%
  mutate(repetitions = 1,
         tsteps_btwn_repetitions = -1,
         agemin = 0, agemax = 5,
         distribution_name = 'CONSTANT_DISTRIBUTION',
         run_col = 'run') %>%
  mutate(coverage_levels = 0.8) %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/RTSS/NGA_rtss_80.csv'))


