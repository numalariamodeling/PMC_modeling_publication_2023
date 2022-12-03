#library('rdhs')
#dhs_tags()

##---------------
## PfPR
##---------------
dhs_data(tagIds = 36, countryIds = c("NG"), breakdown = "subnational", surveyYearStart = 2018) %>%
  dhs_get_name1() %>%
  filter(Indicator %in% c("Malaria prevalence according to microscopy", "Malaria prevalence according to RDT")) %>%
  dplyr::select(SurveyYear, Indicator, Value, NAME_1) %>%
  mutate(Indicator = gsub("Malaria prevalence according to ", "", Indicator)) %>%
  pivot_wider(names_from = "Indicator", values_from = "Value") %>%
  pivot_wider(names_from = "SurveyYear", values_from = c("RDT", "microscopy")) %>%
  add_nga_admin() %>%
  filter(southern_nga == 'Southern States') %>%
  mutate(RDT_2018 = RDT_2018 / 100,
         microscopy_2018 = microscopy_2018 / 100) %>%
  dplyr::select(NAME_1, RDT_2018, microscopy_2018) %>%
  rename(ADM1_NAME = NAME_1) %>%
  mutate(State = ADM1_NAME, ADM1_NAME = ifelse(ADM1_NAME == 'Akwa Ibom', 'Akwa lbom', ADM1_NAME)) %>%
  arrange(State) %>%
  fwrite(file = file.path("dhs_pfpr_df.csv"))


##---------------
## CM
##---------------
#'ML_FEVR_C_FEV',
dhs_data(indicatorIds = c('ML_FEVR_C_FEV', 'ML_FEVT_C_ADV', 'ML_FEVT_C_ACT'), countryIds = c("NG"), breakdown = "subnational", surveyYearStart = 2018) %>%
  dhs_get_name1() %>%
  dplyr::select(SurveyYear, Indicator, Value, NAME_1) %>%
  mutate(Indicator = gsub("Children with fever who took a combination with artemisinin", "GotACT", Indicator)) %>%
  mutate(Indicator = gsub("Children with fever for whom advice or treatment was sought", "SoughtCare", Indicator)) %>%
  mutate(Indicator = gsub("Children under 5 with fever in the last two weeks", "HadFever", Indicator)) %>%
  pivot_wider(names_from = "Indicator", values_from = "Value") %>%
  add_nga_admin() %>%
  filter(southern_nga == 'Southern States') %>%
  mutate(HadFever = HadFever / 100,
         SoughtCare = SoughtCare / 100,
         GotACT = GotACT / 100) %>%
  dplyr::select(NAME_1, HadFever, SoughtCare, GotACT) %>%
  rename(ADM1_NAME = NAME_1) %>%
  mutate(State = ADM1_NAME, ADM1_NAME = ifelse(ADM1_NAME == 'Akwa Ibom', 'Akwa lbom', ADM1_NAME)) %>%
  arrange(State) %>%
  fwrite(file = file.path("dhs_act_df.csv"))

fread(file.path(pmc_path, 'simulation_inputs/nigeria/projection_csvs/CM', 'HS_by_LGA_v5.csv')) %>%
  mutate(NAME_1 = State) %>%
  filter(year == 2018) %>%
  group_by(NAME_1) %>%
  summarize(U5_coverage = mean(U5_coverage)) %>%
  add_nga_admin() %>%
  filter(southern_nga == 'Southern States') %>%
  dplyr::select(NAME_1, U5_coverage) %>%
  mutate(State = NAME_1, NAME_1 = ifelse(NAME_1 == 'Akwa Ibom', 'Akwa lbom', NAME_1)) %>%
  arrange(State) %>%
  fwrite(file = file.path("dhs_act_df.csv"))

##---------------
## EPI
##---------------
dhs_data(tagIds = 32, countryIds = c("NG"), breakdown = "subnational", surveyYearStart = 2018) %>%
  dhs_get_name1() %>%
  dplyr::select(SurveyYear, Indicator, Value, NAME_1) %>%
  mutate(Indicator = gsub(" vaccination", "", Indicator),
         Indicator = tolower(gsub(' ', '_', Indicator))) %>%
  filter(!(grepl('hepatitis', Indicator)) &
           !(grepl('haemophilus', Indicator)) &
           !(grepl('polio', Indicator))) %>%
  pivot_wider(names_from = "Indicator", values_from = "Value") %>%
  pivot_longer(cols = -c('SurveyYear', 'NAME_1', 'received_all_8_basics', 'received_nos',
                         'number_of_children_12-23_months', 'number_of_children_12-23_months_(unweighted)'),
               names_to = 'vacc_round') %>%
  mutate(vacc_round = gsub('_received', '', vacc_round),
         vacc_age = case_when(vacc_round == 'bcg' ~ 0,
                              vacc_round == 'dpt_1' ~ round(2 * (365 / 12), 0),
                              vacc_round == 'dpt_2' ~ round(2.5 * (365 / 12), 0),
                              vacc_round == 'dpt_3' ~ round(3.5 * (365 / 12), 0),
                              vacc_round == 'measles' ~ round(9 * (365 / 12), 0)),
         pmc_3_recommen = ifelse(vacc_round %in% c('dpt_2', 'dpt_3', 'measles'), 'yes', 'no')) %>%
  separate(vacc_round, into = c('vacc', 'round'), sep = '_', remove = F) %>%
  add_nga_admin() %>%
  group_by(geo_zone, NAME_1) %>%
  mutate(mealses_incr = ifelse(value[vacc_round == 'measles'] > value[vacc_round == 'dpt_3'], 'yes', 'no'),
         dtp1_to_3_drop = round(value[vacc_round == 'dpt_3'] / value[vacc_round == 'dpt_1'], 2),
         dtp2_to_measles_drop = round(value[vacc_round == 'measles'] / value[vacc_round == 'dpt_2'], 2)) %>%
  rename(ADM1_NAME = NAME_1) %>%
  mutate(State = ADM1_NAME, ADM1_NAME = ifelse(ADM1_NAME == 'Akwa Ibom', 'Akwa lbom', ADM1_NAME)) %>%
  arrange(State) %>%
  fwrite(file = file.path("dhs_vacc_df.csv"))

##---------------
## EPI pattern
##---------------
exploreEPIpattern <- F
if (exploreEPIpattern) {
  dhs_vacc_df <- fread(file.path("dhs_vacc_df.csv")) %>% filter(southern_nga == 'Southern States')

  p1 <- ggplot(data = dhs_vacc_df) +
    geom_line(aes(x = vacc_age / 30, y = value, group = ADM1_NAME, col = mealses_incr)) +
    geom_point(aes(x = vacc_age / 30, y = value, group = ADM1_NAME, col = mealses_incr)) +
    facet_wrap(~geo_zone) +
    scale_x_continuous(breaks = seq(0, 12, 3)) +
    scale_y_continuous(breaks = seq(0, 100, 20)) +
    labs(x = 'age (months)', col = 'Coverage increase 9mth') +
    theme_cowplot()
  print(p1)

  summary(dhs_vacc_df$dtp1_to_3_drop)
  summary(dhs_vacc_df$dtp2_to_measles_drop)

  p2 <- dhs_vacc_df %>%
    mutate(dtp1_to_3_drop_grp = case_when(dtp1_to_3_drop <= 0.75 ~ '<=0.75',
                                          dtp1_to_3_drop > 0.75 ~ '>0.75')) %>%
    ggplot() +
    geom_line(aes(x = vacc_age / 30, y = value, group = ADM1_NAME, col = mealses_incr)) +
    geom_point(aes(x = vacc_age / 30, y = value, group = ADM1_NAME, col = mealses_incr)) +
    facet_wrap(~dtp1_to_3_drop_grp) +
    scale_x_continuous(breaks = seq(0, 12, 3)) +
    scale_y_continuous(breaks = seq(0, 100, 20)) +
    labs(x = 'age (months)', col = 'Coverage increase 9mth') +
    theme_cowplot()
  print(p2)

  p3 <- dhs_vacc_df %>%
    mutate(dtp2_to_measles_drop_grp = case_when(dtp2_to_measles_drop <= 0.9 ~ '<=0.9',
                                                dtp2_to_measles_drop > 0.9 & dtp2_to_measles_drop < 1 ~ '>0.9',
                                                dtp2_to_measles_drop > 1 ~ '>1')) %>%
    filter(vacc_round %in% c('dpt_2', 'dpt_3', 'measles')) %>%
    group_by(ADM1_NAME) %>%
    mutate(value_scl = value - value[vacc_round == 'dpt_2']) %>%
    ggplot() +
    geom_line(aes(x = vacc_age / 30, y = value_scl, group = ADM1_NAME, col = mealses_incr)) +
    geom_point(aes(x = vacc_age / 30, y = value_scl, group = ADM1_NAME, col = mealses_incr)) +
    # geom_smooth(aes(x = vacc_age / 30, y = value, group = mealses_incr),se=FALSE) +
    facet_wrap(~dtp2_to_measles_drop_grp) +
    scale_x_continuous(breaks = seq(0, 12, 3)) +
    scale_y_continuous(breaks = seq(0, 100, 20)) +
    labs(x = 'age (months)', col = 'Coverage increase 9mth') +
    theme_cowplot()
  print(p3)

  p4 <- dhs_vacc_df %>%
    mutate(dtp2_to_measles_drop_grp = case_when(dtp2_to_measles_drop <= 0.9 ~ '<=0.9',
                                                dtp2_to_measles_drop > 0.9 & dtp2_to_measles_drop < 1 ~ '>0.9',
                                                dtp2_to_measles_drop > 1 ~ '>1')) %>%
    filter(vacc_round %in% c('dpt_2', 'dpt_3', 'measles')) %>%
    group_by(vacc_age, vacc_round, dtp2_to_measles_drop_grp, mealses_incr) %>%
    summarize(value = mean(value)) %>%
    ungroup() %>%
    mutate(value_scl = value - value[vacc_round == 'dpt_2']) %>%
    ggplot() +
    geom_line(aes(x = vacc_round, y = value_scl, col = mealses_incr, group = interaction(mealses_incr, dtp2_to_measles_drop_grp))) +
    geom_point(aes(x = vacc_round, y = value_scl, col = mealses_incr, group = interaction(mealses_incr, dtp2_to_measles_drop_grp))) +
    scale_y_continuous(breaks = seq(0, 100, 20)) +
    labs(x = 'vacc_round', col = '') +
    #facet_wrap(~dtp2_to_measles_drop_grp) +
    theme_cowplot()

  print(p4)

}

##---------------
## EPI to PMC downscaling
##---------------
dhs_vacc_df <- fread(file.path("dhs_vacc_df.csv"))
sclalefactors <- EPI_to_pmc_cov(get_scaling_factors = T)
dhs_pmc_cov_df <- dhs_vacc_df %>%
  rename(vacc_cov = value) %>%
  filter(pmc_3_recommen == 'yes') %>%
  filter(southern_nga == 'Southern States') %>%
  mutate(vacc_cov = vacc_cov / 100,
         vacc_cov_adj = case_when(vacc_round == 'dpt_2' ~ vacc_cov * sclalefactors[1],
                                  vacc_round == 'dpt_3' ~ vacc_cov * sclalefactors[2],
                                  vacc_round == 'measles' ~ vacc_cov * sclalefactors[3]),
         vacc_cov_adj = ifelse(vacc_cov_adj < 0.2, 0.2, vacc_cov_adj)) %>%
  group_by(ADM1_NAME, State) %>%
  mutate(mean_vacc_cov = mean(vacc_cov),
         mean_vacc_cov_adj = mean(vacc_cov_adj)) %>%
  fwrite(file.path('dhs_epi_dat.csv'))


########################
touchpoints <- list()
touchpoints[['3tp']] <- round(c(10 / 4, 14 / 4, 9) * (365 / 12), 0) # WHO recommended
touchpoints[['5tp2ndyr']] <- sort(round(c(touchpoints[['3tp']], c(12, 15) * (365 / 12)), 0)) # explorative
touchpoints[['7tp']] <- sort(round(c(touchpoints[['3tp']], c(6, 12, 15, 18) * (365 / 12)), 0)) # explorative 12, 15, 18 as in Greenwood et al 2021
pmc_mode = '7tp'

## Adjust for additional doses
dhs_cov_dat <- fread(file.path('dhs_epi_dat.csv')) %>%
  dplyr::select(State, vacc_age, vacc_cov, vacc_cov_adj) %>%
  group_by(vacc_age) %>%
  mutate(round = cur_group_id(),
         epi_based = 1)

additional_vacc_ages <- touchpoints[[pmc_mode]][!(touchpoints[[pmc_mode]] %in% unique(dhs_cov_dat$vacc_age))]
i = 0
for (add_age in additional_vacc_ages) {
  i = i + 1
  tmp_df <- dhs_cov_dat %>%
    filter(round == 3) %>%
    mutate(round = 3 + i,
           vacc_age = add_age,
           epi_based = 0)
  dhs_cov_dat <- dhs_cov_dat %>% bind_rows(tmp_df)
}
dhs_cov_dat <- dhs_cov_dat %>%
  arrange(vacc_age) %>%
  group_by(vacc_age) %>%
  mutate(round_new = cur_group_id()) %>%
  group_by(State) %>%
  mutate(round_fix = ifelse(round_new != round, round_new, NA),
         vacc_cov = ifelse(vacc_age < max(touchpoints[['3tp']]) & !is.na(round_fix), ((vacc_cov[round_fix - 1] + vacc_cov[round_fix + 1]) / 2), vacc_cov),
         vacc_cov_adj = ifelse(vacc_age < max(touchpoints[['3tp']]) & !is.na(round_fix), ((vacc_cov_adj[round_fix - 1] + vacc_cov_adj[round_fix + 1]) / 2), vacc_cov_adj)) %>%
  group_by(vacc_age) %>%
  mutate(
    mean_vacc_cov = mean(vacc_cov),
    mean_vacc_cov_adj = mean(vacc_cov_adj))

target_dat <- dhs_cov_dat %>%
  dplyr::select(vacc_age, mean_vacc_cov) %>%
  unique() %>%
  pivot_longer(cols = -vacc_age) %>%
  mutate(value = 0.8, name = 'target')

pplot <- dhs_cov_dat %>%
  ungroup() %>%
  dplyr::select(vacc_age, mean_vacc_cov, mean_vacc_cov_adj) %>%
  unique() %>%
  pivot_longer(cols = -vacc_age) %>%
  #bind_rows(target_dat) %>%
  mutate(name = gsub('mean_vacc_', '', name)) %>%
  ggplot() +
  geom_hline(yintercept = 80) +
  geom_point(aes(x = vacc_age / (365 / 12), y = value * 100, fill = name), shape = 21, size = 1) +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_brewer(palette = 'Dark2') +
  labs(fill = '')
f_save_plot(pplot, plot_name = paste0('fig_5_coverage'), width = 3, height = 2, plot_dir)

pplot <- dhs_cov_dat %>%
  dplyr::select(State, vacc_age, vacc_cov, vacc_cov_adj) %>%
  pivot_longer(cols = c(-State, -vacc_age)) %>%
  mutate(name = gsub('vacc_', '', name)) %>%
  ggplot() +
  geom_hline(yintercept = 80) +
  geom_point(aes(x = vacc_age / (365 / 12), y = value * 100, fill = name), shape = 21, size = 1) +
  scale_x_continuous(breaks = seq(0, 24, 3)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_brewer(palette = 'Dark2') +
  labs(fill = '') +
  facet_wrap(~State)


f_save_plot(pplot, plot_name = paste0('fig_5_coverage_state'), width = 3, height = 2, plot_dir)
fwrite(dhs_cov_dat, file.path('dhs_epi_dat_extended.csv'))

#### Seasonality


##-----------------------------------------------
###### Create input df
##-----------------------------------------------
pfpr_df <- fread(file = file.path("dhs_pfpr_eir_df.csv")) %>%
  dplyr::select(State, RDT_2018, microscopy_2018, predEIR) %>%
  rename(annual_EIR = predEIR)

cm_df <- fread(file = file.path("dhs_act_df.csv")) %>%
  dplyr::select(State, U5_coverage) %>%
  mutate(State = gsub('Akwa lbom', 'Akwa Ibom', State)) %>%
  rename(cm_coverage = U5_coverage)

fread(file = file.path("dhs_act_df.csv")) %>%
  dplyr::select(State, U5_coverage) %>%
  mutate(State = gsub('Akwa lbom', 'Akwa Ibom', State)) %>%
  mutate(adult_coverage = 0.625,
         severe_cases = 0.8,
         simday = 0,
         duration = -1,
         run_col = 'run') %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/CM/NGA_CM.csv'))


epi_df <- fread(file = file.path("dhs_epi_dat_extended.csv")) %>%
  dplyr::select(State, vacc_age, round_new, vacc_cov_adj) %>%
  rename(round = round_new, pmc_coverage = vacc_cov_adj)
#mutate(round = paste0('pmc_', round))
#pivot_wider(names_from = round_new, values_from = vacc_cov_adj)


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
#mutate(round = paste0('pmc_', round))
#pivot_wider(names_from = round_new, values_from = vacc_cov_adj)


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
      mutate(pmc_coverage =0.8) %>%
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
    mutate(pmc_coverage =0.8) %>%
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
    mutate(pmc_coverage =0.8) %>%
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
    mutate(coverage_levels =0.8) %>%
  fwrite(file.path(pmc_path, 'simulation_inputs/generic_cohort_ds/RTSS/NGA_rtss_80.csv'))


