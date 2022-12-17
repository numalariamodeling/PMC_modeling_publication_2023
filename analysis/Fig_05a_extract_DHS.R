## Fig_05a_extract_DHS.R
source(file.path('analysis', '_config.R'))
source(file.path('analysis', '_fig05_helper_functions.R'))
library('rdhs')


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
  fwrite(file = file.path("data_files", "dhs_pfpr_df.csv"))


##---------------
## CM
##---------------
#'ML_FEVR_C_FEV',
dhs_data(indicatorIds = c('ML_FEVR_C_FEV', 'ML_FEVT_C_ADV', 'ML_FEVT_C_ACT'),
         countryIds = c("NG"), breakdown = "subnational", surveyYearStart = 2018) %>%
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
  fwrite(file = file.path("data_files", "dhs_act_df.csv"))

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
  fwrite(file = file.path("data_files", "dhs_act_df.csv"))

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
  fwrite(file = file.path("data_files", "dhs_vacc_df.csv"))
