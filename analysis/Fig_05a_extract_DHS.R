##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## Fig_05a_extract_DHS.R
##---------------------

source(file.path('analysis', '_config.R'))
source(file.path('analysis', '_fig05_helper_functions.R'))

##---------------
## PfPR
##---------------
dhs_data(tagIds = 36, countryIds = c("NG"), breakdown = "subnational", , surveyYearStart = 2018, surveyYearEnd= 2018) %>%
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
  fwrite(file = file.path("data_files", "ndhs_2018_pfpr_df.csv"))


##---------------
## CM
##---------------
## Used from Ozodiegwu et al 2022
### see https://github.com/numalariamodeling/hbhi-nigeria-publication-2021/blob/main/hbhi-dhs-tools/1_variables_scripts/CM/CM_DHS_estimates.R


##---------------
## EPI
##---------------
dhs_data(tagIds = 32, countryIds = c("NG"), breakdown = "subnational", surveyYearStart = 2018, surveyYearEnd= 2018) %>%
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
  fwrite(file = file.path("data_files", "ndhs_2018_vacc_df.csv"))

