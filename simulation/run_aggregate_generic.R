pckg <- c("data.table", "wrapr", "dplyr", "tidyr", "lubridate", "zoo")
a <- lapply(pckg, require, character.only = TRUE)
rm(a)


load_generic_simdat_week <- function(exp_name, simout_dir, sweepVars = NULL) {
  if (is.null(sweepVars))sweepVars <- c('pmc_mode', 'pmc_coverage', 'rtss_coverage', 'cm_coverage')

  grpVARS <- c('age', 'week', 'year', sweepVars,
               'Annual_EIR', 'seasonality', 'Cohort_birth_month', 'Run_Number')

  dat <- fread(file.path(simout_dir, exp_name, 'All_Age_weekly_Cases.csv'))

  if ('ipti_mode' %in% colnames(dat)) {
    dat$pmc_mode <- dat$ipti_mode
    dat$pmc_coverage <- dat$ipti_coverage
  }
  if (('date' %in% colnames(dat)))dat <- dat %>% mutate(date = as.Date(date))
  if (!('year' %in% colnames(dat)))dat <- dat %>% mutate(year = year(date) - 2019)
  print(summary(dat$year))

  if (!('age' %in% colnames(dat)))dat <- dat %>% mutate(age = as.numeric(date - min(date)))
  print(summary(dat$age))

  if (!('week' %in% colnames(dat)))dat <- dat %>% mutate(week = week(date))
  print(summary(dat$week))

  dat <- dat %>%
    group_by(year) %>%
    filter(date < max(date)) %>%    # need to remove last day of year due to aggregation issue otherwise (dump in cases)
    #filter(age > min(age)) %>%
    ungroup() %>%
    mutate(clinical_cases_pppa = New_Clinical_Cases / (Statistical_Population / 52),
           severe_cases_pppa = New_Severe_Cases / (Statistical_Population / 52),
           clinical_cases = New_Clinical_Cases / (Statistical_Population / 52) * 1000,
           severe_cases = New_Severe_Cases / (Statistical_Population / 52) * 1000,
           pmc_rtss_cov = paste0(pmc_coverage, '-', rtss_coverage))

  exp_key <- grpVARS[!(grpVARS %in% c('pmc_mode', 'pmc_coverage', 'rtss_coverage', 'pmc_rtss_cov'))]
  dat <- f_impact_measures(dat, exp_key)

  return(dat)
}

load_generic_simdat_month <- function(exp_name, simout_dir, sweepVars = NULL) {
  if (is.null(sweepVars))sweepVars <- c('pmc_mode', 'pmc_coverage', 'rtss_coverage', 'cm_coverage')

  grpVARS <- c('age', 'month', 'year', sweepVars,
               'Annual_EIR', 'seasonality', 'Cohort_birth_month', 'Run_Number')

  dat <- fread(file.path(simout_dir, exp_name, 'All_Age_monthly_Cases.csv')) %>%
    rename_with(~gsub(" ", "_", .x))

  if ('ipti_mode' %in% colnames(dat)) {
    dat$pmc_mode <- dat$ipti_mode
    dat$pmc_coverage <- dat$ipti_coverage
  }
  if (('date' %in% colnames(dat)))dat <- dat %>% mutate(date = as.Date(date))
  if (!('year' %in% colnames(dat)))dat <- dat %>% mutate(year = year(date) - 2019)
  print(summary(dat$year))

  if (!('age' %in% colnames(dat)))dat <- dat %>% mutate(age = as.numeric(date - min(date)))
  print(summary(dat$age))

  if (!('month' %in% colnames(dat)))dat <- dat %>% mutate(month = month(date))
  print(summary(dat$month))

  dat <- dat %>%
    mutate(clinical_cases_pppa = New_Clinical_Cases / (Statistical_Population / 12),
           severe_cases_pppa = New_Severe_Cases / (Statistical_Population / 12),
           clinical_cases = New_Clinical_Cases / (Statistical_Population / 12) * 1000,
           severe_cases = New_Severe_Cases / (Statistical_Population / 12) * 1000,
           pmc_rtss_cov = paste0(pmc_coverage, '-', rtss_coverage))

  exp_key <- grpVARS[!(grpVARS %in% c('pmc_mode', 'pmc_coverage', 'rtss_coverage', 'pmc_rtss_cov'))]
  dat <- f_impact_measures(dat, exp_key)

  return(dat)
}

load_generic_simdat_year <- function(exp_name, simout_dir, sweepVars = NULL) {

  if (is.null(sweepVars))sweepVars <- c('pmc_mode', 'pmc_coverage', 'rtss_coverage', 'cm_coverage')

  grpVARS <- c('year', sweepVars, 'Annual_EIR', 'seasonality', 'Cohort_birth_month', 'Run_Number')


  dat <- fread(file.path(simout_dir, exp_name, 'All_Age_yearly_Cases.csv'))

  if ('ipti_mode' %in% colnames(dat)) {
    dat$pmc_mode <- dat$ipti_mode
    dat$pmc_coverage <- dat$ipti_coverage
  }

  dat <- dat %>%
    rename_with(~gsub(" ", "_", .x)) %>%
    dplyr::group_by_at(.vars = grpVARS) %>%
    mutate(age = year * 365,
           clinical_cases_pppa = New_Clinical_Cases / Statistical_Population,
           severe_cases_pppa = New_Severe_Cases / Statistical_Population,
           clinical_cases = New_Clinical_Cases / Statistical_Population * 1000,
           severe_cases = New_Severe_Cases / Statistical_Population * 1000,
           pmc_rtss_cov = paste0(pmc_coverage, '-', rtss_coverage))

  grpVARS <- c(grpVARS, 'age')
  exp_key <- grpVARS[!(grpVARS %in% c('pmc_mode', 'pmc_coverage', 'rtss_coverage', 'pmc_rtss_cov'))]
  dat <- f_impact_measures(dat, exp_key)

  return(dat)

}

load_generic_simdat_agegroup <- function(exp_name, simout_dir, sweepVars = NULL, max_years = c(1, 2, 5, 10)) {

  if (is.null(sweepVars))sweepVars <- c('pmc_mode', 'pmc_coverage', 'rtss_coverage', 'cm_coverage')

  grpVARS <- c(sweepVars, 'Annual_EIR', 'seasonality', 'Cohort_birth_month', 'Run_Number')

  dat <- fread(file.path(simout_dir, exp_name, 'All_Age_yearly_Cases.csv'))

  if ('ipti_mode' %in% colnames(dat)) {
    dat$pmc_mode <- dat$ipti_mode
    dat$pmc_coverage <- dat$ipti_coverage
  }

  dat <- dat %>%
    rename_with(~gsub(" ", "_", .x)) %>%
    mutate(pmc_rtss_cov = paste0(pmc_coverage, '-', rtss_coverage))

  if (!('year' %in% colnames(dat)))dat <- dat %>% mutate(year = year(date) - 2020)
  print(summary(dat$year))

  if (max(dat$year) < max(max_years)) max_years <- max_years[max_years <= max(dat$year)]
  grpVARS <- c(grpVARS, 'pmc_rtss_cov')
  age_groups <- paste0("U", max_years)
  df_list <- list()
  for (i_max in max_years) {
    tdf <- dat %>%
      filter(year < i_max) %>%
      ungroup() %>%
      mutate(n_age = n_distinct(age)) %>%
      dplyr::group_by_at(.vars = c(grpVARS, 'n_age')) %>%
      dplyr::summarise(PfHRP2_Prevalence = mean(PfHRP2_Prevalence, na.rm = TRUE),
                       New_Clinical_Cases = sum(New_Clinical_Cases, na.rm = TRUE),
                       New_Severe_Cases = sum(New_Severe_Cases, na.rm = TRUE),
                       Statistical_Population = mean(Statistical_Population, na.rm = TRUE)) %>%
      mutate(clinical_cases_pppa = New_Clinical_Cases / (Statistical_Population * n_age),
             severe_cases_pppa = New_Severe_Cases / (Statistical_Population * n_age),
             clinical_cases = New_Clinical_Cases / (Statistical_Population * n_age) * 1000,
             severe_cases = New_Severe_Cases / (Statistical_Population * n_age) * 1000,
             clinical_cases_cum = New_Clinical_Cases / Statistical_Population * 1000,
             severe_cases_cum = New_Severe_Cases / Statistical_Population * 1000) %>%
      mutate(age_group = paste0('U', i_max)) %>%
      dplyr::ungroup()
    df_list[[length(df_list) + 1]] <- tdf
  }

  dat <- bind_rows(df_list)
  dat$age_group <- factor(dat$age_group, levels = age_groups, labels = age_groups)

  grpVARS <- c(grpVARS, 'age_group')
  exp_key <- grpVARS[!(grpVARS %in% c('Statistical_Population', 'pmc_mode', 'pmc_coverage', 'rtss_coverage', 'pmc_rtss_cov'))]
  dat <- f_impact_measures(dat, exp_key)

  return(dat)

}

f_impact_measures <- function(dat, exp_key) {
  #exp_key <- grpVARS[!(grpVARS %in% c('pmc_mode', 'pmc_coverage', 'rtss_coverage'))]
  if ('0-0' %in% unique(dat$pmc_rtss_cov)) {
    dat <- data.table(dat, key = exp_key)
    dat[, PE_clinical_incidence := 1 - clinical_cases / clinical_cases[pmc_rtss_cov == '0-0'], by = exp_key]
    dat[, PE_severe_incidence := 1 - severe_cases / severe_cases[pmc_rtss_cov == '0-0'], by = exp_key]
    dat[, clinical_cases_averted := clinical_cases[pmc_rtss_cov == '0-0'] - clinical_cases, by = exp_key]
    dat[, severe_cases_averted := severe_cases[pmc_rtss_cov == '0-0'] - severe_cases, by = exp_key]
    if ('clinical_cases_cum' %in% colnames(dat)) {
      dat[, clinical_cases_averted_cum := clinical_cases_cum[pmc_rtss_cov == '0-0'] - clinical_cases_cum, by = exp_key]
      dat[, severe_cases_averted_cum := severe_cases_cum[pmc_rtss_cov == '0-0'] - severe_cases_cum, by = exp_key]
    }
  }
  return(dat)
}


aggregate_dat <- function(dat, grp_vars, aggr_vars, outcome_vars) {
  if ('0-0' %in% unique(dat$pmc_rtss_cov)) {
    outcome_vars <- outcome_vars[outcome_vars %in% colnames(dat)]
  }

  dat <- dat %>%
    ungroup() %>%
    dplyr::select_at(.vars = c(grp_vars, aggr_vars, outcome_vars)) %>%
    tidyr::pivot_longer(col = outcome_vars) %>%
    dplyr::group_by_at(.vars = c(grp_vars, 'name')) %>%
    dplyr::summarise(mean_val = mean(value),
                     median_val = median(value),
                     low_val = quantile(value, probs = 0.05, na.rm = TRUE),
                     up_val = quantile(value, probs = 0.95, na.rm = TRUE),
                     min_val = min(value),
                     max_val = max(value)) %>%
    ungroup()

  return(dat)
}


##______________________________ run script
args <- commandArgs(TRUE)
exp_name <- args[1]
print(exp_name)

projectpath = os.getcwd()
simout_dir <- file.path(projectpath, 'simulation_output')
processout_dir <- file.path(projectpath, 'postprocessed_output', exp_name)

if (!dir.exists(processout_dir))dir.create(processout_dir)

sweepVars_access <- c('rtss_target_group', 'cm_target_group', 'pmc_target_group', 'frac_high_access')
sweepVars <- c('pmc_mode', 'pmc_coverage', 'rtss_coverage', 'cm_coverage') #,sweepVars_access
grp_vars <- c('age', 'year', 'Annual_EIR', sweepVars, 'pmc_rtss_cov')
aggr_vars <- c('Run_Number', 'Cohort_birth_month', 'seasonality')
outcome_vars <- c('Statistical_Population', 'PfHRP2_Prevalence',
                  'New_Clinical_Cases', 'New_Severe_Cases',
                  'clinical_cases', 'severe_cases',
                  'clinical_cases_pppa', 'severe_cases_pppa',
                  'PE_clinical_incidence', 'PE_severe_incidence',
                  'clinical_cases_averted', 'severe_cases_averted')


if (file.exists(file.path(simout_dir, exp_name, 'All_Age_yearly_Cases.csv'))) {
  cat('running load_generic_simdat_year')
  load_generic_simdat_year(exp_name, simout_dir, sweepVars) %>%
    aggregate_dat(grp_vars, aggr_vars, outcome_vars) %>%
    fwrite(file.path(processout_dir, 'simdat_aggr_year.csv'))
}

if (file.exists(file.path(simout_dir, exp_name, 'All_Age_monthly_Cases.csv'))) {
  cat('running load_generic_simdat_month')
  load_generic_simdat_month(exp_name, simout_dir, sweepVars) %>%
    aggregate_dat(grp_vars, aggr_vars, outcome_vars) %>%
    fwrite(file.path(processout_dir, 'simdat_aggr_month.csv'))
}

if (file.exists(file.path(simout_dir, exp_name, 'All_Age_weekly_Cases.csv'))) {
  cat('running load_generic_simdat_week')
  load_generic_simdat_week(exp_name, simout_dir, sweepVars) %>%
    aggregate_dat(grp_vars, aggr_vars, outcome_vars) %>%
    fwrite(file.path(processout_dir, 'simdat_aggr_week.csv'))
}


if (file.exists(file.path(simout_dir, exp_name, 'All_Age_yearly_Cases.csv'))) {
  cat('running load_generic_simdat_agegroup')
  grp_vars <- c('age_group', 'Annual_EIR', sweepVars, 'pmc_rtss_cov')
  outcome_vars <- unique(c(outcome_vars, 'clinical_cases_cum', 'severe_cases_cum',
                           'clinical_cases_averted_cum', 'severe_cases_averted_cum'))
  load_generic_simdat_agegroup(exp_name, simout_dir, sweepVars) %>%
    aggregate_dat(grp_vars, aggr_vars, outcome_vars) %>%
    fwrite(file.path(processout_dir, 'simdat_aggr_agegroup.csv'))
}
