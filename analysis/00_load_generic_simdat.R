load_generic_simdat_month <- function(exp_name, simout_dir) {

  if (exp_name == 'generic_PMC_RTSS_EIR_IIV_SP') {
    dat <- fread(file.path(simout_dir, exp_name, 'All_Age_monthly_Cases.csv'))
  }
  if (exp_name == 'generic_PMC_RTSS_EIR_IIV_SDX_PYR') {
    dat <- fread(file.path(simout_dir, exp_name, 'All_Age_monthly_Cases_U5.csv')) %>%
      bind_rows(fread(file.path(simout_dir, exp_name, 'All_Age_monthly_Cases_5to10.csv')))
  }

  dat <- dat %>%
    rename_with(~gsub(" ", "_", .x)) %>%
    mutate(date = as.Date(date),
           year = year(date) - 2020) %>%
    filter(seasonality %in% default_seasons) %>%
    dplyr::group_by(date, year, ipti_mode, ipti_coverage, rtss_coverage, Annual_EIR, seasonality, Cohort_birth_month, Run_Number) %>%
    summarize(PfHRP2_Prevalence = mean(PfHRP2_Prevalence),
              New_Clinical_Cases = sum(New_Clinical_Cases),
              New_Severe_Cases = sum(New_Severe_Cases),
              Statistical_Population = mean(Statistical_Population)) %>%
    mutate(clinical_cases = New_Clinical_Cases / Statistical_Population * 1000,
           severe_cases = New_Severe_Cases / Statistical_Population * 1000,
           ipti_rtss_cov = paste0(ipti_coverage, '-', rtss_coverage))

  exp_key <- c('date', 'year', 'Annual_EIR', 'Run_Number', 'Cohort_birth_month', 'seasonality')  #'ipti_mode',
  dat <- data.table(dat, key = exp_key)
  dat[, PE_clinical_incidence := 1 - clinical_cases / clinical_cases[ipti_rtss_cov == '0-0'], by = exp_key]
  dat[, PE_severe_incidence := 1 - severe_cases / severe_cases[ipti_rtss_cov == '0-0'], by = exp_key]
  dat[, clinical_cases_averted := clinical_cases[ipti_rtss_cov == '0-0'] - clinical_cases, by = exp_key]
  dat[, severe_cases_averted := severe_cases[ipti_rtss_cov == '0-0'] - severe_cases, by = exp_key]
  dat[, clinical_incidence_red := clinical_cases[ipti_rtss_cov == '0-0'] - clinical_cases, by = exp_key]
  dat[, severe_incidence_red := severe_cases[ipti_rtss_cov == '0-0'] - severe_cases, by = exp_key]

  dat <- dat %>% mutate(ipti_mode_fct = ifelse(ipti_rtss_cov == '0-0', 'None',
                                               ifelse(ipti_rtss_cov == '0-0.8', 'RTS,S',
                                                      ifelse(ipti_rtss_cov == '0.8-0', 'PMC-3', 'PMC-3 + RTS,S'))))

  dat$ipti_mode_fct <- factor(dat$ipipti_mode_fctti_mode,
                              levels = rev(c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')),
                              labels = rev(c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')))
  return(dat)
}

load_generic_simdat_year <- function(exp_name, simout_dir, do_smooth) {

  if (exp_name != 'generic_PMC_RTSS_EIR_IIV_SDX_PYR') {
    dat <- fread(file.path(simout_dir, exp_name, 'All_Age_monthly_Cases.csv'))
  }
  if (exp_name == 'generic_PMC_RTSS_EIR_IIV_SDX_PYR') {
    dat <- fread(file.path(simout_dir, exp_name, 'All_Age_monthly_Cases_U5.csv')) %>%
      bind_rows(fread(file.path(simout_dir, exp_name, 'All_Age_monthly_Cases_5to10.csv')))
  }

  dat <- dat %>%
    rename_with(~gsub(" ", "_", .x)) %>%
    mutate(date = as.Date(date) - years(2020),
           year = year(date))

  if (do_smooth & 'date' %in% colnames(dat)) dat <- dat %>% mutate(date = lubridate::round_date(date, unit = '12 month'))

  dat <- dat %>%
    filter(seasonality %in% default_seasons) %>%
    dplyr::group_by(year, ipti_mode, ipti_coverage, rtss_coverage, Annual_EIR, seasonality, Cohort_birth_month, Run_Number) %>%
    summarize(PfHRP2_Prevalence = mean(PfHRP2_Prevalence),
              New_Clinical_Cases = sum(New_Clinical_Cases),
              New_Severe_Cases = sum(New_Severe_Cases),
              Statistical_Population = mean(Statistical_Population)) %>%
    mutate(clinical_cases = New_Clinical_Cases / Statistical_Population * 1000,
           severe_cases = New_Severe_Cases / Statistical_Population * 1000,
           ipti_rtss_cov = paste0(ipti_coverage, '-', rtss_coverage))

  exp_key <- c('year', 'Annual_EIR', 'Run_Number', 'Cohort_birth_month', 'seasonality')
  dat <- data.table(dat, key = exp_key)
  dat[, PE_clinical_incidence := 1 - clinical_cases / clinical_cases[ipti_rtss_cov == '0-0'], by = exp_key]
  dat[, PE_severe_incidence := 1 - severe_cases / severe_cases[ipti_rtss_cov == '0-0'], by = exp_key]
  dat[, clinical_cases_averted := clinical_cases[ipti_rtss_cov == '0-0'] - clinical_cases, by = exp_key]
  dat[, severe_cases_averted := severe_cases[ipti_rtss_cov == '0-0'] - severe_cases, by = exp_key]
  dat[, clinical_incidence_red := clinical_cases[ipti_rtss_cov == '0-0'] - clinical_cases, by = exp_key]
  dat[, severe_incidence_red := severe_cases[ipti_rtss_cov == '0-0'] - severe_cases, by = exp_key]

  dat <- dat %>% mutate(ipti_mode_fct = ifelse(ipti_rtss_cov == '0-0', 'None',
                                               ifelse(ipti_rtss_cov == '0-0.8', 'RTS,S',
                                                      ifelse(ipti_rtss_cov == '0.8-0', 'PMC-3', 'PMC-3 + RTS,S'))))

  dat$ipti_mode_fct <- factor(dat$ipti_mode_fct,
                              levels = rev(c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')),
                              labels = rev(c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')))

  return(dat)

}

load_generic_simdat_week <- function(exp_name, simout_dir) {

  dat <- fread(file.path(simout_dir, exp_name, 'All_Age_weekly_Cases.csv')) %>%
    rename_with(~gsub(" ", "_", .x)) %>%
    mutate(date = as.Date(date)) %>%
    # mutate( date =as.Date( lubridate::round_date(date, unit = '2 month'))) %>%
    #mutate( dateM =as.Date( lubridate::round_date(date, unit = '2 month')),
    #       date = ifelse(date <= as.Date('2020-12-01'), date, dateM)) %>%
    mutate(date = as.Date(date),
           week = week(date),
           year = year(date) - 2020) %>%
    group_by(year) %>%
    filter(week < max(week)) %>%
    ungroup() %>%
    filter(seasonality %in% default_seasons) %>%
    dplyr::group_by(date, year, ipti_mode, ipti_coverage, rtss_coverage, Annual_EIR, seasonality, Cohort_birth_month, Run_Number) %>%
    summarize(PfHRP2_Prevalence = mean(PfHRP2_Prevalence),
              New_Clinical_Cases = sum(New_Clinical_Cases),
              New_Severe_Cases = sum(New_Severe_Cases),
              Statistical_Population = mean(Statistical_Population)) %>%
    dplyr::group_by(date, ipti_mode, ipti_coverage, rtss_coverage, Annual_EIR, seasonality, Cohort_birth_month, Run_Number) %>%
    mutate(clinical_cases = New_Clinical_Cases / Statistical_Population * 1000,
           severe_cases = New_Severe_Cases / Statistical_Population * 1000,
           ipti_rtss_cov = paste0(ipti_coverage, '-', rtss_coverage))

  exp_key <- c('date', 'year', 'Annual_EIR', 'Run_Number', 'Cohort_birth_month', 'seasonality')
  dat <- data.table(dat, key = exp_key)
  dat[, PE_clinical_incidence := 1 - clinical_cases / clinical_cases[ipti_rtss_cov == '0-0'], by = exp_key]
  dat[, PE_severe_incidence := 1 - severe_cases / severe_cases[ipti_rtss_cov == '0-0'], by = exp_key]
  dat[, clinical_cases_averted := clinical_cases[ipti_rtss_cov == '0-0'] - clinical_cases, by = exp_key]
  dat[, severe_cases_averted := severe_cases[ipti_rtss_cov == '0-0'] - severe_cases, by = exp_key]
  dat[, clinical_incidence_red := clinical_cases[ipti_rtss_cov == '0-0'] - clinical_cases, by = exp_key]
  dat[, severe_incidence_red := severe_cases[ipti_rtss_cov == '0-0'] - severe_cases, by = exp_key]

  dat <- dat %>% mutate(ipti_mode_fct = ifelse(ipti_rtss_cov == '0-0', 'None',
                                               ifelse(ipti_rtss_cov == '0-0.8', 'RTS,S',
                                                      ifelse(ipti_rtss_cov == '0.8-0', 'PMC-3', 'PMC-3 + RTS,S'))))

  dat$ipti_mode_fct <- factor(dat$ipti_mode_fct,
                              levels = rev(c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')),
                              labels = rev(c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')))

  return(dat)
}

load_generic_simdat_agegroup <- function(exp_name, simout_dir) {

  if (exp_name == 'generic_PMC_RTSS_EIR_IIV_SP') {
    cases_df <- fread(file.path(simout_dir, exp_name, 'All_Age_monthly_Cases.csv'))
  }
  if (exp_name == 'generic_PMC_RTSS_EIR_IIV_SDX_PYR') {
    cases_df <- fread(file.path(simout_dir, exp_name, 'All_Age_monthly_Cases_U5.csv')) %>%
      bind_rows(fread(file.path(simout_dir, exp_name, 'All_Age_monthly_Cases_5to10.csv')))
  }

  # TODO use per week
  cases_df <- cases_df %>%
    rename_with(~gsub(" ", "_", .x)) %>%
    mutate(date = as.Date(date),
           year = year(date) - 2020) %>%
    dplyr::group_by(year, ipti_mode, ipti_coverage, rtss_coverage, Annual_EIR, seasonality, Cohort_birth_month, Run_Number) %>%
    summarize(PfHRP2_Prevalence = mean(PfHRP2_Prevalence),
              New_Clinical_Cases = sum(New_Clinical_Cases),
              New_Severe_Cases = sum(New_Severe_Cases),
              Statistical_Population = mean(Statistical_Population)) %>%
    mutate(clinical_cases = New_Clinical_Cases / Statistical_Population * 1000,
           severe_cases = New_Severe_Cases / Statistical_Population * 1000,
           ipti_rtss_cov = paste0(ipti_coverage, '-', rtss_coverage))

}

aggregate_agegroups <- function(dat, max_years, grp_vars) {
  # max_years <-c(1,2,5,10)
  #grp_vars <- c('Annual_EIR', 'ipti_mode', 'ipti_coverage', 'rtss_coverage', 'Run_Number', 'Cohort_birth_month', 'seasonality', 'ipti_rtss_cov')
  age_groups <- paste0("U", max_years)
  df_list <- list()
  for (i_max in max_years) {
    # get cases per 1000 across all included ages
    tdf = dat %>%
      filter(year < i_max) %>%
      ungroup() %>%
      dplyr::group_by_at(grp_vars) %>%
      dplyr::summarise(PfHRP2_Prevalence = mean(PfHRP2_Prevalence, na.rm = TRUE),
                       New_Clinical_Cases = sum(New_Clinical_Cases, na.rm = TRUE),
                       New_Severe_Cases = sum(New_Severe_Cases, na.rm = TRUE),
                       Statistical_Population = mean(Statistical_Population, na.rm = TRUE)) %>%
      mutate(clinical_cases = New_Clinical_Cases / Statistical_Population * 1000,
             severe_cases = New_Severe_Cases / Statistical_Population * 1000) %>%
      mutate(age_group = paste0('U', i_max)) %>%
      dplyr::ungroup()
    df_list[[length(df_list) + 1]] <- tdf
  }

  df_aggr <- bind_rows(df_list)
  table(df_aggr$age_group)
  df_aggr$age_group <- factor(df_aggr$age_group,
                              levels = age_groups,
                              labels = age_groups)


  exp_key <- c('age_group', 'Annual_EIR', 'Run_Number', 'Cohort_birth_month', 'seasonality')
  df_aggr <- data.table(df_aggr, key = exp_key)
  df_aggr[, PE_clinical_incidence := 1 - clinical_cases / clinical_cases[ipti_rtss_cov == '0-0'], by = exp_key]
  df_aggr[, PE_severe_incidence := 1 - severe_cases / severe_cases[ipti_rtss_cov == '0-0'], by = exp_key]
  df_aggr[, clinical_cases_averted := clinical_cases[ipti_rtss_cov == '0-0'] - clinical_cases, by = exp_key]
  df_aggr[, severe_cases_averted := severe_cases[ipti_rtss_cov == '0-0'] - severe_cases, by = exp_key]
  df_aggr[, clinical_incidence_red := clinical_cases[ipti_rtss_cov == '0-0'] - clinical_cases, by = exp_key]
  df_aggr[, severe_incidence_red := severe_cases[ipti_rtss_cov == '0-0'] - severe_cases, by = exp_key]

  return(df_aggr)
}

aggregate_dat <- function(dat, grp_vars, aggr_vars, outcome_var = NULL) {
  if (is.null(outcome_vars))outcome_vars <- c('Statistical_Population', 'PfHRP2_Prevalence',
                                              'New_Clinical_Cases', 'New_Severe_Cases',
                                              'clinical_cases', 'severe_cases',
                                              'PE_clinical_incidence', 'PE_severe_incidence',
                                              'clinical_cases_averted', 'severe_cases_averted')

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
