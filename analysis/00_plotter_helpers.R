library(dplyr)
library(data.table)

len <- function(x) { return(length(x)) }

f_save_plot <- function(pplot, plot_name, plot_dir, width = 14, height = 8, units = 'in', device_format = c('pdf', 'png')) {
  if ('png' %in% device_format) {
    ggsave(paste0(plot_name, ".png"), plot = pplot, path = plot_dir,
           width = width, height = height, units = units, device = "png")
  }
  if ('pdf' %in% device_format) {
    if (!dir.exists(file.path(plot_dir, "pdf"))) { dir.create(file.path(plot_dir, "pdf")) }
    ggsave(paste0(plot_name, ".pdf"), plot = pplot, path = file.path(plot_dir, "pdf"),
           width = width, height = height, units = units, device = "pdf", useDingbats = FALSE)
  }
}

plot_combine <- function(plist, ncol = 1, legend_position = 'right', rel_dims = 1, leg_dim = 0.25, labels = c('', '')) {
  plegend <- get_legend(plist[[1]])

  for (i in c(1:length(plist))) {
    plist[[i]] <- plist[[i]] + theme(legend.position = 'None')
  }

  if (sum(lengths(rel_dims)) == 1)rel_dims <- rep(rel_dims, length(plist))
  pplot <- cowplot::plot_grid(plotlist = plist, ncol = ncol, labels = labels, rel_widths = rel_dims)
  if (tolower(legend_position) == 'right')pplot <- plot_grid(pplot, plegend, ncol = 2, rel_widths = c(1, leg_dim))
  if (tolower(legend_position) == 'left')pplot <- plot_grid(plegend, pplot, ncol = 1, rel_heights = c(leg_dim, 1))
  if (tolower(legend_position) == 'bottom')pplot <- plot_grid(pplot, plegend, ncol = 1, rel_heights = c(1, leg_dim))
  if (tolower(legend_position) == 'top')pplot <- plot_grid(plegend, pplot, ncol = 1, rel_heights = c(leg_dim, 1))
  if (tolower(legend_position) == 'none')pplot <- pplot
  return(pplot)
}


f_getCustomTheme <- function(fontscl = 1) {

  customTheme <- theme(
    strip.text.x = element_text(size = 12 * fontscl, face = "bold"),
    strip.text.y = element_text(size = 12 * fontscl, face = "bold"),
    strip.background = element_blank(),
    plot.title = element_text(size = 14 * fontscl, vjust = -1, hjust = 0, color='black'),
    plot.subtitle = element_text(size = 12 * fontscl, color='black'),
    plot.caption = element_text(size = 9 * fontscl, color='black'),
    legend.title = element_text(size = 12 * fontscl, color='black'),
    legend.text = element_text(size = 12 * fontscl, color='black'),
    axis.title.x = element_text(size = 12 * fontscl, color='black'),
    axis.text.x = element_text(size = 12 * fontscl, color='black'),
    axis.title.y = element_text(size = 12 * fontscl, color='black'),
    axis.text.y = element_text(size = 12 * fontscl, color='black'),
    axis.ticks = element_line(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )

  return(customTheme)
}


f_add_variables <- function(cur_df) {

  if ('age_group' %in% colnames(cur_df)) {
    cur_df$age_group <- factor(cur_df$age_group,
                               levels = age_groups, labels = age_groups)
  }

  cur_df$clinical_cases_per100000 <- cur_df$clinical_cases * 100
  cur_df$severe_cases_per100000 <- cur_df$severe_cases * 100
  cur_df$cases_averted <- cur_df$cases_averted_per100000 / 100
  cur_df$severe_cases_averted <- cur_df$severe_cases_averted_per100000 / 100
  cur_df$rtss_cases_averted <- cur_df$rtss_cases_averted_per100000 / 100
  cur_df$rtss_severe_cases_averted <- cur_df$rtss_severe_cases_averted_per100000 / 100


  cur_df$ipti_mode_age <- factor(cur_df$ipti_mode,
                                 levels = c('3tp', '5tp', '5tp2ndyr', '6tp2ndyr', '7tp2ndyr'),
                                 labels = c('10w, 14w, 9m',
                                            '+ 6m, 12m',
                                            '+ 12m, 15m',
                                            '+ 6m, 12m, 15m',
                                            '+ 6m, 12m, 15m, 18m'))

  cur_df$ipti_mode_fct <- factor(cur_df$ipti_mode,
                                 levels = c('3tp', '5tp', '5tp2ndyr', '6tp2ndyr', '7tp2ndyr'),
                                 labels = c('PMC-3',
                                            'PMC-5a',
                                            'PMC-5b',
                                            'PMC-6',
                                            'PMC-7'))

  cur_df$ipti_rtss_grp <- as.character(cur_df$ipti_mode_fct)
  cur_df$ipti_rtss_grp[cur_df$rtss_coverage == 0 &
                         cur_df$ipti_coverage == 0] <- 'None'
  cur_df$ipti_rtss_grp[cur_df$rtss_coverage != 0] <- paste0(cur_df$ipti_rtss_grp[cur_df$rtss_coverage != 0], ' + RTS,S')
  cur_df$ipti_rtss_grp[cur_df$rtss_coverage != 0 &
                         cur_df$ipti_coverage == 0] <- 'RTS,S'

  cur_df$ipti_rtss_grp2 <- 'PMC'
  cur_df$ipti_rtss_grp2[cur_df$rtss_coverage == 0 &
                          cur_df$ipti_coverage == 0] <- 'None'
  cur_df$ipti_rtss_grp2[cur_df$rtss_coverage != 0] <- paste0(cur_df$ipti_rtss_grp2[cur_df$rtss_coverage != 0], ' + RTS,S')
  cur_df$ipti_rtss_grp2[cur_df$rtss_coverage != 0 &
                          cur_df$ipti_coverage == 0] <- 'RTS,S'

  return(cur_df)
}

load_Age_monthly_Cases <- function(simout_dir, exp_name, exp_sweeps = NULL, add_PE_perAge = TRUE, max_years = c(1, 2, 5, 10),
                                   exp_names_counterfactual = NULL, aggregateAge = TRUE, keep_birth_month = FALSE,
                                   keep_Run_Numbers = FALSE, fname = 'All_Age_monthly_Cases.csv') {
  # modified from
  # https://github.com/numalariamodeling/rtss-scenarios/blob/main/simulation/generic/app/scripts/process_sim_outputs.R

  cases_df <- fread(file.path(simout_dir, exp_name, fname)) %>% rename_with(~gsub(" ", "_", .x))

  # When using All_Age_Annual_Cases as returned from IPTi postprocessing, year had already been generated
  if (fname != 'All_Age_Annual_Cases.csv') {
    cases_df <- cases_df %>%
      mutate(date = as.Date(date),
             year = lubridate::year(date),
             year = year - min(year, na.rm = TRUE))
  }

  exp_sweeps <- c(exp_sweeps, c('Scenario_id', 'Annual_EIR', 'seasonality', 'cm_coverage', 'ipti_coverage',
                                'rtss_coverage', 'smc_coverage', 'Cohort_birth_month', 'rtss_mode', 'minBoostAge',
                                'cm_target_group', 'smc_target_group', 'rtss_target_group', 'ipti_target_group'))
  if (any(!(exp_sweeps %in% colnames(cases_df)))) exp_sweeps <- exp_sweeps[-which(!(exp_sweeps %in% colnames(cases_df)))]

  # aggregate to the year level
  if (aggregateAge) {
    cases_df_annual = cases_df %>%
      dplyr::group_by_at(c(exp_sweeps, 'Run_Number', 'year')) %>%
      dplyr::summarise(total_cases = sum(New_Clinical_Cases),
                       total_severe_cases = sum(New_Severe_Cases),
                       received_treatment = sum(Received_Treatment),
                       received_severe_treatment = sum(Received_Severe_Treatment),
                       received_NMF_treatment = sum(Received_NMF_Treatment),
                       average_pfpr = mean(PfHRP2_Prevalence),
                       average_pop = mean(Statistical_Population)) %>%
      dplyr::ungroup()
  }else {
    ## keep same names
    cases_df_annual = cases_df %>%
      mutate(year = as.numeric(date - min(date)) / 365) %>%
      dplyr::group_by_at(c(exp_sweeps, 'Run_Number', 'year')) %>%
      dplyr::summarise(total_cases = sum(New_Clinical_Cases),
                       total_severe_cases = sum(New_Severe_Cases),
                       received_treatment = sum(Received_Treatment),
                       received_severe_treatment = sum(Received_Severe_Treatment),
                       received_NMF_treatment = sum(Received_NMF_Treatment),
                       average_pfpr = mean(PfHRP2_Prevalence),
                       average_pop = mean(Statistical_Population)) %>%
      dplyr::ungroup()
  }

  # get averages across runs and optional across birth cohorts; calculate rates per 1000 people per year
  if (!keep_birth_month & ('Cohort_birth_month' %in% exp_sweeps)) exp_sweeps <- exp_sweeps[-which(exp_sweeps == 'Cohort_birth_month')]
  if (keep_Run_Numbers)exp_sweeps <- c(exp_sweeps, 'Run_Number')
  cases_df_annual_ave = cases_df_annual %>%
    dplyr::group_by_at(c(exp_sweeps, 'year')) %>%
    dplyr::summarise(clinical_cases = mean(total_cases),
                     severe_cases = mean(total_severe_cases),
                     received_treatment = mean(received_treatment),
                     received_severe_treatment = mean(received_severe_treatment),
                     received_NMF_treatment = mean(received_NMF_treatment),
                     pfpr = mean(average_pfpr),
                     pop = mean(average_pop)
    ) %>%
    dplyr::ungroup() %>%
    mutate(clinical_cases = clinical_cases / pop * 1000,
           severe_cases = severe_cases / pop * 1000,
           received_treatment = received_treatment / pop * 1000,
           received_severe_treatment = received_severe_treatment / pop * 1000,
           received_NMF_treatment = received_NMF_treatment / pop * 1000
    )

  if (add_PE_perAge) {
    cases_against_reference <- get_PE_and_averted_estimates(dat = cases_df_annual_ave, exp_sweeps = exp_sweeps, max_years = max_years)
    cases_age_aggregated = cases_against_reference[[1]]
    cases_each_year = cases_against_reference[[2]]
    cases_age_aggregated_sum = cases_against_reference[[3]]
    simdat <- list(cases_df_annual, cases_df_annual_ave, cases_age_aggregated, cases_each_year, cases_age_aggregated_sum)
  }else {
    simdat <- list(cases_df_annual, cases_df_annual_ave)
  }
  return(simdat)

}


add_PfPR_2_10 = function(prev_df_ave, pfpr_colname, exp_sweeps) {
  # cut out ages under 2 and over 10
  df = prev_df_ave[prev_df_ave$year >= 2 & prev_df_ave$year < 10,]

  # aggregate across ages within each scenario
  if (any(!(exp_sweeps %in% colnames(df)))) exp_sweeps <- exp_sweeps[-which(!(exp_sweeps %in% colnames(df)))]
  df_ave = df %>%
    dplyr::group_by_at(c(exp_sweeps)) %>%
    dplyr::summarise(PfPR_2_10 = mean(get(pfpr_colname))) %>%
    dplyr::ungroup()

  # merge the 2-10 PfPR values into the main dataframe
  merged_df = merge(prev_df_ave, df_ave, by = exp_sweeps)
  return(merged_df)
}


get_PE_and_averted_estimates <- function(dat, exp_sweeps, max_years = c(1, 2, 5, 10), reference = 'cm_coverage') {
  # calculate PE and cases averted relative to matched CM-only scenario and to the matched no-RTS,S scenario

  # reference scenarios are those with one intervention only (i.e., CM-only)
  intervention_coverages <- c('cm_coverage', 'rtss_coverage', 'smc_coverage', 'ipti_coverage')
  if (any(!(intervention_coverages %in% colnames(dat)))) intervention_coverages <- intervention_coverages[-which(!(intervention_coverages %in% colnames(dat)))]
  intervention_coverages <- intervention_coverages[-which(intervention_coverages == reference)]

  # remove non-reference scenarios and associated columns
  remove_cols <- c(intervention_coverages, 'Scenario_id', 'rtss_mode', 'rtss_target_group', 'minBoostAge', 'smc_target_group', 'ipti_mode', 'ipti_target_group', 'frac_high_access')
  if (any(!remove_cols %in% colnames(dat))) remove_cols <- remove_cols[-which(!(remove_cols %in% colnames(dat)))]
  df_references = dat %>%
    filter(across(intervention_coverages, ~. == 0)) %>%
    dplyr::select(c(all_of(exp_sweeps), 'year', 'clinical_cases', 'severe_cases', 'pfpr',
                    'received_treatment', 'received_severe_treatment', 'received_NMF_treatment')) %>%
    dplyr::select(-all_of(remove_cols)) %>%
    rename(clinical_cases_ref = clinical_cases,
           severe_cases_ref = severe_cases,
           received_treatment_ref = received_treatment,
           received_severe_treatment_ref = received_severe_treatment,
           received_NMF_treatment_ref = received_NMF_treatment)
  df_references = add_PfPR_2_10(prev_df_ave = df_references, pfpr_colname = 'pfpr', exp_sweeps = exp_sweeps) %>%
    rename(pfpr_2_10_ref = PfPR_2_10) %>%
    dplyr::select(-pfpr)


  # reference scenarios are those without RTS,S
  # remove non-reference scenarios and associated columns
  remove_cols <- c('rtss_coverage', 'Scenario_id', 'rtss_mode', 'rtss_target_group', 'minBoostAge')
  if (any(!remove_cols %in% colnames(dat))) remove_cols <- remove_cols[-which(!(remove_cols %in% colnames(dat)))]
  df_no_rtss_references = dat %>%
    filter(rtss_coverage == 0) %>%
    dplyr::select(c(all_of(exp_sweeps), 'year', 'clinical_cases', 'severe_cases', 'pfpr',
                    'received_treatment', 'received_severe_treatment', 'received_NMF_treatment')) %>%
    dplyr::select(-all_of(remove_cols)) %>%
    rename(clinical_cases_no_rtss = clinical_cases,
           severe_cases_no_rtss = severe_cases,
           received_treatment_no_rtss = received_treatment,
           received_severe_treatment_no_rtss = received_severe_treatment,
           received_NMF_treatment_no_rtss = received_NMF_treatment)
  df_no_rtss_references = add_PfPR_2_10(prev_df_ave = df_no_rtss_references, pfpr_colname = 'pfpr', exp_sweeps = exp_sweeps) %>%
    rename(pfpr_2_10_no_rtss = PfPR_2_10) %>%
    dplyr::select(-pfpr)

  # merge the reference values into the data frame
  cases_scenarios_references <- dat %>% left_join(df_references)
  cases_scenarios_references <- cases_scenarios_references %>% left_join(df_no_rtss_references)


  # = = = = = 
  # calculate protective efficacy and cases averted for each age year
  # = = = = = 
  # calculate key simulation summary results
  cases_scenarios_each_year <- cases_scenarios_references %>%
    mutate(
      # calculate protective efficacy relative to CM-only (no RTS,S/SMC/IPTi) for each age group
      protective_efficacy = 1 - clinical_cases / clinical_cases_ref,
      protective_efficacy_severe = 1 - severe_cases / severe_cases_ref,
      # calculate protective efficacy relative to no-RTS,S for each age group
      rtss_protective_efficacy = 1 - clinical_cases / clinical_cases_no_rtss,
      rtss_protective_efficacy_severe = 1 - severe_cases / severe_cases_no_rtss,
      # calculate burden relative to no RTS,S/SMC/IPTi for each age group
      relative_burden = (clinical_cases - clinical_cases_ref) / clinical_cases_ref,
      relative_burden_severe = (severe_cases - severe_cases_ref) / severe_cases_ref,
      # calculate burden relative to no RTS,S for each age group
      rtss_relative_burden = (clinical_cases - clinical_cases_no_rtss) / clinical_cases_no_rtss,
      rtss_relative_burden_severe = (severe_cases - severe_cases_no_rtss) / severe_cases_no_rtss,
      # cases and severe cases averted per 100,000 children, relative to CM-only scenario
      cases_averted_per100000 = (clinical_cases_ref - clinical_cases) * 100,
      severe_cases_averted_per100000 = (severe_cases_ref - severe_cases) * 100,
      # cases and severe cases averted per 100,000 children, relative to no-RTS,S scenario
      rtss_cases_averted_per100000 = (clinical_cases_no_rtss - clinical_cases) * 100,
      rtss_severe_cases_averted_per100000 = (severe_cases_no_rtss - severe_cases) * 100,

      ### treatments
      protective_efficacy_treatment = 1 - received_treatment / received_treatment_ref,
      protective_efficacy_treatment_severe = 1 - received_severe_treatment / received_severe_treatment_ref,
      # calculate protective efficacy relative to no-RTS,S for each age group
      rtss_protective_efficacy_treatment = 1 - received_treatment / received_treatment_no_rtss,
      rtss_protective_efficacy_treatment_severe = 1 - received_severe_treatment / received_severe_treatment_no_rtss,
      # calculate burden relative to no RTS,S/SMC/IPTi for each age group
      relative_treatment = (received_treatment - received_treatment_ref) / received_treatment_ref,
      relative_treatment_severe = (received_severe_treatment - received_severe_treatment_ref) / received_severe_treatment_ref,
      # calculate burden relative to no RTS,S for each age group
      rtss_relative_treatment = (received_treatment - received_treatment_no_rtss) / received_treatment_no_rtss,
      rtss_relative_treatment_severe = (received_severe_treatment - received_severe_treatment_no_rtss) / received_severe_treatment_no_rtss,
      # cases and severe cases averted per 100,000 children, relative to CM-only scenario
      treatment_saved_per100000 = (received_treatment_ref - received_treatment) * 100,
      severe_treatment_saved_per100000 = (received_severe_treatment_ref - received_severe_treatment) * 100,
      # cases and severe cases averted per 100,000 children, relative to no-RTS,S scenario
      rtss_treatment_saved_per100000 = (received_treatment_no_rtss - received_treatment) * 100,
      rtss_severe_treatment_saved_per100000 = (received_severe_treatment_no_rtss - received_severe_treatment) * 100,

      age_group = paste(year, '-', (year + 1))
    )


  # = = = = =
  # calculate protective efficacy and cases averted for each age group
  # = = = = =
  # aggregate dataframes to U1, U2, U5, and U10
  cases_age_aggregated_list <- list()
  for (i_max in max_years) {
    # get cases per 1000 across all included ages
    tdf = cases_scenarios_references %>%
      filter(year < i_max) %>%
      ungroup() %>%
      dplyr::group_by_at(exp_sweeps) %>%
      dplyr::summarise(clinical_cases = mean(clinical_cases),
                       clinical_cases_ref = mean(clinical_cases_ref),
                       clinical_cases_no_rtss = mean(clinical_cases_no_rtss),
                       severe_cases = mean(severe_cases),
                       severe_cases_ref = mean(severe_cases_ref),
                       severe_cases_no_rtss = mean(severe_cases_no_rtss),
                       pfpr_2_10_no_rtss = mean(pfpr_2_10_no_rtss),
                       pfpr_2_10_ref = mean(pfpr_2_10_ref),
                       received_treatment = mean(received_treatment),
                       received_treatment_ref = mean(received_treatment_ref),
                       received_treatment_no_rtss = mean(received_treatment_no_rtss),
                       received_severe_treatment = mean(received_severe_treatment),
                       received_severe_treatment_ref = mean(received_severe_treatment_ref),
                       received_severe_treatment_no_rtss = mean(received_severe_treatment_no_rtss),
      ) %>%
      mutate(age_group = paste0('U', i_max)) %>%
      dplyr::ungroup()

    # calculate key simulation summary results
    tdf <- tdf %>%
      mutate(
        # calculate protective efficacy relative to CM-only (no RTS,S/SMC/IPTi) for each age group
        protective_efficacy = 1 - clinical_cases / clinical_cases_ref,
        protective_efficacy_severe = 1 - severe_cases / severe_cases_ref,
        # calculate protective efficacy relative to no-RTS,S for each age group
        rtss_protective_efficacy = 1 - clinical_cases / clinical_cases_no_rtss,
        rtss_protective_efficacy_severe = 1 - severe_cases / severe_cases_no_rtss,
        # calculate burden relative to no RTS,S/SMC/IPTi for each age group
        relative_burden = (clinical_cases - clinical_cases_ref) / clinical_cases_ref,
        relative_burden_severe = (severe_cases - severe_cases_ref) / severe_cases_ref,
        # calculate burden relative to no RTS,S for each age group
        rtss_relative_burden = (clinical_cases - clinical_cases_no_rtss) / clinical_cases_no_rtss,
        rtss_relative_burden_severe = (severe_cases - severe_cases_no_rtss) / severe_cases_no_rtss,
        # cases and severe cases averted per 100,000 children, relative to CM-only scenario
        cases_averted_per100000 = (clinical_cases_ref - clinical_cases) * 100,
        severe_cases_averted_per100000 = (severe_cases_ref - severe_cases) * 100,
        # cases and severe cases averted per 100,000 children, relative to no-RTS,S scenario
        rtss_cases_averted_per100000 = (clinical_cases_no_rtss - clinical_cases) * 100,
        rtss_severe_cases_averted_per100000 = (severe_cases_no_rtss - severe_cases) * 100,
        ### treatments
        protective_efficacy_treatment = 1 - received_treatment / received_treatment_ref,
        protective_efficacy_treatment_severe = 1 - received_severe_treatment / received_severe_treatment_ref,
        # calculate protective efficacy relative to no-RTS,S for each age group
        rtss_protective_efficacy_treatment = 1 - received_treatment / received_treatment_no_rtss,
        rtss_protective_efficacy_treatment_severe = 1 - received_severe_treatment / received_severe_treatment_no_rtss,
        # calculate burden relative to no RTS,S/SMC/IPTi for each age group
        relative_treatment = (received_treatment - received_treatment_ref) / received_treatment_ref,
        relative_treatment_severe = (received_severe_treatment - received_severe_treatment_ref) / received_severe_treatment_ref,
        # calculate burden relative to no RTS,S for each age group
        rtss_relative_treatment = (received_treatment - received_treatment_no_rtss) / received_treatment_no_rtss,
        rtss_relative_treatment_severe = (received_severe_treatment - received_severe_treatment_no_rtss) / received_severe_treatment_no_rtss,
        # cases and severe cases averted per 100,000 children, relative to CM-only scenario
        treatment_saved_per100000 = (received_treatment_ref - received_treatment) * 100,
        severe_treatment_saved_per100000 = (received_severe_treatment_ref - received_severe_treatment) * 100,
        # cases and severe cases averted per 100,000 children, relative to no-RTS,S scenario
        rtss_treatment_saved_per100000 = (received_treatment_no_rtss - received_treatment) * 100,
        rtss_severe_treatment_saved_per100000 = (received_severe_treatment_no_rtss - received_severe_treatment) * 100)
    cases_age_aggregated_list[[length(cases_age_aggregated_list) + 1]] <- tdf
  }

  cases_age_aggregated <- bind_rows(cases_age_aggregated_list)

  # = = = = =
  # calculate protective efficacy and cases averted for each age group (sum)
  # = = = = =
  # aggregate dataframes to U1, U2, U5, and U10
  cases_age_aggregated_sum_list <- list()
  for (i_max in max_years) {
    # get cases per 1000 across all included ages
    tdf = cases_scenarios_references %>%
      filter(year < i_max) %>%
      ungroup() %>%
      dplyr::group_by_at(exp_sweeps) %>%
      dplyr::summarise(clinical_cases = sum(clinical_cases),
                       clinical_cases_ref = sum(clinical_cases_ref),
                       clinical_cases_no_rtss = sum(clinical_cases_no_rtss),
                       severe_cases = sum(severe_cases),
                       severe_cases_ref = sum(severe_cases_ref),
                       severe_cases_no_rtss = sum(severe_cases_no_rtss),
                       pfpr_2_10_no_rtss = mean(pfpr_2_10_no_rtss),
                       pfpr_2_10_ref = mean(pfpr_2_10_ref),
                       received_treatment = sum(received_treatment),
                       received_treatment_ref = sum(received_treatment_ref),
                       received_treatment_no_rtss = sum(received_treatment_no_rtss),
                       received_severe_treatment = sum(received_severe_treatment),
                       received_severe_treatment_ref = sum(received_severe_treatment_ref),
                       received_severe_treatment_no_rtss = sum(received_severe_treatment_no_rtss),
      ) %>%
      mutate(age_group = paste0('U', i_max)) %>%
      dplyr::ungroup()

    # calculate key simulation summary results
    tdf <- tdf %>%
      mutate(
        # calculate protective efficacy relative to CM-only (no RTS,S/SMC/IPTi) for each age group
        protective_efficacy = 1 - clinical_cases / clinical_cases_ref,
        protective_efficacy_severe = 1 - severe_cases / severe_cases_ref,
        # calculate protective efficacy relative to no-RTS,S for each age group
        rtss_protective_efficacy = 1 - clinical_cases / clinical_cases_no_rtss,
        rtss_protective_efficacy_severe = 1 - severe_cases / severe_cases_no_rtss,
        # calculate burden relative to no RTS,S/SMC/IPTi for each age group
        relative_burden = (clinical_cases - clinical_cases_ref) / clinical_cases_ref,
        relative_burden_severe = (severe_cases - severe_cases_ref) / severe_cases_ref,
        # calculate burden relative to no RTS,S for each age group
        rtss_relative_burden = (clinical_cases - clinical_cases_no_rtss) / clinical_cases_no_rtss,
        rtss_relative_burden_severe = (severe_cases - severe_cases_no_rtss) / severe_cases_no_rtss,
        # cases and severe cases averted per 100,000 children, relative to CM-only scenario
        cases_averted_per100000 = (clinical_cases_ref - clinical_cases) * 100,
        severe_cases_averted_per100000 = (severe_cases_ref - severe_cases) * 100,
        # cases and severe cases averted per 100,000 children, relative to no-RTS,S scenario
        rtss_cases_averted_per100000 = (clinical_cases_no_rtss - clinical_cases) * 100,
        rtss_severe_cases_averted_per100000 = (severe_cases_no_rtss - severe_cases) * 100,
        ### treatments
        protective_efficacy_treatment = 1 - received_treatment / received_treatment_ref,
        protective_efficacy_treatment_severe = 1 - received_severe_treatment / received_severe_treatment_ref,
        # calculate protective efficacy relative to no-RTS,S for each age group
        rtss_protective_efficacy_treatment = 1 - received_treatment / received_treatment_no_rtss,
        rtss_protective_efficacy_treatment_severe = 1 - received_severe_treatment / received_severe_treatment_no_rtss,
        # calculate burden relative to no RTS,S/SMC/IPTi for each age group
        relative_treatment = (received_treatment - received_treatment_ref) / received_treatment_ref,
        relative_treatment_severe = (received_severe_treatment - received_severe_treatment_ref) / received_severe_treatment_ref,
        # calculate burden relative to no RTS,S for each age group
        rtss_relative_treatment = (received_treatment - received_treatment_no_rtss) / received_treatment_no_rtss,
        rtss_relative_treatment_severe = (received_severe_treatment - received_severe_treatment_no_rtss) / received_severe_treatment_no_rtss,
        # cases and severe cases averted per 100,000 children, relative to CM-only scenario
        treatment_saved_per100000 = (received_treatment_ref - received_treatment) * 100,
        severe_treatment_saved_per100000 = (received_severe_treatment_ref - received_severe_treatment) * 100,
        # cases and severe cases averted per 100,000 children, relative to no-RTS,S scenario
        rtss_treatment_saved_per100000 = (received_treatment_no_rtss - received_treatment) * 100,
        rtss_severe_treatment_saved_per100000 = (received_severe_treatment_no_rtss - received_severe_treatment) * 100)
    cases_age_aggregated_sum_list[[length(cases_age_aggregated_sum_list) + 1]] <- tdf
  }

  cases_age_aggregated_sum <- bind_rows(cases_age_aggregated_sum_list)


  return(list(cases_age_aggregated, cases_scenarios_each_year, cases_age_aggregated_sum))
}


calc_high_low_access_coverages = function(coverage_all, high_access_frac) {
  if ((high_access_frac < 1) & (coverage_all >= high_access_frac)) {
    coverage_high = 1
    coverage_low = (coverage_all - high_access_frac) / (1 - high_access_frac)
  } else {
    coverage_high = coverage_all / high_access_frac
    coverage_low = 0
  }
  return(c(coverage_high, coverage_low))
}


calc_frac_no_inter = function(inter_names, coverages, targeting, frac_in_high, createVenn = FALSE) {
  sample_pop = 10000
  high_access_nums = 1:round(sample_pop * frac_in_high)
  low_access_nums = round(sample_pop * frac_in_high + 1):sample_pop

  # assign each intervention to people in the population
  xx = list()
  for (ii in 1:length(inter_names)) {
    if (targeting[ii] == 'high') {
      high_low_cov = calc_high_low_access_coverages(coverage_all = coverages[ii], high_access_frac = frac_in_high)
      # sample from high-access
      cov_ids = sample(high_access_nums, size = high_low_cov[1] * length(high_access_nums), replace = FALSE)
      # sample from low-access
      cov_ids = c(cov_ids, sample(low_access_nums, size = high_low_cov[2] * length(low_access_nums), replace = FALSE))
    } else if (targeting[ii] == 'low') {
      low_high_cov = calc_high_low_access_coverages(coverage_all = coverages[ii], high_access_frac = (1 - frac_in_high))
      # sample from low-access
      cov_ids = sample(low_access_nums, size = low_high_cov[1] * length(low_access_nums), replace = FALSE)
      # sample from high-access
      cov_ids = c(cov_ids, sample(high_access_nums, size = low_high_cov[2] * length(high_access_nums), replace = FALSE))
    } else { # random
      cov_ids = sample(1:sample_pop, size = coverages[ii] * sample_pop, replace = FALSE)
    }
    xx[[ii]] = cov_ids
  }
  # add venn diagrams of children reached by each intervention
  if (createVenn) {
    library('ggVennDiagram')
    ggVennDiagram(xx,
                  label = "percent")
  }

  # get fraction not receiving any intervention
  inc_ids = unique(unlist(xx))
  frac_no_inter = (sample_pop - length(inc_ids)) / sample_pop


  # get fraction receiving either RTS,S or SMC
  inc_ids = unique(unlist(xx[which(inter_names %in% c('smc', 'rtss'))]))
  frac_rtss_or_smc = (length(inc_ids)) / sample_pop

  # get fraction receiving either RTS,S or SMC
  inc_ids = unique(unlist(xx[which(inter_names %in% c('cm', 'rtss'))]))
  frac_rtss_or_cm = (length(inc_ids)) / sample_pop

  # get fraction not previously covered by an intervention that are now covered by RTSS
  smc_cm_ids = unique(unlist(xx[which(inter_names != 'rtss')]))
  no_smc_cm_ids = which(!(1:sample_pop %in% smc_cm_ids))
  now_rtss = intersect(xx[[which(inter_names == 'rtss')]], no_smc_cm_ids)
  frac_newly_covered = length(now_rtss) / length(no_smc_cm_ids)


  # get fraction not previously covered by SMC that are now covered by RTSS
  smc_ids = unique(xx[[which(inter_names == 'smc')]])
  no_smc_ids = which(!(1:sample_pop %in% smc_ids))
  now_rtss = intersect(xx[[which(inter_names == 'rtss')]], no_smc_ids)
  frac_newly_covered_not_smc = length(now_rtss) / length(no_smc_ids)


  # get fraction not previously covered by CM that are now covered by RTSS
  cm_ids = unique(xx[[which(inter_names == 'cm')]])
  no_cm_ids = which(!(1:sample_pop %in% cm_ids))
  now_rtss = intersect(xx[[which(inter_names == 'rtss')]], no_cm_ids)
  frac_newly_covered_not_cm = length(now_rtss) / length(no_cm_ids)


  return(c(frac_no_inter, frac_newly_covered, frac_newly_covered_not_smc, frac_rtss_or_smc, frac_rtss_or_cm))
}


