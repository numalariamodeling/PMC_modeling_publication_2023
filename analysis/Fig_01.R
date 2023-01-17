## Fig_01.R
source(file.path('analysis', '_config.R'))

exp_name_pmc <- 'generic_single_PMC_vaccSP_IIV'
exp_name_rtss <- 'rtss_validation_KintampoPhase3_wBooster'

pmc_single <- TRUE
rtss_single <- TRUE
lexis_plot <- TRUE

##---------------------------------------
## A - PMC single dose
##---------------------------------------
if (pmc_single) {
  mc_efficacy_weekly <- as.data.frame(cbind(c(0:8), c(78.64, 78.30, 77.73, 74.77, 66.93, 51.02, 27.49, 7.60, 0.00)))
  colnames(mc_efficacy_weekly) <- c('week', 'PE')
  mc_efficacy_weekly$rel_red <- mc_efficacy_weekly$PE / 100
  mc_efficacy_weekly$age_days <- mc_efficacy_weekly$week * 7

  cases_df <- fread(file.path(simout_dir, exp_name_pmc, 'simdat_aggr_week.csv')) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    filter(cm_coverage == '0.6' & Annual_EIR == 32)

  ## Offset inlcuded in simulations 10 days
  t_offset_post = 7 # additional offset (see S1_pmc_single.R)
  pdat <- cases_df %>%
    filter(pmc_coverage == 1 & name == 'PE_clinical_incidence') %>%
    filter(age >= 76) %>%
    mutate(age = age - 76 - t_offset_post,
           time = age / 7)

  pplot <- ggplot(data = pdat) +
    geom_hline(yintercept = 0, size = 0.7) +
    geom_ribbon(data = pdat, aes(x = time, ymin = low_val, ymax = up_val), alpha = 0.3, fill = pmc_cols[1]) +
    geom_line(data = pdat, aes(x = time, y = mean_val, group = 1), col = pmc_cols[1], size = 1.1) +
    geom_point(data = mc_efficacy_weekly, aes(x = week, y = rel_red)) +
    scale_y_continuous(lim = c(-0.2, 1), breaks = seq(-0.2, 1, 0.2), labels = percent) +
    scale_x_continuous(lim = c(0, 12), breaks = c(0:12), labels = c(0:12)) +
    scale_color_manual(values = pmc_cols[1]) +
    scale_fill_manual(values = pmc_cols[1]) +
    labs(title = '',
         y = 'Modeled incidence\nreduction (%)',
         x = 'Weeks after single PMC dose') +
    customTheme_nogrid +
    theme(legend.position = 'None')

  print(pplot)
  f_save_plot(pplot, plot_name = paste0('Fig1A'), width = 6, height = 4, plot_dir = plot_dir)

  rm(cases_df, pplot, pdat, pdat1, pdat2, pdat3)
} # pmc_single


##---------------------------------------
## B  -  RTS,S single dose
##---------------------------------------
if (rtss_single) {
  vaccine_date <- as.Date('2021-01-01')

  reference_filepath <- file.path('data_files', 'Penny_et_al_2016_supplement_Kintampo.csv')
  ref_df <- fread(reference_filepath) %>%
    rename(time_group = bar, PE_ref = mean) %>%
    dplyr::select('time_group', 'PE_ref')

  cases_df <- fread(file.path('simulation_output', exp_name_rtss, 'All_Age_monthly_Cases.csv')) %>%
    rename_with(~gsub(' ', '_', .x)) %>%
    mutate(date = as.Date(date)) %>%
    filter(date >= vaccine_date) %>%  ## after vaccination date
    mutate(elapsed_months = round(as.numeric(difftime(date, vaccine_date, units = 'days')) / (365 / 12), 0),
           time_group = floor(elapsed_months / 3) + 1) %>%
    dplyr::group_by(Scenario_id, Run_Number, time_group, Annual_EIR, rtss_coverage) %>%
    dplyr::summarise(total_cases = sum(New_Clinical_Cases),
                     middle_month = mean(date)) %>%
    arrange(Run_Number, time_group, middle_month, Annual_EIR) %>%
    dplyr::group_by(Run_Number, time_group, middle_month, Annual_EIR) %>%
    mutate(PE_sim = 1 - (total_cases / total_cases[rtss_coverage == 0]),
           time_since_3rd_dose = as.numeric((middle_month - vaccine_date) / 365)) %>%
    arrange(middle_month)

  ## Aggregate runs
  cases_df <- cases_df %>%
    filter(rtss_coverage == 1) %>%
    dplyr::group_by(Annual_EIR, time_group, time_since_3rd_dose) %>%
    dplyr::summarize(mean_val = mean(PE_sim),
                     low_val = quantile(PE_sim, probs = 0.05, na.rm = TRUE),
                     up_val = quantile(PE_sim, probs = 0.95, na.rm = TRUE)) %>%
    arrange(time_since_3rd_dose)


  pplot <- cases_df %>%
    left_join(ref_df) %>%
    ggplot() +
    geom_hline(yintercept = 0, size = 0.7, linetype = 'dashed') +
    geom_ribbon(aes(x = time_since_3rd_dose, ymin = low_val, ymax = up_val), alpha = 0.3, fill = rtss_col) +
    geom_line(aes(x = time_since_3rd_dose, y = mean_val), col = rtss_col, size = 1.1) +
    geom_point(aes(x = time_since_3rd_dose, y = PE_ref)) +
    scale_y_continuous(lim = c(-0.50, 1),
                       breaks = seq(-0.50, 1, 0.25),
                       labels = seq(-0.50, 1, 0.25) * 100) +
    scale_x_continuous(breaks = c(0:5), labels = c(0:5)) +
    labs(title = '',
         y = 'Modeled incidence\nreduction(%)',
         x = 'Years after third RTS,S dose ') +
    customTheme

  print(pplot)
  f_save_plot(pplot, plot_name = 'Fig1B', width = 4, height = 3, plot_dir = plot_dir)

} # (rtss_single) {


##---------------------------------------
## C  -  Lexis plot
##---------------------------------------
if (lexis_plot) {

  y <- seq(0, 6, 0.01)
  x <- as.Date('0-01-01') + (y * 365)

  y_pmc_all_m <- c(2.5, 3.5, 6, 9, 12, 15, 18)
  y_pmc_all <- round((y_pmc_all_m * (365 / 12)) / 365, 2)
  y_pmc <- as.data.frame(cbind('y' = y_pmc_all, 'pmc_y' = y_pmc_all, 'pmc_m' = y_pmc_all_m))

  y_rtss_all_m <- c(6, 7.5, 9, 24)
  y_rtss_all <- round((y_rtss_all_m * (365 / 12)) / 365, 2)
  y_rtss <- as.data.frame(cbind('y' = y_rtss_all, 'rtss_y' = y_rtss_all, 'rtss_m' = y_rtss_all_m))

  dat <- as.data.frame(x)
  dat$y <- y

  for (i in c(1:11)) {
    dat[paste0('x', i)] = dat$x + days(30 * i)
  }

  dat$index = 1
  dat <- dat %>%
    left_join(y_pmc) %>%
    left_join(y_rtss) %>%
    as.data.frame()

  pts <- 1.5
  pts_shp <- 21

  pplot <- suppressWarnings(
    ggplot(data = subset(dat, x <= min(dat$x) + years(5))) +
      geom_line(aes(x = x, y = y), col = '#999999') +
      geom_point(aes(x = x, y = pmc_y), fill = pmc_cols[1], size = pts, shape = pts_shp) +
      geom_point(aes(x = x, y = rtss_y), fill = rtss_col, size = pts, shape = pts_shp) +
      scale_x_date(lim = c(as.Date('0-01-01'), as.Date('8-01-01')),
                   date_breaks = '1 year', date_labels = '%Y',
                   expand = c(0, 0), date_minor_breaks = '6 month') +
      scale_y_continuous(breaks = seq(0, 5, 0.5), labels = seq(0, 5, 0.5),
                         minor_breaks = seq(0, 5, 0.25), expand = c(0, 0)) +
      labs(x = 'Simulation time (years)', y = 'Age (years)') +
      customTheme +
      theme(panel.grid = element_line(color = '#191919', size = 0.1),
            panel.grid.minor = element_line(color = '#191919', size = 0.1))
  )

  for (i in c(1:11)) {
    pplot <- suppressWarnings(
      pplot +
        geom_line(aes_string(x = (paste0('x', i)), y = 'y'), col = '#999999') +
        geom_point(aes_string(x = (paste0('x', i)), y = 'pmc_y'), fill = pmc_cols[1], size = pts, shape = pts_shp) +
        geom_point(aes_string(x = (paste0('x', i)), y = 'rtss_y'), fill = rtss_col, size = pts, shape = pts_shp)
    )
  }

  print(pplot)
  f_save_plot(pplot, plot_name = 'Fig1D', width = 4, height = 3, plot_dir = plot_dir)

  rm(dat, pplot)

} # (lexis_plot)
