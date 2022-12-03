source(file.path('analysis', '00_config.R'))


exp_name_pmc <- 'generic_single_PMC_vaccSP_IIV'
exp_name_rtss <- 'generic_single_RTSS_vaccSP_IIV'
exp_name_rtss_eir <- 'generic_RTSScov_EIR_constant_vaccSP_IIV'

pmc_single = TRUE
pmc_single_supp = TRUE
rtss_single = TRUE
lexis_plot = TRUE
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

  ## Offset inlcuded in simulations
  t_offset_sim = 10
  pdat1 <- cases_df %>%
    filter(pmc_coverage == 1 & name == 'PE_clinical_incidence') %>%
    filter(age >= 76 - t_offset_sim & age <= 76 - t_offset_sim + (12 * 7)) %>%
    mutate(age = age - (76 - t_offset_sim),
           time = age / 7)

  ## Offset post none
  t_offset2 = 0
  pdat2 <- cases_df %>%
    filter(pmc_coverage == 1 & name == 'PE_clinical_incidence') %>%
    filter(age >= 70 & age <= 76 + t_offset2 + (12 * 7)) %>%
    mutate(age = age - (76 + t_offset2),
           time = age / 7)

  ## Offset post 0
  t_offset3 = 7
  pdat3 <- cases_df %>%
    filter(pmc_coverage == 1 & name == 'PE_clinical_incidence') %>%
    filter(age >= 70 & age <= 76 + t_offset3 + (12 * 7)) %>%
    mutate(age = age - (76 + t_offset3),
           time = age / 7)


  pplot <- ggplot() +
    geom_hline(yintercept = 0, size = 0.7) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_ribbon(data = pdat1, aes(x = time, ymin = min_val, ymax = max_val),
                alpha = 0.3, fill = pmc_cols[1]) +
    geom_line(data = pdat1, aes(x = time, y = mean_val, linetype = ' None'),
              size = 1.1, col = pmc_cols[1]) +
    geom_ribbon(data = pdat2, aes(x = time, ymin = min_val, ymax = max_val),
                alpha = 0.2, fill = pmc_cols[1]) +
    geom_line(data = pdat2, aes(x = time, y = mean_val, linetype = '10 days'),
              size = 1.1, col = pmc_cols[1]) +
    geom_ribbon(data = pdat3, aes(x = time, ymin = min_val, ymax = max_val),
                alpha = 0.2, fill = pmc_cols[1]) +
    geom_line(data = pdat3, aes(x = time, y = mean_val, linetype = '17 days'),
              size = 1.1, col = pmc_cols[1]) +
    geom_point(data = mc_efficacy_weekly, aes(x = week, y = rel_red)) +
    scale_y_continuous(lim = c(-0.2, 1), breaks = seq(-0.2, 1, 0.2), labels = percent) +
    scale_x_continuous(breaks = c(-2:12), labels = c(-2:12)) +
    scale_linetype_manual(values = c('solid', 'dashed', 'dotted')) +
    labs(title = '', linetype = 'Offset',
         y = 'Modeled incidence reduction',
         x = 'Time after single PMC dose (weeks)') +
    customTheme_nogrid +
    theme(legend.key.width = unit(1, "cm")) +
    geom_text(aes(x = 0.8, y = 0.24, label = '  7    10  days'), size = 3) +
    geom_segment(aes(x = 1.8, y = 0.20, xend = -0.7, yend = 0.20),
                 arrow = arrow(length = unit(0.3, "cm")))

  print(pplot)
  f_save_plot(pplot, plot_name = paste0('Fig1A_supp'), width = 6, height = 3.5, plot_dir = plot_dir)


  pplot <- ggplot(data = pdat1) +
    geom_hline(yintercept = 0, size = 0.7) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_ribbon(data = pdat, aes(x = time, ymin = min_val, ymax = max_val), alpha = 0.3, fill = pmc_cols[1]) +
    geom_line(data = pdat, aes(x = time, y = mean_val, group = 1), col = pmc_cols[1], size = 1.1) +
    geom_point(data = mc_efficacy_weekly, aes(x = week, y = rel_red)) +
    scale_y_continuous(lim = c(-0.2, 1), breaks = seq(-0.2, 1, 0.2), labels = percent) +
    scale_x_continuous(breaks = c(-2:12), labels = c(-2:12)) +
    scale_color_manual(values = pmc_cols[1]) +
    scale_fill_manual(values = pmc_cols[1]) +
    labs(title = '',
         y = 'Modeled incidence reduction',
         x = 'Time after single PMC dose (weeks)') +
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
  fdir <- file.path('projects', 'rtss_scenarios', 'generic_setting', 'reference_Kintampo',
                    'kintampo_trial_summary_3month.csv')

  reference_filepath = file.path(drive, fdir)
  ref_df = fread(reference_filepath) %>%
    dplyr::filter(comparison == 'boost') %>%
    dplyr::select('time.group', 'time.out', 'PE') %>%
    dplyr::mutate(time_since_3rd_dose = time.out / 12)


  cases_df <- fread(file.path(simout_dir, exp_name_rtss, 'simdat_aggr_month.csv')) %>%
    filter(name == 'PE_clinical_incidence' & rtss_coverage == 1)


  pplot <- cases_df %>%
    filter(age > 274) %>%
    mutate(time = (age - min(age)) / 365) %>%
    ggplot() +
    geom_hline(yintercept = 0, size = 0.7) +
    geom_ribbon(aes(x = time, ymin = min_val, ymax = max_val), alpha = 0.3, fill = rtss_col) +
    geom_line(aes(x = time, y = mean_val), col = rtss_col, size = 1.1) +
    scale_y_continuous(lim = c(-0.25, 1), breaks = seq(-0.25, 1, 0.25), labels = percent) +
    scale_x_continuous(breaks = c(0:5), labels = c(0:5)) +
    labs(title = '',
         y = 'Modeled incidence reduction',
         x = 'Time after third RTS,S dose (years) ') +
    customTheme


  pplot <- pplot +
    geom_point(data = ref_df, aes(x = time_since_3rd_dose, y = PE)) +
    labs(caption = 'Generic simulation, not matched to Kintampo study site')

  print(pplot)
  f_save_plot(pplot, plot_name = 'Fig1B', width = 4, height = 3, plot_dir = plot_dir)

  #### RTS,S by EIR
  cases_df <- fread(file.path(simout_dir, exp_name_rtss_eir, 'simdat_aggr_month.csv')) %>%
    filter(name %in% c('PE_clinical_incidence') &
             rtss_coverage == 1 &
             cm_coverage == 0.6)


  pplot <- cases_df %>%
    mutate(age_yr = age / 365,
           time = (age - 274) / 365) %>%
    ggplot() +
    geom_hline(yintercept = 0, size = 0.7) +
    geom_ribbon(aes(x = age_yr, ymin = min_val, ymax = max_val, group = Annual_EIR, fill = as.factor(Annual_EIR)),
                alpha = 0.3) +
    geom_line(aes(x = age_yr, y = mean_val, group = Annual_EIR, col = as.factor(Annual_EIR)), size = 1.1) +
    scale_y_continuous(lim = c(-0.10, 1), breaks = seq(-0.10, 1, 0.20), labels = seq(-0.10, 1, 0.20) * 100) +
    scale_x_continuous(breaks = c(0:5), labels = c(0:5)) +
    labs(title = '',
         y = 'Modeled incidence reduction (%)',
         x = 'Age (years) ',
         color = 'annual EIR', fill = 'annual EIR') +
    scale_color_manual(values = getPalette) +
    scale_fill_manual(values = getPalette) +
    customTheme_nogrid +
    guides(colour = guide_legend(reverse = T),
           fill = guide_legend(reverse = T))


  print(pplot)
  f_save_plot(pplot, plot_name = 'Fig1B_supp', width = 8, height = 4, plot_dir = plot_dir)


  rm(cases_df, pplot)
} # (rtss_single) {


##---------------------------------------
## C  -  Lexis plot
##---------------------------------------
if (lexis_plot) {

  x <- seq(as.Date('0-01-01'), as.Date('13-01-01'), 'weeks')
  y <- round(as.numeric(rev((max(x) - (x)) / 365)), 2)

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
  pplot <- ggplot(data = subset(dat, x <= min(dat$x) + years(5))) +
    geom_line(aes(x = x, y = y), col = '#999999') +
    geom_point(aes(x = x, y = pmc_y), fill = 'deepskyblue3', size = pts, shape = pts_shp) +
    geom_point(aes(x = x, y = rtss_y), fill = '#55B748', size = pts, shape = pts_shp) +
    scale_x_date(lim = c(as.Date('0-01-01'), as.Date('8-01-01')),
                 date_breaks = '1 year', date_labels = '%Y',
                 expand = c(0, 0), date_minor_breaks = '6 month') +
    scale_y_continuous(breaks = seq(0, 5, 0.5), labels = seq(0, 5, 0.5),
                       minor_breaks = seq(0, 5, 0.25), expand = c(0, 0)) +
    labs(x = 'Simulation time (years)', y = 'Age (years)') +
    customTheme +
    theme(panel.grid = element_line(color = '#191919', size = 0.1),
          panel.grid.minor = element_line(color = '#191919', size = 0.1))

  for (i in c(1:11)) {
    pplot <- pplot +
      geom_line(aes_string(x = (paste0('x', i)), y = 'y'), col = '#999999') +
      geom_point(aes_string(x = (paste0('x', i)), y = 'pmc_y'), fill = 'deepskyblue3', size = pts, shape = pts_shp) +
      geom_point(aes_string(x = (paste0('x', i)), y = 'rtss_y'), fill = '#55B748', size = pts, shape = pts_shp)
  }

  print(pplot)
  f_save_plot(pplot, plot_name = 'Fig1D', width = 4, height = 3, plot_dir = plot_dir)

  rm(dat, pplot)

} # (lexis_plot) {

##---------------------------------------
## Supp  -  PMC efficacy by EIR and reference
##---------------------------------------
if (pmc_single_supp) {


  trial_PEs <- fread(file.path(pmc_path, "data", "trial_settings.csv")) %>%
    mutate(n_rounds = nchar(gsub("[^0-9]+", "", ipti_touchpoints_mth)))

  esu_PE <- fread(file.path(pmc_path, "data", "IPTi_effectiveness.csv")) %>%
    mutate(country_abbr = 'AFR') %>%
    dplyr::filter(Author == 'Esu', year == 2021, Drug %in% c("total", "SP") & Outcome == "clinical malaria") %>%
    dplyr::mutate(author_year = paste(country_abbr, year, Author, sep = "_"),
                  study_PE_mean = (1 - effectsize) * 100,
                  study_PE_up = (1 - low_ci) * 100,
                  study_PE_lo = (1 - up_ci) * 100,
                  n_rounds = NA) %>%
    dplyr::select(author_year, n_rounds, Drug, Outcome, study_PE_mean, study_PE_up, study_PE_lo)

  ref_df <- trial_PEs %>%
    filter(!is.na(study_PE_mean)) %>%
    rename(Drug = ipti_drug, Outcome = outcome) %>%
    mutate(author_year = paste(country_abbr, year, author, sep = "_")) %>%
    select(author_year, n_rounds, Drug, Outcome, study_PE_mean, study_PE_up, study_PE_lo) %>%
    bind_rows(esu_PE) %>%
    mutate(Drug = ifelse(author_year == 'Aponte_2009', 'total', Drug))

  days_offset <- 7
  pdat <- fread(file.path(simout_dir, exp_name_pmc, 'simdat_aggr_week.csv')) %>%
    filter(name == 'PE_clinical_incidence' &
             rtss_coverage == 0 &
             pmc_coverage == 1 &
             cm_coverage == 0.6) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    filter(age > 68 & age <= 76 + 70) %>%  #  mutate(time_res = '2-15month')
    mutate(week = (age - 76 - days_offset) / 7)

  pdat$eir_grp <- cut(pdat$Annual_EIR, c(-Inf, 10, 60, Inf))
  pdat$eir_grp <- factor(pdat$eir_grp, levels = c('(-Inf,10]', '(10,60]', '(60, Inf]'),
                         labels = c('EIR <10', 'EIR 10-60', 'EIR 60-256'))
  table(pdat$eir_grp, pdat$Annual_EIR, exclude = NULL)

  pdat <- pdat %>%
    dplyr::group_by(week, name, pmc_mode, eir_grp) %>%
    dplyr::summarise(mean_val = mean(mean_val, na.rm = T),
                     low_val = min(low_val, na.rm = T), #min(value, na.rm = T),
                     up_val = max(up_val, na.rm = T)) %>% # max(value, na.rm = T)) %>%
    ungroup() %>%
    mutate(study_PE_mean = mean_val * 100, study_PE_up = low_val * 100, study_PE_lo = up_val * 100,
           Drug = 'EMOD PMC-SP', author_year = case_when(eir_grp == 'EIR <10' ~ 'EMOD_lowEIR',
                                                         eir_grp == 'EIR 10-60' ~ 'EMOD_moderateEIR',
                                                         eir_grp == 'EIR 60-256' ~ 'EMOD_highEIR'),
           Outcome = 'clinical_malaria')


  p1 <- ggplot(data = pdat, aes(x = week, y = mean_val)) +
    geom_hline(yintercept = 0) +
    geom_ribbon(aes(ymin = low_val, ymax = up_val, fill = eir_grp), alpha = 0.3) +
    geom_line(aes(col = eir_grp, linetype = 'simulation')) +
    scale_y_continuous(lim = c(-1, 1)) +
    scale_color_manual(values = colorRampPalette(brewer.pal(9, "YlOrRd"))(9)[c(3, 5, 8)]) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(9, "YlOrRd"))(9)[c(3, 5, 8)]) +
    geom_point(data = mc_efficacy_weekly, aes(x = week, y = rel_red, shape = 'reference')) +
    scale_x_continuous(breaks = seq(-2, 12, 1), labels = seq(-2, 12, 1)) +
    labs(x = 'Time since single dose (weeks)', y = 'efficacy against clinical malaria',
         caption = 'Data estimates from "PMC_effect_PM3.docx" Malaria Consortium, 20220',
         linetype = '', shape = '') +
    customTheme +
    guides(colour = "none", fill = "none") +
    facet_wrap(~eir_grp)


  exp_name <- 'generic_PMCcov_EIR_PMC3_constant_vaccSP_IIV'
  pdat <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>% #simdat_aggr_month.csv
    filter(name == 'PE_clinical_incidence' &
             rtss_coverage == 0 &
             pmc_coverage == 1 &
             cm_coverage == 0.6) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    filter(age_group == 'U1')

  pdat$eir_grp <- cut(pdat$Annual_EIR, c(-Inf, 10, 60, Inf))
  pdat$eir_grp <- factor(pdat$eir_grp, levels = c('(-Inf,10]', '(10,60]', '(60, Inf]'),
                         labels = c('EIR <10', 'EIR 10-60', 'EIR 60-128'))
  table(pdat$eir_grp, pdat$Annual_EIR, exclude = NULL)


  pdat <- pdat %>%
    dplyr::group_by(name, pmc_mode, eir_grp) %>%
    dplyr::summarise(
      low_val = min(low_val, na.rm = T),
      up_val = max(up_val, na.rm = T),
      mean_val = mean(mean_val, na.rm = T)) %>%
    ungroup() %>%
    mutate(n_rounds = 3,
           study_PE_mean = mean_val * 100, study_PE_up = low_val * 100, study_PE_lo = up_val * 100,
           Drug = 'EMOD PMC-SP', author_year = case_when(eir_grp == 'EIR <10' ~ ' EMOD_lowEIR',
                                                         eir_grp == 'EIR 10-60' ~ ' EMOD_moderateEIR',
                                                         eir_grp == 'EIR 60-128' ~ ' EMOD_highEIR'), Outcome = 'clinical_malaria') %>%
    dplyr::select(colnames(ref_df)) %>%
    arrange(author_year)

  ref_df <- ref_df %>% bind_rows(pdat)
  ref_df <- ref_df %>% mutate(n_rounds = ifelse(is.na(n_rounds) | n_rounds == 0, 1, n_rounds))
  table(ref_df$n_rounds)

  p2 <- ggplot() +
    geom_hline(yintercept = 22) +
    geom_pointrange(data = subset(ref_df, !is.na(author_year)), aes(x = author_year, y = study_PE_mean,
                                                                    ymin = study_PE_lo, ymax = study_PE_up, col = as.factor(n_rounds), shape = as.factor(Drug))) +
    scale_y_continuous(lim = c(0, 100), breaks = seq(0, 100, 10)) +
    labs(x = '', y = 'Effect size', color = 'N PMC rounds', shape = 'Drug') +
    scale_color_brewer(palette = 'Dark2', labels = c('pooled', 3, 4, 5)) +
    coord_flip() +
    customTheme


  pplot <- plot_grid(p1, p2, ncol = 1, labels = c('A', 'B'))
  print(pplot)
  f_save_plot(pplot, plot_name = 'S1_Fig1_pmc_efficacy', width = 9, height = 6.5, plot_dir = plot_dir)

} # (pmc_single_supp) {
