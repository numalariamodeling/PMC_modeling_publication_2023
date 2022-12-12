source(file.path('analysis', '00_config.R'))

pmc_offset <- TRUE
pmc_single_supp <- TRUE

if (pmc_offset) {
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
    geom_ribbon(data = pdat1, aes(x = time, ymin = low_val, ymax = up_val),
                alpha = 0.3, fill = pmc_cols[1]) +
    geom_line(data = pdat1, aes(x = time, y = low_val, linetype = ' None'),
              size = 1.1, col = pmc_cols[1]) +
    geom_ribbon(data = pdat2, aes(x = time, ymin = low_val, ymax = up_val),
                alpha = 0.2, fill = pmc_cols[1]) +
    geom_line(data = pdat2, aes(x = time, y = low_val, linetype = '10 days'),
              size = 1.1, col = pmc_cols[1]) +
    geom_ribbon(data = pdat3, aes(x = time, ymin = low_val, ymax = up_val),
                alpha = 0.2, fill = pmc_cols[1]) +
    geom_line(data = pdat3, aes(x = time, y = low_val, linetype = '17 days'),
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

} # pmc_offset


## Supp  -  PMC efficacy by EIR and reference
if (pmc_single_supp) {

  trial_PEs <- fread(file.path("reference_dat", "trial_settings.csv")) %>%
    mutate(n_rounds = nchar(gsub("[^0-9]+", "", ipti_touchpoints_mth)))

  esu_PE <- fread(file.path("reference_dat", "IPTi_effectiveness.csv")) %>%
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
  pdat$eir_grp <- factor(pdat$eir_grp,
                         levels = c('(-Inf,10]', '(10,60]', '(60, Inf]'),
                         labels = c('EIR <10', 'EIR 10-60', 'EIR 60-256'))
  table(pdat$eir_grp, pdat$Annual_EIR, exclude = NULL)

  pdat <- pdat %>%
    dplyr::group_by(week, name, pmc_mode, eir_grp) %>%
    dplyr::summarise(mean_val = mean(mean_val, na.rm = T),
                     low_val = min(low_val, na.rm = T),
                     up_val = max(up_val, na.rm = T)) %>%
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
  pdat <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
    filter(name == 'PE_clinical_incidence' &
             rtss_coverage == 0 &
             pmc_coverage == 1 &
             cm_coverage == 0.6) %>%
    rename_with(~gsub('ipti', 'pmc', .x)) %>%
    filter(age_group == 'U1')

  pdat$eir_grp <- cut(pdat$Annual_EIR, c(-Inf, 10, 60, Inf))
  pdat$eir_grp <- factor(pdat$eir_grp,
                         levels = c('(-Inf,10]', '(10,60]', '(60, Inf]'),
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
                                                         eir_grp == 'EIR 60-128' ~ ' EMOD_highEIR'),
           Outcome = 'clinical_malaria') %>%
    dplyr::select(colnames(ref_df)) %>%
    arrange(author_year)

  ref_df <- ref_df %>% bind_rows(pdat)
  ref_df <- ref_df %>% mutate(n_rounds = ifelse(is.na(n_rounds) | n_rounds == 0, 1, n_rounds))
  table(ref_df$n_rounds)

  p2 <- ggplot() +
    geom_hline(yintercept = 22) +
    geom_pointrange(data = subset(ref_df, !is.na(author_year)),
                    aes(x = author_year, y = study_PE_mean, ymin = study_PE_lo, ymax = study_PE_up,
                        col = as.factor(n_rounds), shape = as.factor(Drug))) +
    scale_y_continuous(lim = c(0, 100), breaks = seq(0, 100, 10)) +
    labs(x = '', y = 'Effect size', color = 'N PMC rounds', shape = 'Drug') +
    scale_color_brewer(palette = 'Dark2', labels = c('pooled', 3, 4, 5)) +
    coord_flip() +
    customTheme

  pplot <- plot_grid(p1, p2, ncol = 1, labels = c('A', 'B'))
  print(pplot)
  f_save_plot(pplot, plot_name = 'S1_Fig1_pmc_efficacy', width = 9, height = 6.5, plot_dir = plot_dir)

} # (pmc_single_supp)
