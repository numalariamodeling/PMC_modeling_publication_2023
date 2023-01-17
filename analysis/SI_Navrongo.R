##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## SI_Navrongo.R
##---------------------
source(file.path('analysis', '_config.R'))

study_start <- as.Date("2000-09-30")
study_end <- as.Date("2004-06-01")
study_years <- as.numeric((study_end - study_start) / 365)
ipti_touchpoints <- c(91, 122, 274, 365)
ipti_touchpoints_wks <- ceiling(ipti_touchpoints / 7)

##---------------------------------------
## Reference
##---------------------------------------
refs <- 'ipti3_mean'
ref_wk = 'week3'

exp_name <- 'navrongo_ipti_test_decayshape'
expsubdir <- '_navrongo'

simout_dir <- file.path(ipti_path, 'simulation_output', expsubdir, exp_name)
output_dir <- simout_dir


##---------------------------------------
## Custom functions
##---------------------------------------
f_load_eir_sweeps <- function(aggregate_years = TRUE, keep_runs = FALSE, sweeep_vars = NULL) {

  if (is.null(sweeep_vars))sweeep_vars <- c('EIR.scale.factor', 'cm_cov_U5', 'study_group')
  grp_channels <- c('date', 'month', sweeep_vars)
  summary_reports <- list.files(simout_dir, pattern = 'U1_PfPR_ClinicalIncidence')
  eir_sweep_list <- list()
  for (report in summary_reports) {
    eir_sweep_list[[length(eir_sweep_list) + 1]] <- read.csv(file.path(simout_dir, report)) %>%
      mutate(study_group = gsub("_", "", gsub(".csv", "", gsub("U1_PfPR_ClinicalIncidence", "", report))))
  }
  eir_sweep <- eir_sweep_list %>% bind_rows()

  if (!('month') %in% colnames(eir_sweep)) {
    ## if month not in dataframe, simulation from cohort setup, starting end of September
    age_cutoffs <- c(-Inf, 30, 60, 90, 121, 152, 183, 212, 243, 274, 305, 336, 370)
    eir_sweep$month <- cut(eir_sweep$age, breaks = age_cutoffs, labels = c(c(10:12), c(1:9)))
    eir_sweep$month <- as.numeric(eir_sweep$month)
    if (max(eir_sweep$year) == 1)eir_sweep$year <- 2001

    eir_sweep <- eir_sweep %>%
      dplyr::group_by_at(.vars = c('Run_Number', 'month', 'year', all_of(sweeep_vars))) %>%
      dplyr::summarize(Cases.U1 = sum(Cases.U1),
                       Pop.U1 = mean(Pop.U1)) %>%
      dplyr::mutate(total_cases_U1 = Cases.U1 * Pop.U1 * 1 / 365)

  }

  eir_sweep <- eir_sweep %>%
    dplyr::group_by_at(.vars = c('Run_Number', 'month', 'year', all_of(sweeep_vars))) %>%
    dplyr::mutate(date = as.Date(paste(year, month, '01', sep = '-'))) %>%
    dplyr::mutate(total_cases_U1 = Cases.U1 * Pop.U1 * 30 / 365)


  if (aggregate_years) {
    eir_sweep <- eir_sweep %>%
      dplyr::group_by_at(.vars = c('Run_Number', 'month', all_of(sweeep_vars))) %>%
      dplyr::summarize(PfPR.U1 = mean(PfPR.U1),
                       Cases.U1 = mean(Cases.U1),
                       total_cases_U1 = sum(total_cases_U1),
                       Pop.U1 = sum(Pop.U1))

  }


  if (!keep_runs) {
    grp_channels <- grp_channels[!(grp_channels == 'Run_Number')]
    eir_sweep <- eir_sweep %>%
      ungroup() %>%
      dplyr::group_by_at(.vars = grp_channels) %>%
      dplyr::summarize(PfPR.U1 = mean(PfPR.U1),
                       Cases.U1 = mean(Cases.U1),
                       total_cases_U1 = mean(total_cases_U1),
                       Pop.U1 = mean(Pop.U1))
  }

  eir_sweep <- eir_sweep %>%
    mutate(Cases.U1 = Cases.U1 / study_years,
           lPI = quantile(Cases.U1, probs = 0.025, na.rm = T) / study_years,
           hPI = quantile(Cases.U1, probs = 0.975, na.rm = T) / study_years,
           risk = total_cases_U1 / Pop.U1 / study_years,
           incidence = total_cases_U1 / Pop.U1 * 1000) %>%
    rename(cases_ppa_U1 = Cases.U1)

  return(eir_sweep)
}

f_load_ref_df <- function(fig_suffix = 'S8') {

  #from Cairns et al 2011, cohort starting in sep
  #ref = c(0.2203, 0.1868, 0.0909, 0.0747, 0.0792, 0.0518, 0.0605, 0.0146, 0.0107, 0.0150, 0.1160, 0.1795)
  #refm = c(0.1757, 0.1835, 0.1558, 0.1127, 0.0594, 0.0314, 0.0190, 0.0131, 0.0176, 0.0342, 0.0800, 0.1467)
  ref_df_S8 <- fread(file.path(ipti_path, 'data', 'Cairns_2011', 'figS8.csv'))
  colnames(ref_df_S8) <- tolower(colnames(ref_df_S8))
  ref_df_S8 <- ref_df_S8 %>%
    dplyr::mutate(date = as.Date(paste(year, month, '01', sep = '-'))) %>%
    dplyr::group_by(month, treatment_grp, type) %>%
    dplyr::summarise(malaria_risk = mean(malaria_risk))

  ## Malaria incidence from Chandramohan et al 2005, Figure 3, placebo only
  ref_df_fig3 <- fread(file.path(ipti_path, 'data', 'Cairns_2011', 'fig3.csv'))
  colnames(ref_df_fig3) <- tolower(colnames(ref_df_fig3))
  ref_df_fig3 <- ref_df_fig3 %>% mutate(date = as.Date(paste(year, month, '01', sep = '-')))


  #prevalence from Chandramohan et al 2005
  chandramohan_age <- c(273, 61, 4, 37, 63, 93, 122, 152, 184, 213, 243, 274, 307, 341)
  chandramohan_ref <- c(0.086, 0.528, 20, 25, 14, 7.9, 6, 6.3, 3.8, 2.8, 2.2, 1, 12, 23)
  chandramohan_df <- as.data.frame(cbind('age' = chandramohan_age, 'ref' = chandramohan_ref))

  if (fig_suffix == ' chandramohan')ref_df = chandramohan_df
  if (fig_suffix == 'S8')ref_df = ref_df_S8
  if (fig_suffix == 'fig3')ref_df = ref_df_fig3
  return(ref_df)
}

f_load_pd_sweeps_aggr <- function(keep_runs = FALSE, monthly = F) {
  ##load simulation outputs with EIR sweeps
  # Cases.U1 =clinical symptoms per person per year

  summary_reports <- list.files(simout_dir, pattern = 'Agebin_aggr_PfPR_ClinicalIncidence')
  pd_sweep_list <- list()
  for (report in summary_reports) {
    pd_sweep_list[[length(pd_sweep_list) + 1]] <- read.csv(file.path(simout_dir, report)) %>%
      mutate(treatment_grp = gsub("_", "", gsub(".csv", "", gsub("Agebin_aggr_PfPR_ClinicalIncidence", "", report))))
  }
  pd_sweep <- pd_sweep_list %>%
    bind_rows()

  if (!('C50_SP' %in% colnames(pd_sweep)))pd_sweep$C50_SP = -9
  if (!('kmax_SP' %in% colnames(pd_sweep)))pd_sweep$kmax_SP = -9
  ## Get total cases and total population
  pd_sweep <- pd_sweep %>%
    dplyr::filter(agebin <= 2) %>%
    dplyr::group_by(Run_Number, agebin, cm_cov_U5, treatment_grp, C50_SP, kmax_SP) %>%
    dplyr::mutate(total_cases = Cases * Pop / study_years,
                  total_population = Pop * study_years,
                  age_mth = floor(agebin * 365 / 30))


  #Aggregate runs (optional)
  if (!keep_runs) {
    grp_channels <- c('agebin', 'cm_cov_U5', 'treatment_grp', 'C50_SP', 'kmax_SP')
    pd_sweep <- pd_sweep %>%
      ungroup() %>%
      dplyr::group_by_at(.vars = grp_channels) %>%
      dplyr::summarize(PfPR = mean(PfPR),
                       Cases = mean(Cases),
                       total_cases = mean(total_cases),
                       total_population = mean(total_population),
                       Pop = mean(Pop))
  }

  if (monthly) {
    grp_channels <- c('age_mth', 'cm_cov_U5', 'treatment_grp', 'C50_SP', 'kmax_SP', 'Run_Number')
    pd_sweep <- pd_sweep %>%
      ungroup() %>%
      dplyr::group_by_at(.vars = grp_channels) %>%
      mutate(PfPR = mean(PfPR),
             Cases = sum(Cases),
             total_cases = mean(total_cases),
             total_population = mean(total_population),
             Pop = mean(Pop))

  }

  #if(!monthly){
  pd_sweep <- pd_sweep %>%
    mutate(age_days = agebin * 365,
           age_weeks = agebin * 365 / 7,
           #age_mth = floor(agebin * 365 / 30),
           Cases = Cases / study_years,
           risk = total_cases / total_population,
           incidence = total_cases / total_population * 1000) %>%
    rename(cases_ppa = Cases)
  # }


  return(pd_sweep)
}

f_agebin_protection <- function(cases_averted.df, start_wk = 9, end_wk = 60, ipti_set = 'week2', ref_set = 'ca_ref_vec', SAVE = F) {
  #separate out by treatment groups to create placebo comparison
  inc_trt <- as.data.frame(cases_averted.df[which(cases_averted.df$treatment_grp == 'treatment'),])
  inc_plc <- as.data.frame(cases_averted.df[which(cases_averted.df$treatment_grp == 'placebo'),])
  inc_plc$cases_ppa_plc <- inc_plc$cases_ppa

  #create smoothed placebo inc values across all sims
  plcdat <- group_by(inc_plc, agebin)
  plcdat <- as.data.frame(summarise(plcdat, avgcases_ppa_plc = mean(cases_ppa),
                                    plcincidence = mean(incidence)))
  inc_plc <- merge(inc_plc, plcdat, by = 'agebin')

  #recombine treatment and placebo (sim only and smoothed)
  inc_red.df <- merge(inc_trt, inc_plc[, c('age_weeks', 'cases_ppa_plc', 'avgcases_ppa_plc', 'Run_Number', 'plcincidence')], by = c('age_weeks', 'Run_Number'))
  inc_red.df$Inc_reduction_percent <- 0
  inc_red.df$Inc_reduction_percent <- ((inc_red.df$avgcases_ppa_plc - inc_red.df$cases_ppa) / inc_red.df$avgcases_ppa_plc) * 100

  #set up confidence intervals and set to have only 1 seed remain
  inc_redc.df <- inc_red.df %>%
    dplyr::group_by(agebin, age_weeks, age_days, cm_cov_U5, treatment_grp) %>%
    dplyr::summarise(lPI = quantile(Inc_reduction_percent, probs = 0.025, na.rm = T),
                     hPI = quantile(Inc_reduction_percent, probs = 0.975, na.rm = T),
                     meanred = mean(Inc_reduction_percent, na.rm = T),
                     meaninc = mean(incidence, na.rm = T),
                     meanincp = mean(plcincidence, na.rm = T))
  # dplyr::mutate(lCI=min(Inc_reduction_percent),
  #               hCI=max(Inc_reduction_percent),
  #               meanred=mean(Inc_reduction_percent))

  #inc_conf<-inc_redc.df[which(inc_redc.df$Run_Number==1),]
  #inc_conf<-inc_redc.df
  return(inc_redc.df)
}


##---------------------------------------
### Chandramohan et al 2005, Fig 3 malaria incidence over time in placebo group
##---------------------------------------
ref_df <- f_load_ref_df(fig_suffix = 'fig3')
ref_df_pfpr <- fread(file.path(ipti_path, 'data', 'Chandramohan_2005', 'pfpr.csv')) %>%
  mutate(date = as.Date(date, format = '%m/%d/%Y'), month = month(date)) %>%
  select(-date)
sim_df <- f_load_eir_sweeps(aggregate_years = FALSE) %>% left_join(ref_df_pfpr)

pplot_incidence <- ggplot(data = subset(sim_df, study_group == 'placebo')) +
  geom_line(aes(x = date, y = cases_ppa_U1, col = as.factor(round(EIR.scale.factor, 3)))) +
  geom_point(data = subset(ref_df, age_grp == 'infants'), aes(x = date, y = malaria_incidence_ppa, shape = age_grp)) +
  ylim(0, 4) +
  labs(x = '', y = 'clinical cases pppa', color = 'simulation\n(eir scale factor)', shape = 'data') +
  scale_color_manual(values = 'deepskyblue3') +
  scale_x_date(date_breaks = '3 months', date_labels = '%b\n%Y') +
  customTheme_nogrid

pplot_pfpr <- ggplot(data = subset(sim_df, study_group == 'placebo')) +
  geom_line(aes(x = date, y = PfPR.U1, col = as.factor(round(EIR.scale.factor, 3)))) +
  geom_pointrange(aes(x = date, y = pfpr, ymin = pfpr_low, ymax = pfpr_up)) +
  scale_color_manual(values = 'deepskyblue3') +
  labs(x = '', y = expr(italic(Pf) * 'PR')) +
  theme(legend.position = 'None') +
  scale_x_date(date_breaks = '3 months', date_labels = '%b\n%Y') +
  customTheme_nogrid

pplot <- plot_combine(list(pplot_incidence, pplot_pfpr), labels = c('A', 'B'))
print(pplot)

f_save_plot(pplot, paste0('fig_SI_navrongo_clinicalcases_EIRsweep'),
            file.path(plot_dir), width = 12, height = 6, units = 'in', device_format = device_format)


##---------------------------------------
### navrongo_ipti_test_decayshape_adj_eff95
##---------------------------------------
exp_name <- 'navrongo_ipti_test_decayshape_adj_eff95'
simout_dir <- file.path(ipti_path, 'simulation_output/_navrongo', exp_name)
output_dir <- simout_dir

ca_ref.df <- fread(file.path(ipti_path, 'data', 'Cairns_2008', 'table1.csv')) %>%
  rename(weekipti = week_since_ipti) %>%
  as.data.frame()
ca_ref.df$week2 <- seq(ipti_touchpoints_wks[2], ipti_touchpoints_wks[2] + 11, 1)
ca_ref.df$week3 <- seq(ipti_touchpoints_wks[3], ipti_touchpoints_wks[3] + 11, 1)
ca_ref.df$week4 <- seq(ipti_touchpoints_wks[4], ipti_touchpoints_wks[4] + 11, 1)

ca_ref.df_long <- ca_ref.df %>%
  pivot_longer(cols = -c('weekipti', 'week2', 'week3', 'week4')) %>%
  separate(name, into = c('ipti_dose', 'statistic'), sep = '_') %>%
  mutate(ipti_dose = ifelse(ipti_dose == 'iptiall', 'combined', gsub("ipti", "PMC-", ipti_dose))) %>%
  pivot_wider(names_from = statistic, values_from = value)

cases_averted.df <- f_load_pd_sweeps_aggr(keep_runs = T)

conf.df <- f_agebin_protection(cases_averted.df) %>%
  ungroup() %>%
  dplyr::select(age_days, meanred, lPI, hPI)

conf.df$week_after_dose <- (conf.df$age_days - ipti_touchpoints[2]) / 7
conf.df_2 <- conf.df %>%
  mutate(ipti_dose = 'PMC-2') %>%
  dplyr::select(-age_days)

conf.df$week_after_dose <- (conf.df$age_days - ipti_touchpoints[3]) / 7
conf.df_3 <- conf.df %>%
  mutate(ipti_dose = 'PMC-3') %>%
  dplyr::select(-age_days)

conf.df$week_after_dose <- (conf.df$age_days - ipti_touchpoints[4]) / 7
conf.df_4 <- conf.df %>%
  mutate(ipti_dose = 'PMC-4') %>%
  dplyr::select(-age_days)

conf.df_perdose <- rbind(conf.df_2, conf.df_3, conf.df_4) %>%
  filter(week_after_dose > -1) %>%
  mutate(week_after_dose = round(week_after_dose, 0))
summary(conf.df_perdose$week_after_dose)

conf.df_combined <- conf.df_perdose %>%
  group_by(week_after_dose) %>%
  summarize(meanred = mean(meanred),
            lPI = mean(lPI),
            hPI = mean(hPI)) %>%
  mutate(ipti_dose = 'combined')
conf.df_perdose <- rbind(conf.df_perdose, conf.df_combined)

ipti_levels <- c('PMC-2', 'PMC-3', 'PMC-4', 'combined')
conf.df_perdose$ipti_dose <- factor(conf.df_perdose$ipti_dose, levels = ipti_levels, labels = ipti_levels)
ca_ref.df_long$ipti_dose <- factor(ca_ref.df_long$ipti_dose, levels = ipti_levels, labels = ipti_levels)

pplot1 <- ggplot(data = subset(ca_ref.df_long, ipti_dose != 'combined')) +
  geom_ribbon(data = subset(conf.df_perdose, ipti_dose != 'combined'), aes(x = week_after_dose, y = meanred, ymin = lPI, ymax = hPI), fill = 'deepskyblue3', alpha = 0.4) +
  geom_line(data = subset(conf.df_perdose, ipti_dose != 'combined'), aes(x = week_after_dose, y = meanred), col = 'deepskyblue3') + #col = 'deepskyblue3') +
  geom_hline(yintercept = 0) +
  #ylim(-200, 100) +
  scale_x_continuous(breaks = seq(0, 12, 1), lim = c(-1, 12)) +
  labs(x = "Time after dose (weeks)", y = "% reduction clinical incidence") +
  geom_point(aes(x = weekipti, y = mean), col = 'black') +
  geom_pointrange(aes(x = weekipti, y = mean, ymin = cilow, ymax = ciup), col = 'black', size = 0.1) +
  facet_wrap(~ipti_dose, ncol = 3) +
  customTheme_nogrid +
  theme(legend.position = 'None')


pplot2 <- ggplot(data = subset(ca_ref.df_long, ipti_dose == 'combined')) +
  geom_ribbon(data = subset(conf.df_perdose, ipti_dose == 'combined'), aes(x = week_after_dose, y = meanred, ymin = lPI, ymax = hPI), fill = 'deepskyblue3', alpha = 0.4) +
  geom_line(data = subset(conf.df_perdose, ipti_dose == 'combined'), aes(x = week_after_dose, y = meanred), col = 'deepskyblue3') + #col = 'deepskyblue3') +
  geom_hline(yintercept = 0) +
  #ylim(-200, 100) +
  scale_x_continuous(breaks = seq(0, 12, 1), lim = c(-1, 12)) +
  labs(x = "Time after dose (weeks)", y = "% reduction clinical incidence") +
  geom_point(aes(x = weekipti, y = mean), col = 'black') +
  geom_pointrange(aes(x = weekipti, y = mean, ymin = cilow, ymax = ciup), col = 'black', size = 0.1) +
  facet_wrap(~ipti_dose, ncol = 3) +
  customTheme_nogrid +
  theme(legend.position = 'None')


## overall impact
PE_ref_df <- as.data.frame(NA) %>% mutate(PE_mean = 0.243,
                                          PE_low = 0.218,
                                          PE_up = 0.300,
                                          source = 'data')

pplot3 <- cases_averted.df %>%
  filter(agebin <= 1.0) %>%
  group_by(Run_Number, EIR.scale.factor, cm_cov_U5, cm_cov_adults, ipti_cov, ipti_touchpoints, treatment_grp) %>%
  summarize(total_cases = sum(total_cases),
            Pop = mean(Pop)) %>%
  pivot_wider(names_from = treatment_grp, values_from = c(total_cases, Pop)) %>%
  mutate(PE = 1 - (total_cases_treatment / total_cases_placebo)) %>%
  group_by(EIR.scale.factor, cm_cov_U5, cm_cov_adults, ipti_cov, ipti_touchpoints) %>%
  summarize(PE_mean = mean(PE),
            PE_min = min(PE),
            PE_max = max(PE),
            PE_low = quantile(PE, probs = 0.05, na.rm = T),
            PE_up = quantile(PE, probs = 0.95, na.rm = T)) %>%
  mutate(source = 'sim', Uage = 'children U1') %>%
  bind_rows(PE_ref_df) %>%
  ggplot() +
  geom_hline(yintercept = c(0, 0.5)) +
  geom_pointrange(aes(x = source, y = PE_mean, ymin = PE_low, ymax = PE_up, col = EIR.scale.factor)) +
  labs(subtitle = 'children U1', x = '', y = '% reduction in clinical cases') +
  scale_y_continuous(breaks = seq(0, 0.5, 0.1), labels = seq(0, 0.5, 0.1) * 100, expand = c(0, 0)) +
  customTheme +
  theme(legend.position = 'None', panel.grid.major.x = element_blank())

pplot <- plot_grid(pplot2, pplot3, labels = c('B', 'C'), align = 'hv', rel_widths = c(1, 0.4))
pplot <- plot_grid(pplot1, pplot, nrow = 2, labels = c('A', NA, NA))
print(pplot)

f_save_plot(pplot, paste0('fig_SI_navrongo_PMC'),
            file.path(plot_dir), width = 12, height = 6, units = 'in', device_format = device_format)




