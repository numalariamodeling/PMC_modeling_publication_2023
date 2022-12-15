source(file.path('analysis', '_config.R'))
source(file.path('_fig05_helper_functions.R'))
library('rdhs')

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
  mutate(vacc_age = as.numeric(vacc_age / (365 / 12))) %>%
  group_by(vacc_age) %>%
  summarize(mean_vacc_cov = mean(vacc_cov),
            mean_vacc_cov_adj = mean(vacc_cov_adj),
            min_vacc_cov = min(vacc_cov),
            min_vacc_cov_adj = min(vacc_cov_adj),
            max_vacc_cov = max(vacc_cov),
            max_vacc_cov_adj = max(vacc_cov_adj)) %>%
  ggplot() +
  geom_hline(yintercept = 0.8, col = '#db2b27') +
  geom_errorbar(aes(x = vacc_age, y = mean_vacc_cov, ymin = min_vacc_cov, ymax = max_vacc_cov, col = 'EPI'), width = 0) +
  geom_errorbar(aes(x = vacc_age + 0.2, y = mean_vacc_cov_adj, ymin = min_vacc_cov_adj, ymax = max_vacc_cov_adj, col = 'PMC'), width = 0) +
  geom_point(aes(x = vacc_age, y = mean_vacc_cov, fill = 'EPI'), size = 1.7, shape = 21) +
  geom_point(aes(x = vacc_age + 0.2, y = mean_vacc_cov_adj, fill = 'PMC'), size = 1.7, shape = 21) +
  scale_x_continuous(lim=c(0, 20),breaks = seq(0, 24, 3)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_color_manual(values = c('#353535', '#12719e')) +
  scale_fill_manual(values = c('#353535', '#12719e')) +
  labs(x = 'Age at EPI touchpoint', y = 'Coverage', color = '', fill = '') +
  customTheme
f_save_plot(pplot, plot_name = paste0('fig_5_coverage'), width = 5.5, height = 3, plot_dir)





cov <- seq(0, 1, 0.01)
rounds <- c(1:3)
dat <- as.data.frame(expand.grid(cov, rounds))
colnames(dat) <- c('vacc_cov', 'round')
sclalefactors <- EPI_to_pmc_cov(get_scaling_factors = T)

pplot <- dat %>%
  mutate(vacc_cov_adj = case_when(round == 1 ~ vacc_cov * sclalefactors[1],
                                  round == 2 ~ vacc_cov * sclalefactors[2],
                                  round == 3 ~ vacc_cov * sclalefactors[3]),
         vacc_dose = case_when(round == 1 ~ 'dtp-2', round == 2 ~ 'dtp-3', round == 3 ~ 'measles')) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_line(aes(x = vacc_cov, y = vacc_cov_adj, col = vacc_dose)) +
  scale_x_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  labs(title = '', x = 'EPI coverage (%)', y = 'PMC coverage (%)',
       color = 'EPI touchpoint',
       caption = paste0('Scaling factors per dose (1-3): ',
                        paste0(round(sclalefactors, 3), collapse = ', '))) +
  scale_color_brewer(palette = "Dark2")

f_save_plot(pplot, 'EPI_to_PMC_scaling', file.path(pmc_path, 'simulation_inputs'), width = 5, height = 3.5)
