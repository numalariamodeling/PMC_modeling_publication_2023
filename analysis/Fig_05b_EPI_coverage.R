## Fig_05b_EPI_coverage.R
source(file.path('analysis', '_config.R'))
source(file.path('analysis', '_fig05_helper_functions.R'))
library('rdhs')

##---------------
## EPI to PMC downscaling
##---------------
dhs_vacc_df <- fread(file.path("data_files", "dhs_vacc_df.csv"))
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
  fwrite(file.path("data_files", 'dhs_epi_dat.csv'))


##---------------
## Adjust for additional doses
##---------------

pmc_mode <- '7tp'
dhs_cov_dat <- fread(file.path("data_files", 'dhs_epi_dat.csv')) %>%
  dplyr::select(State, vacc_age, vacc_cov, vacc_cov_adj) %>%
  group_by(vacc_age) %>%
  mutate(round = cur_group_id(),
         epi_based = 1)

additional_vacc_ages <- touchpoints[[pmc_mode]][!(touchpoints[[pmc_mode]] %in% unique(dhs_cov_dat$vacc_age))]
i <- 0
for (add_age in additional_vacc_ages) {
  i <- i + 1
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
         vacc_cov = ifelse(vacc_age < max(touchpoints[['3tp']]) & !is.na(round_fix),
                           ((vacc_cov[round_fix - 1] + vacc_cov[round_fix + 1]) / 2),
                           vacc_cov),
         vacc_cov_adj = ifelse(vacc_age < max(touchpoints[['3tp']]) & !is.na(round_fix),
                               ((vacc_cov_adj[round_fix - 1] + vacc_cov_adj[round_fix + 1]) / 2),
                               vacc_cov_adj)) %>%
  group_by(vacc_age) %>%
  mutate(
    mean_vacc_cov = mean(vacc_cov),
    mean_vacc_cov_adj = mean(vacc_cov_adj))

target_dat <- dhs_cov_dat %>%
  dplyr::select(vacc_age, mean_vacc_cov) %>%
  unique() %>%
  pivot_longer(cols = -vacc_age) %>%
  mutate(value = 0.8, name = 'target')

fwrite(dhs_cov_dat, file.path("data_files", 'dhs_epi_dat_extended.csv'))

##---------------
## Figure
##---------------
pplot <- fread(file.path("data_files", 'dhs_epi_dat_extended.csv')) %>%
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
  scale_x_continuous(lim = c(0, 20), breaks = seq(0, 24, 3)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_color_manual(values = c('#353535', '#12719e')) +
  scale_fill_manual(values = c('#353535', '#12719e')) +
  labs(x = 'Age at EPI touchpoint', y = 'Coverage', color = '', fill = '') +
  customTheme

print(pplot)
f_save_plot(pplot, plot_name = paste0(plot_dir, 'Fig05_F'), width = 5.5, height = 3, plot_dir)
