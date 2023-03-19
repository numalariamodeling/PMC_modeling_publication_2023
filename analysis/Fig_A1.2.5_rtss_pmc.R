##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## Fig_A1.2.5_rtss_pmc.R
##---------------------

exp_name <- 'generic_PMCmode_RTSS_vaccSP_IIV'

dat_aggr <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_week.csv'))
dat_aggr$pmc_mode_fct <- factor(dat_aggr$pmc_mode,
                                levels = c('3tp', '5tp', '7tp2ndyr'),
                                labels = c('PMC-3', 'PMC-5', 'PMC-7'))

dat_aggr %>%
  group_by(Annual_EIR, pmc_rtss_cov) %>%
  unique() %>%
  tally() %>%
  pivot_wider(names_from = pmc_rtss_cov, values_from = n)

age_y <- 2
pdat1 <- subset(dat_aggr, age <= y * 365 &
  name == 'clinical_cases' &
  pmc_coverage == 0.8 &
  pmc_mode %in% c('3tp', '5tp', '7tp2ndyr'))

pdat2 <- dat_aggr %>%
  filter(rtss_coverage == 0 &
           pmc_coverage == 0 &
           age <= age_y * 365 &
           name == 'clinical_cases') %>%
  dplyr::select(-pmc_mode, -pmc_mode_fct, -rtss_coverage, -pmc_coverage)

pdat3 <- subset(dat_aggr, age <= age_y * 365 &
  name == 'clinical_cases' &
  rtss_coverage != 0 &
  pmc_coverage == 0) %>%
  dplyr::select(-pmc_mode, -pmc_mode_fct)


pplot <- ggplot(data = pdat1) +
  geom_vline(xintercept = c(365)) +
  geom_line(data = pdat2,
            aes(x = age, y = median_val), col = 'darkgrey') +
  geom_ribbon(data = pdat2,
              aes(x = age, ymin = low_val, ymax = up_val),
              alpha = 0.4, fill = 'darkgrey') +
  geom_ribbon(aes(x = age, ymin = low_val, ymax = up_val,
                  group = interaction(pmc_mode_fct, rtss_coverage),
                  fill = pmc_mode_fct), alpha = 0.4) +
  geom_ribbon(data = pdat3,
              aes(x = age, ymin = low_val, ymax = up_val), alpha = 0.3, fill = rtss_col) +
  geom_line(aes(x = age, y = median_val, group = interaction(pmc_mode_fct, rtss_coverage),
                col = pmc_mode_fct, linetype = as.factor(rtss_coverage)), size = 1) +
  geom_line(data = pdat3,
            aes(x = age, y = median_val),
            col = rtss_col, size = 1) +
  facet_wrap(~pmc_mode_fct, nrow = 1, scales = 'free_x') +
  scale_linetype_manual(values = c('solid', 'dashed')) +
  scale_color_manual(values = c(pmc_cols[c(1, 2, 5)], rtss_col)) +
  scale_fill_manual(values = c(pmc_cols[c(1, 2, 5)], rtss_col)) +
  scale_x_continuous(lim = c(0, age_y * 365),
                     breaks = seq(0, age_y * 365, 30 * 6),
                     labels = seq(0, age_y * 12, 6) / 12) +
  labs(col = '', fill = '', x = 'Age (years)',
       y = 'Annual cases \nper 1000 population') +
  theme(panel.grid.major = element_blank(),
        legend.position = 'None') +
  customTheme_nogrid

print(pplot)
f_save_plot(pplot, paste0('Fig A1.2.5'),
            file.path(plot_dir), width = 8, height = 3.5, units = 'in', device_format = device_format)


