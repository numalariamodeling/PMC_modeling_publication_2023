##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## Fig_A1.2.4_agegroups.R
##---------------------
source(file.path('analysis', '_config.R'))

exp_name <- 'generic_PMC_RTSS_EIR_vaccSP_IIV'

dat_aggr <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
  filter(Annual_EIR < eir_max) %>%
  mutate(pmc_mode_fct = ifelse(pmc_rtss_cov == '0-0', 'None',
                               ifelse(pmc_rtss_cov == '0-0.8', 'RTS,S',
                                      ifelse(pmc_rtss_cov == '0.8-0', 'PMC-3', 'PMC-3 + RTS,S'))))

dat_aggr$pmc_mode_fct <- factor(dat_aggr$pmc_mode_fct,
                                levels = rev(c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')),
                                labels = rev(c('None', 'PMC-3', 'RTS,S', 'PMC-3 + RTS,S')))

age_groups <- paste0("U", c(1, 2, 5, 10))
dat_aggr$age_group <- factor(dat_aggr$age_group, levels = age_groups, labels = age_groups)

table(dat_aggr$pmc_mode_fct, dat_aggr$pmc_mode)

pdat <- dat_aggr %>% filter(pmc_mode == '3tp' &
                              pmc_mode_fct != 'None' &
                              name %in% c('clinical_cases_averted'))

pmc_col_alpha <- '#4EA3D1'
rtss_col_alpha <- '#FBB040'

pplot1 <- ggplot(data = pdat,
                 aes(x = Annual_EIR, y = median_val)) +
  geom_hline(yintercept = 0, alpha = 0) +
  geom_ribbon(aes(ymin = 0, ymax = median_val, fill = as.factor(pmc_mode_fct)), alpha = 0.3, show.legend = F) +
  geom_line(aes(col = as.factor(pmc_mode_fct))) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~age_group, nrow = 1) +
  scale_color_manual(values = c(pmc_cols[1], rtss_col, rtss_pmc_cols[1])) +
  scale_fill_manual(values = c(pmc_col_alpha, rtss_col_alpha, rtss_pmc_cols[1])) +
  scale_shape_manual(values = c(21, 22)) +
  labs(subtitle = '',
       x = 'annual EIR (log scale)',
       y = 'Clinical cases averted\nper 1000 population',
       color = '', fill = '') +
  customTheme_nogrid +
  guides(panel.spacing = unit(1.3, "lines")) +
  theme(panel.grid = element_blank()) +
  scale_x_log10(breaks = unique(pdat$Annual_EIR), labels = unique(pdat$Annual_EIR))

pdat <- dat_aggr %>% filter(pmc_mode_fct != 'None' & name %in% c('severe_cases_averted'))
pplot2 <- ggplot(data = pdat,
                 aes(x = Annual_EIR, y = median_val)) +
  geom_hline(yintercept = 0, alpha = 0) +
  geom_line(aes(col = as.factor(pmc_mode_fct))) +
  geom_ribbon(aes(ymin = 0, ymax = median_val, fill = as.factor(pmc_mode_fct)), alpha = 0.3, show.legend = F) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~age_group, nrow = 1) +
  scale_color_manual(values = c(pmc_cols[1], rtss_col, rtss_pmc_cols[1])) +
  scale_fill_manual(values = c(pmc_col_alpha, rtss_col_alpha, rtss_pmc_cols[1])) +
  scale_shape_manual(values = c(21, 22)) +
  labs(subtitle = '',
       x = 'annual EIR (log scale)',
       y = 'Severe cases averted\nper 1000 population',
       color = '', fill = '') +
  customTheme_nogrid +
  guides(colour = guide_legend(reverse = T),
         fill = guide_legend(reverse = T),
         panel.spacing = unit(1.3, "lines")) +
  theme(panel.grid = element_blank()) +
  scale_x_log10(breaks = unique(pdat$Annual_EIR), labels = unique(pdat$Annual_EIR))

pplot <- plot_combine(list(pplot1, pplot2))

print(pplot)
f_save_plot(pplot, 'Fig A1.2.4',
            file.path(plot_dir), width = 12, height = 6, units = 'in', device_format = device_format)

