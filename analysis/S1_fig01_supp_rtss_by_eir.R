#### RTS,S by EIR
source(file.path('analysis', '_config.R'))

exp_name <- 'generic_RTSScov_EIR_constant_vaccSP_IIV'

cases_df <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_month.csv')) %>%
  filter(name %in% c('PE_clinical_incidence') &
           rtss_coverage == 1 &
           cm_coverage == 0.6)

pplot <- cases_df %>%
  mutate(age_yr = age / 365,
         time = (age - 274) / 365) %>%
  ggplot() +
  geom_hline(yintercept = 0, size = 0.7) +
  geom_ribbon(aes(x = age_yr, ymin = low_val, ymax = up_val,
                  group = Annual_EIR, fill = as.factor(Annual_EIR)), alpha = 0.3) +
  geom_line(aes(x = age_yr, y = mean_val, group = Annual_EIR, col = as.factor(Annual_EIR)), size = 1.1) +
  scale_y_continuous(lim = c(-0.10, 1), breaks = seq(-0.10, 1, 0.20), labels = seq(-0.10, 1, 0.20) * 100) +
  scale_x_continuous(breaks = c(0:5), labels = c(0:5)) +
  labs(title = '',
       y = 'Modeled incidence reduction (%)',
       x = 'Age (years) ',
       color = 'annual EIR',
       fill = 'annual EIR') +
  scale_color_manual(values = getPalette) +
  scale_fill_manual(values = getPalette) +
  customTheme_nogrid +
  guides(colour = guide_legend(reverse = T),
         fill = guide_legend(reverse = T))

print(pplot)
f_save_plot(pplot, plot_name = 'Fig1B_supp', width = 8, height = 4, plot_dir = plot_dir)
