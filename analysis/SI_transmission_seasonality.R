##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## SI_transmission_seasonality.R
##---------------------

source(file.path('analysis', '_config.R'))


dat <- fread(file.path('scenario_csvs/Seasonality/SouthernNGA_extractedMonthlyEIR.csv'))
dat2 <- fread(file.path('scenario_csvs/Seasonality/seasonality_eir_NGApmc_multipliers.csv'))

dat <- dat %>%
  pivot_longer(cols = -c(month))

pplot <- ggplot(data = dat) +
  geom_hline(yintercept = 0, alpha = 0) +
  geom_line(aes(x = month, y = value, group = name), col = 'dodgerblue3') +
  scale_x_continuous(breaks = c(1:12)) +
  labs(x = 'Month', y = 'Monthly EIR multiplier\n(simulation estimates)', col = 'State') +
  f_getCustomTheme()

pplot <- plot_grid(pplot, labels=c('C'))
f_save_plot(pplot, paste0('SI_fig_seasonality'),
            file.path(plot_dir), width = 7, height =3.5, units = 'in', device_format = device_format)