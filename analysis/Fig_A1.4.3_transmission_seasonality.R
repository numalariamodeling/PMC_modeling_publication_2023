##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## Fig_A1.4.3_transmission_seasonality.R
##---------------------

baseline_pfpr <- fread(file.path(simout_dir, exp_name, 'baseline_pfpr.csv'))
pfpr_df <- fread(file.path("data_files", "ndhs_2018_pfpr_eir_df.csv")) %>%
  mutate(State = ifelse(State == 'Akwa Ibom', 'Akwa lbom', State),
         ADM1_NAME = State)

proj_str <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
admin1_sp <- shapefile(file.path("data_files/nga_shp/NGA_States_ipti_n20.shp"))
admin1_sp.f <- as.MAPshp(admin1_sp)

pplot_pfpr_eir <- ggplot(data = subset(baseline_pfpr, statistic == 'mean_val')) +
  geom_smooth(data = baseline_pfpr_sub, aes(x = Annual_EIR, y = pfpr_U5), col = 'dodgerblue3') +
  geom_point(data = pfpr_df, aes(x = predEIR, y = RDT_2018, group = State), fill = 'darkorange', col = 'darkgrey', size = 1.7, shape = 21) +
  scale_x_log10(breaks = unique(baseline_pfpr$Annual_EIR)) +
  scale_y_continuous(lim = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.2)) +
  labs(x = 'annual EIR (log scale)', y = expr(italic(Pf) * 'PR U5 (RDT)')) +
  customTheme

pplot_pfprmap <- admin1_sp.f %>%
  left_join(pfpr_df, by = "ADM1_NAME") %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group, fill = RDT_2018 * 100), color = "white", size = 0.35) +  #microscopy_2018
  geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = NA, size = 0.35) +
  theme(legend.position = "right") +
  labs(fill = 'PfPR U5\n(RDT, \nNDHS 2018)') %>%
    scale_fill_gradientn(colours = rev(RdYlBupalette), lim = c(0, 60), labels = comma) +
  customTheme +
  theme_map()

pplot <- plot_grid(pplot_pfpr_eir, pplot_pfprmap, ncol = 2, labels = c('A', 'B'))
print(pplot)
f_save_plot(pplot, paste0('Fig A1.4.3AB'), file.path(plot_dir), width = 10, height = 3.5, units = 'in', device_format = device_format)


### Note Fig A1.4.3 A and B generated in ??
source(file.path('analysis', '_config.R'))

dat <- fread(file.path('simulation_inputs/Seasonality/SouthernNGA_extractedMonthlyEIR.csv'))
dat2 <- fread(file.path('simulation_inputs/Seasonality/seasonality_eir_NGApmc_multipliers.csv'))

dat <- dat %>%
  pivot_longer(cols = -c(month))

pplot <- ggplot(data = dat) +
  geom_hline(yintercept = 0, alpha = 0) +
  geom_line(aes(x = month, y = value, group = name), col = 'dodgerblue3') +
  scale_x_continuous(breaks = c(1:12)) +
  labs(x = 'Month', y = 'Monthly EIR multiplier\n(simulation estimates)', col = 'State') +
  f_getCustomTheme()

pplot <- plot_grid(pplot, labels = c('C'))
f_save_plot(pplot, paste0('Fig A1.4.3C'),
            file.path(plot_dir), width = 7, height = 3.5, units = 'in', device_format = device_format)