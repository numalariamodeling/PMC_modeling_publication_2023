##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## Fig_A1.2.2_eir-pfpr-incidence_curves.R
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


##WHO, Global Malaria Programme. A framework for malaria elimination. Geneva, World Health Organization; 2017.
pfpr_WHOcuts <- c(0, 0.01, 0.10, 0.35)
inc_WHOcuts <- c(0, 100, 250, 450, 450)

pdat <- dat_aggr %>% filter(pmc_mode_fct == 'None' &
                              name %in% c('PfHRP2_Prevalence', 'clinical_cases_pppa', 'severe_cases_pppa'))

pdat$name <- factor(pdat$name, levels = c('PfHRP2_Prevalence', 'clinical_cases_pppa', 'severe_cases_pppa'),
                    labels = c('PfHRP2 Prevalence', 'Clinical cases ppa', 'Severe cases ppa'))

pplot <- ggplot() +
  scale_color_brewer(palette = 'Dark2') +
  labs(subtitle = '', x = 'annual EIR (log scale)', color = 'Age group', fill = 'Age group') +
  customTheme_nogrid +
  guides(colour = guide_legend(reverse = T),
         fill = guide_legend(reverse = T),
         panel.spacing = unit(1.3, "lines")) +
  theme(panel.grid = element_blank())

p1 <- pplot +
  geom_hline(yintercept = pfpr_WHOcuts) +
  #geom_rect(ymin=pfpr_WHOcuts[2] , ymax=pfpr_WHOcuts[3], xmin=-Inf, xmax=Inf, fill='orange') +
  geom_line(data = subset(pdat, name == 'PfHRP2 Prevalence'), aes(x = Annual_EIR, y = mean_val, col = age_group)) +
  scale_x_log10(breaks = unique(pdat$Annual_EIR), labels = unique(pdat$Annual_EIR)) +
  labs(y = expr(italic(Pf) * "PR")) +
  scale_y_continuous(lim = c(0, 1), expand = c(0, 0))

p2 <- pplot +
  geom_hline(yintercept = inc_WHOcuts / 1000) +
  geom_line(data = subset(pdat, name == 'Clinical cases ppa'), aes(x = Annual_EIR, y = mean_val, col = age_group)) +
  scale_x_log10(breaks = unique(pdat$Annual_EIR), labels = unique(pdat$Annual_EIR)) +
  labs(y = "Clinical incidence\nper person per year") +
  scale_y_continuous(lim = c(0, 6))

p3 <- pplot +
  geom_line(data = subset(pdat, name == 'Severe cases ppa'), aes(x = Annual_EIR, y = mean_val, col = age_group)) +
  scale_x_log10(breaks = unique(pdat$Annual_EIR), labels = unique(pdat$Annual_EIR)) +
  labs(y = "Severe incidence\nper person per year") +
  scale_y_continuous(lim = c(0, 0.12))


pplot <- plot_combine(list(p1, p2, p3), labels = c('A', 'B', 'C'), ncol = 3)
print(pplot)
f_save_plot(pplot, 'Fig S1.2.2', file.path(plot_dir),
            width = 12, height = 5, units = 'in', device_format = device_format)

fwrite(pdat, file.path('figures/csv', 'pdat_Fig S1.2.2.csv'))

