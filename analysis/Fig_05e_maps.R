source(file.path('analysis', '_config.R'))
source(file.path('analysis', '_fig05_helper_functions.R'))

getPalette <- rev(colorRampPalette(brewer.pal(8, "RdYlBu"))(12))
getPalette2 <- rev(colorRampPalette(brewer.pal(8, "YlGnBu"))(12))
getPalette3 <- rev(colorRampPalette(brewer.pal(8, "BuPu"))(12))

proj_str <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
admin0_sp <- getShp(ISO = "NGA", admin_level = "admin0")
admin1_all_sp <- getShp(ISO = "NGA", admin_level = "admin1")
admin1_sp <- shapefile(file.path(drive, "data/nigeria/nigeria_shapefiles/nga_ipti/NGA_States_ipti_n20.shp"))
admin1_sp_disolv <- shapefile(file.path(drive, "data/nigeria/nigeria_shapefiles/nga_ipti/nga_adm1_dissolved_ipti.shp"))
admin1_all_sp.f <- as.MAPshp(admin1_all_sp)
admin1_sp.f <- as.MAPshp(admin1_sp)
admin1_sp_disolv.f <- as.MAPshp(admin1_sp_disolv)
admin0_sp.f <- as.MAPshp(admin0_sp)

prop_popU2_totalpop <- 0.068
require(readxl)
nga_pop <- read_xlsx(file.path(pmc_path, 'NGA', 'population_size', 'GEOPODE_NGA_population_ages0to100_version2022.xlsx'), skip = 6, trim_ws = T,
                     .name_repair = "universal") %>%
  dplyr::select_if(~!(all(is.na(.)) | all(. == ""))) %>%
  filter(!(is.na(LGA))) %>%
  rename_with(~gsub('[...]', '', .x)) %>%
  filter(State %in% admin1_sp$ADM1_NAME) %>%
  group_by(State) %>%
  summarise(total_pop = sum(Total7)) %>%
  mutate(pop_U2 = round(total_pop * prop_popU2_totalpop, 0),
         State = gsub('Akwa lbom', 'Akwa Ibom', State))
unique(nga_pop$State)


df_out <- fread(file.path(simout_dir, 'NGA_simdat_aggr_agegroup.csv')) %>%
  mutate(State = seasonality) %>%
  filter(statistic == 'mean_val' &
           age_group == 'U2') %>%
  left_join(nga_pop) %>%
  mutate(clinicalinc_ppa = clinical_cases / 1000,
         clinicalcases_pop = clinical_cases / 1000 * pop_U2,
         severeinc_ppa = severe_cases / 1000,
         severecases_pop = severe_cases / 1000 * pop_U2,
         clinical_cases_averted_pop = clinical_cases_averted / 1000 * pop_U2, ,
         severe_cases_averted_pop = severe_cases_averted / 1000 * pop_U2,) %>%
  mutate(ADM1_NAME = gsub('Akwa Ibom', 'Akwa lbom', State))


df_out$scen <- factor(df_out$scen,
                      levels = c('counterfactual', 'PMC-3', 'PMC-5', 'PMC-7', 'RTS,S', 'PMC-3 + RTS,S'),
                      labels = c('counterfactual', 'PMC-3', 'PMC-5', 'PMC-7', 'RTS,S', 'PMC-3 + RTS,S'))

######-------------------------
generate_maps = TRUE
dhs_data_maps = TRUE
sim_result_maps = TRUE
scenarios = TRUE
if (generate_maps) {

  if (dhs_data_maps) {

    dhs_pfpr_dat <- fread(file.path('dhs_pfpr_df.csv'))
    dhs_epicov_dat <- fread(file.path('dhs_epi_dat.csv'))
    dhs_epicov_dat_aggr <- dhs_epicov_dat %>%
      dplyr::select(State, ADM1_NAME, mean_vacc_cov, mean_vacc_cov_adj) %>%
      unique()


    pplot0 <- admin1_sp.f %>%
      left_join(dhs_pfpr_dat, by = "ADM1_NAME") %>%
      ggplot() +
      geom_polygon(data = admin1_sp_disolv.f, aes(x = long, y = lat, group = group), fill = 'deepskyblue2', color = "black", size = 0.35) +
      geom_polygon(data = admin1_sp.f, aes(x = long, y = lat, group = group), color = "white", fill = NA, size = 0.35) +
      #geom_polygon(data = admin1_all_sp.f, aes(x = long, y = lat, group = group), color = "grey", fill = NA, size = 0.35) +
      geom_polygon(data = admin0_sp.f, aes(x = long, y = lat, group = group), color = "black", fill = NA, size = 0.35) +
      theme(legend.position = "right") +
      customTheme +
      theme_map()

    pplot_pfpr <- admin1_sp.f %>%
      left_join(dhs_pfpr_dat, by = "ADM1_NAME") %>%
      ggplot() +
      geom_polygon(aes(x = long, y = lat, group = group, fill = RDT_2018 * 100), color = "white", size = 0.35) +  #microscopy_2018
      geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = NA, size = 0.35) +
      geom_polygon(data = admin1_sp_disolv.f, aes(x = long, y = lat, group = group), color = "black", fill = NA, size = 0.35) +
      theme(legend.position = "right") +
      labs(fill = 'PfPR U5\n(RDT, \nNDHS 2018)') %>%
        scale_fill_gradientn(colours = rev(RdYlBupalette), lim = c(0, 60), labels = comma) +
      customTheme +
      theme_map()

    pplot_epicov_perdose <- admin1_sp.f %>%
      left_join(dhs_epicov_dat, by = "ADM1_NAME") %>%
      ggplot() +
      geom_polygon(aes(x = long, y = lat, group = group, fill = vacc_cov * 100), color = "white", size = 0.35) +
      geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = NA, size = 0.35) +
      geom_polygon(data = admin1_sp_disolv.f, aes(x = long, y = lat, group = group), color = "black", fill = NA, size = 0.35) +
      theme(legend.position = "right") +
      facet_wrap(~vacc_round) +
      labs(fill = 'Coverage') %>%
        scale_fill_gradientn(colours = rev(getPalette2), limits = c(40, 100)) +
      customTheme +
      theme_map()

    pplot_epicov <- admin1_sp.f %>%
      left_join(subset(dhs_epicov_dat, vacc_round == 'dpt_2'), by = "ADM1_NAME") %>%
      ggplot() +
      geom_polygon(aes(x = long, y = lat, group = group, fill = vacc_cov * 100), color = "white", size = 0.35) +
      geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = NA, size = 0.35) +
      geom_polygon(data = admin1_sp_disolv.f, aes(x = long, y = lat, group = group), color = "black", fill = NA, size = 0.35) +
      theme(legend.position = "right") +
      labs(fill = 'DTP-2 coverage') %>%
        scale_fill_gradientn(colours = rev(getPalette2), limits = c(40, 100)) +
      customTheme +
      theme_map()

    pplot_epicov_adj <- admin1_sp.f %>%
      left_join(subset(dhs_epicov_dat, vacc_round == 'dpt_2'), by = "ADM1_NAME") %>%
      ggplot() +
      geom_polygon(aes(x = long, y = lat, group = group, fill = vacc_cov_adj * 100), color = "white", size = 0.35) +
      geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = NA, size = 0.35) +
      geom_polygon(data = admin1_sp_disolv.f, aes(x = long, y = lat, group = group), color = "black", fill = NA, size = 0.35) +
      theme(legend.position = "right") +
      labs(fill = 'PMC-3\n1st dose coverage') %>%
        scale_fill_gradientn(colours = rev(getPalette2), limits = c(40, 100)) +
      customTheme +
      theme_map()

    print(pplot0)
    print(pplot_pfpr)
    print(pplot_epicov_perdose)
    print(pplot_epicov)
    print(pplot_epicov_adj)

    f_save_plot(pplot0, paste0('pplot0'), file.path(plot_dir), width = 4, height = 3, units = 'in', device_format = device_format)
    f_save_plot(pplot_pfpr, paste0('pplot_pfpr'), file.path(plot_dir), width = 4, height = 3, units = 'in', device_format = device_format)
    f_save_plot(pplot_epicov_perdose, paste0('pplot_epicov_perdose'), file.path(plot_dir), height = 4, width = 15, units = 'in', device_format = device_format)
    f_save_plot(pplot_epicov, paste0('pplot_epicov'), file.path(plot_dir), width = 4, height = 3, units = 'in', device_format = device_format)
    f_save_plot(pplot_epicov_adj, paste0('pplot_epicov_adj'), file.path(plot_dir), width = 4, height = 3, units = 'in', device_format = device_format)

  } #dhs_data_maps

  ######-------------------------

  if (sim_result_maps) {

    summary(df_out$clinicalcases_pop)

    plot_dat <- df_out %>% filter(scen == 'counterfactual' & coveragemode == 'operational')

    pplot_clinicalinc_U2 <- admin1_sp.f %>%
      left_join(plot_dat, by = "ADM1_NAME") %>%
      ggplot() +
      geom_polygon(aes(x = long, y = lat, group = group, fill = clinical_cases), color = "white", size = 0.35) +
      geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = NA, size = 0.35) +
      geom_polygon(data = admin1_sp_disolv.f, aes(x = long, y = lat, group = group), color = "black", fill = NA, size = 0.35) +
      theme(legend.position = "right") +
      labs(fill = 'clinicalinc_U2') +
      scale_fill_gradientn(colours = rev(RdYlBupalette), lim = c(0, 3000), labels = comma) +
      customTheme +
      theme_map()

    pplot_clinicalcases_popU2 <- admin1_sp.f %>%
      left_join(plot_dat, by = "ADM1_NAME") %>%
      ggplot() +
      geom_polygon(aes(x = long, y = lat, group = group, fill = clinicalcases_pop), color = "white", size = 0.35) +
      geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = NA, size = 0.35) +
      geom_polygon(data = admin1_sp_disolv.f, aes(x = long, y = lat, group = group), color = "black", fill = NA, size = 0.35) +
      theme(legend.position = "right") +
      labs(fill = 'Clinical cases\n per total population U2') +
      scale_fill_gradientn(colours = rev(RdYlBupalette), limits = c(0, 890000), labels = comma) +
      customTheme +
      theme_map()

    plot_dat <- df_out %>% filter(scen == 'PMC-3' & coveragemode == 'operational')
    summary(plot_dat$PE_clinical_cases) * 100
    summary(plot_dat$PE_severe_cases) * 100

    summary(plot_dat$clinical_cases_averted)
    summary(plot_dat$severe_cases_averted)

    summary(plot_dat$clinical_cases_averted_pop)
    summary(plot_dat$severe_cases_averted_pop)


    pplot_PE_clinical_incidence_U2 <- admin1_sp.f %>%
      left_join(plot_dat, by = "ADM1_NAME") %>%
      ggplot() +
      geom_polygon(aes(x = long, y = lat, group = group, fill = PE_clinical_cases * 100), color = "white", size = 0.35) +
      geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = NA, size = 0.35) +
      geom_polygon(data = admin1_sp_disolv.f, aes(x = long, y = lat, group = group), color = "black", fill = NA, size = 0.35) +
      theme(legend.position = "right") +
      labs(fill = 'PMC-3\n% reduction\nclinical cases U2') +
      scale_fill_gradientn(colours = rev(getPalette3), limits = c(0, 10), labels = comma) +
      customTheme +
      theme_map()
  }

  f_save_plot(pplot0, paste0('pplot0'), file.path(plot_dir), width = 4, height = 3, units = 'in', device_format = device_format)
  f_save_plot(pplot_pfpr, paste0('pplot_pfpr'), file.path(plot_dir), width = 4, height = 3, units = 'in', device_format = device_format)
  f_save_plot(pplot_epicov_perdose, paste0('pplot_epicov_perdose'), file.path(plot_dir), height = 4, width = 15, units = 'in', device_format = device_format)
  f_save_plot(pplot_epicov, paste0('pplot_epicov'), file.path(plot_dir), width = 4, height = 3, units = 'in', device_format = device_format)
  f_save_plot(pplot_epicov_adj, paste0('pplot_epicov_adj'), file.path(plot_dir), width = 4, height = 3, units = 'in', device_format = device_format)

  f_save_plot(pplot_clinicalinc_U2, paste0('pplot_clinicalinc_U2'), file.path(plot_dir), width = 4, height = 3, units = 'in', device_format = device_format)
  f_save_plot(pplot_clinicalcases_popU2, paste0('pplot_clinicalcases_popU2'), file.path(plot_dir), width = 4, height = 3, units = 'in', device_format = device_format)
  f_save_plot(pplot_PE_clinical_incidence_U2, paste0('pplot_PE_clinical_incidence_U2'), file.path(plot_dir), width = 4, height = 3, units = 'in', device_format = device_format)

  pplot <- plot_grid(pplot_pfpr, pplot_clinicalinc_U2, pplot_epicov_adj, pplot_PE_clinical_incidence_U2, ncol = 2, labels = c('A', 'B', 'C', 'D'))
  f_save_plot(pplot, paste0('pmap_combined'), file.path(plot_dir), width = 12, height = 6, units = 'in', device_format = device_format)


  if (scenarios) {
    plot_dat <- df_out %>% filter(scen != 'counterfactual')
    summary(plot_dat$clinical_cases_averted)

    pplot_clinical_cases_averted <- admin1_sp.f %>%
      left_join(plot_dat, by = "ADM1_NAME") %>%
      ggplot() +
      geom_polygon(aes(x = long, y = lat, group = group, fill = clinical_cases_averted), color = "white", size = 0.35) +
      geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = NA, size = 0.35) +
      geom_polygon(data = admin1_sp_disolv.f, aes(x = long, y = lat, group = group), color = "black", fill = NA, size = 0.35) +
      facet_grid(scen ~ coveragemode) +
      theme(legend.position = "right") +
      labs(fill = 'Clinical cases averted\nper 1000 population U2') +
      scale_fill_gradientn(colours = rev(getPalette3), limits = c(0, 1000)) +
      customTheme +
      theme_map()

    pplot_severe_cases_averted <- admin1_sp.f %>%
      left_join(plot_dat, by = "ADM1_NAME") %>%
      ggplot() +
      geom_polygon(aes(x = long, y = lat, group = group, fill = severe_cases_averted), color = "white", size = 0.35) +
      geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = NA, size = 0.35) +
      geom_polygon(data = admin1_sp_disolv.f, aes(x = long, y = lat, group = group), color = "black", fill = NA, size = 0.35) +
      facet_grid(scen ~ coveragemode) +
      theme(legend.position = "right") +
      labs(fill = 'Severe cases averted\nper 1000 population U2') +
      scale_fill_gradientn(colours = rev(getPalette3), limits = c(0, 20)) +
      customTheme +
      theme_map()

    f_save_plot(pplot_clinical_cases_averted, paste0('NGA_map_clinical_incidence_averted'), file.path(plot_dir),
                width = 12, height = 10, units = 'in', device_format = device_format)
    f_save_plot(pplot_severe_cases_averted, paste0('NGA_map_severe_incidence_averted'), file.path(plot_dir),
                width = 12, height = 10, units = 'in', device_format = device_format)

    ###--------------------------------
    plot_dat <- df_out %>% filter(scen != 'counterfactual')
    summary(plot_dat$clinical_cases_averted_pop)
    summary(plot_dat$severe_cases_averted_pop)

    summary(log(plot_dat$clinical_cases_averted_pop))
    summary(log(plot_dat$severe_cases_averted_pop))

    pplot_clinical_cases_averted <- admin1_sp.f %>%
      left_join(plot_dat, by = "ADM1_NAME") %>%
      ggplot() +
      geom_polygon(aes(x = long, y = lat, group = group, fill = (clinical_cases_averted_pop)), color = "white", size = 0.35) +
      geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = NA, size = 0.35) +
      geom_polygon(data = admin1_sp_disolv.f, aes(x = long, y = lat, group = group), color = "black", fill = NA, size = 0.35) +
      facet_grid(scen ~ coveragemode) +
      theme(legend.position = "right") +
      labs(fill = 'Clinical cases averted\nper total population U2') +
      #scale_fill_viridis_c(option='D') +
      scale_fill_gradientn(colours = (YlGnpalette), trans = "log10", labels = comma) +
      customTheme +
      theme_map()

    pplot_severe_cases_averted <- admin1_sp.f %>%
      left_join(plot_dat, by = "ADM1_NAME") %>%
      ggplot() +
      geom_polygon(aes(x = long, y = lat, group = group, fill = severe_cases_averted_pop), color = "white", size = 0.35) +
      geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = NA, size = 0.35) +
      geom_polygon(data = admin1_sp_disolv.f, aes(x = long, y = lat, group = group), color = "black", fill = NA, size = 0.35) +
      facet_grid(scen ~ coveragemode) +
      theme(legend.position = "right") +
      labs(fill = 'Severe cases averted\nper total population U2') +
      scale_fill_gradientn(colours = (YlGnpalette), trans = "log10", labels = comma) +
      customTheme +
      theme_map()

    f_save_plot(pplot_clinical_cases_averted, paste0('NGA_map_clinical_cases_averted'), file.path(plot_dir),
                width = 10, height = 10, units = 'in', device_format = device_format)
    f_save_plot(pplot_severe_cases_averted, paste0('NGA_map_severe_cases_averted'), file.path(plot_dir),
                width = 10, height = 10, units = 'in', device_format = device_format)

  }
} #generate_maps


### for text
df_out <- fread(file.path(simout_dir, 'NGA_simdat_aggr_agegroup.csv')) %>%
  mutate(State = seasonality) %>%
  filter(statistic == 'mean_val' &
           age_group == 'U2') %>%
  left_join(nga_pop) %>%
  mutate(clinicalinc_ppa = clinical_cases / 1000,
         clinicalcases_pop = clinical_cases / 1000 * pop_U2,
         severeinc_ppa = severe_cases / 1000,
         severecases_pop = severe_cases / 1000 * pop_U2,
         clinical_cases_averted_pop = clinical_cases_averted / 1000 * pop_U2, ,
         severe_cases_averted_pop = severe_cases_averted / 1000 * pop_U2,) %>%
  mutate(ADM1_NAME = gsub('Akwa Ibom', 'Akwa lbom', State))


df_out$scen <- factor(df_out$scen,
                      levels = c('counterfactual', 'PMC-3', 'PMC-5', 'PMC-7', 'RTS,S', 'PMC-3 + RTS,S'),
                      labels = c('counterfactual', 'PMC-3', 'PMC-5', 'PMC-7', 'RTS,S', 'PMC-3 + RTS,S'))


plot_dat <- df_out %>% filter(scen == 'PMC-3')
tapply(plot_dat$PE_clinical_cases * 100, plot_dat$coveragemode, summary)
tapply(plot_dat$PE_severe_cases * 100, plot_dat$coveragemode, summary)

tapply(plot_dat$clinical_cases_averted, plot_dat$coveragemode, summary)
tapply(plot_dat$severe_cases_averted, plot_dat$coveragemode, summary)

tapply(plot_dat$clinical_cases_averted_pop, plot_dat$coveragemode, summary)
tapply(plot_dat$severe_cases_averted_pop, plot_dat$coveragemode, summary)

#### All scenarios
plot_dat <- df_out %>%
  group_by(coveragemode, scen) %>%
  filter(scen != 'counterfactual') %>%
  summarize(PE_clinical_cases = mean(PE_clinical_cases),
            PE_severe_cases = mean(PE_severe_cases),
            clinical_cases_averted = mean(clinical_cases_averted),
            severe_cases_averted = mean(severe_cases_averted)) %>%
  pivot_longer(cols = -c(coveragemode, scen)) %>%
  pivot_wider(names_from = coveragemode, values_from = value) %>%
  mutate(diff = target / operational)

tapply(plot_dat$diff, plot_dat$name, summary)

tapply(plot_dat$clinical_cases_averted, plot_dat$coveragemode, summary)
tapply(plot_dat$severe_cases_averted, plot_dat$coveragemode, summary)

tapply(plot_dat$clinical_cases_averted_pop, plot_dat$coveragemode, summary)
tapply(plot_dat$severe_cases_averted_pop, plot_dat$coveragemode, summary)
