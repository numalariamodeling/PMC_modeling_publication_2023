##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## SI_population_demography.R
##---------------------


source(file.path('analysis', '_config.R'))
source(file.path('analysis', '_fig05_helper_functions.R'))

getPalette <- rev(colorRampPalette(brewer.pal(8, "RdYlBu"))(12))
getPalette2 <- rev(colorRampPalette(brewer.pal(8, "YlGnBu"))(12))
getPalette3 <- rev(colorRampPalette(brewer.pal(8, "BuPu"))(12))

prop_popU2_totalpop <- 0.068
require(readxl)

drive <- file.path(Sys.getenv('USERPROFILE'), 'NU Malaria Modeling Dropbox')
ipti_path <- file.path(drive, 'projects', 'ipti_pmc')
pmc_path <- ipti_path

nga_pop <- read_xlsx(file.path(pmc_path, 'NGA', 'population_size', 'GEOPODE_NGA_population_ages0to100_version2022.xlsx'), skip = 6, trim_ws = T,
                     .name_repair = "universal") %>%
  dplyr::select_if(~!(all(is.na(.)) | all(. == ""))) %>%
  filter(!(is.na(LGA))) %>%
  rename_with(~gsub('[...]', '', .x)) %>%
  filter(State %in% c(admin1_sp$ADM1_NAME, 'Akwa Ibom', 'Akwa lbom')) %>%
  group_by(State) %>%
  summarise(total_pop = sum(Total7)) %>%
  mutate(pop_U2 = round(total_pop * prop_popU2_totalpop, 0),
         ADM1_NAME = gsub('Akwa Ibom', 'Akwa lbom', State),
         State = gsub('Akwa lbom', 'Akwa Ibom', State)
  )

unique(nga_pop$State)


pop_map = F
if (pop_map) {

  proj_str <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  admin0_sp <- getShp(ISO = "NGA", admin_level = "admin0")
  admin1_all_sp <- getShp(ISO = "NGA", admin_level = "admin1")
  admin1_sp <- shapefile(file.path(drive, "data/nigeria/nigeria_shapefiles/nga_ipti/NGA_States_ipti_n20.shp"))
  admin1_sp_disolv <- shapefile(file.path(drive, "data/nigeria/nigeria_shapefiles/nga_ipti/nga_adm1_dissolved_ipti.shp"))
  admin1_all_sp.f <- as.MAPshp(admin1_all_sp)
  admin1_sp.f <- as.MAPshp(admin1_sp)
  admin1_sp_disolv.f <- as.MAPshp(admin1_sp_disolv)
  admin0_sp.f <- as.MAPshp(admin0_sp)

  pplot0 <- admin1_sp.f %>%
    left_join(nga_pop, by = "ADM1_NAME") %>%
    ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group, fill = total_pop), color = "white", size = 0.35) +  #microscopy_2018
    geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = NA, size = 0.35) +
    geom_polygon(data = admin1_sp_disolv.f, aes(x = long, y = lat, group = group), color = "black", fill = NA, size = 0.35) +
    theme(legend.position = "right") +
    labs(fill = 'Total population\n(all ages)') +
    scale_fill_viridis_c(labels = comma) +
    customTheme +
    theme_map()

  pplot1 <- admin1_sp.f %>%
    left_join(nga_pop, by = "ADM1_NAME") %>%
    ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group, fill = pop_U2), color = "white", size = 0.35) +  #microscopy_2018
    geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = NA, size = 0.35) +
    geom_polygon(data = admin1_sp_disolv.f, aes(x = long, y = lat, group = group), color = "black", fill = NA, size = 0.35) +
    theme(legend.position = "right") +
    labs(fill = 'Population U2\n(6.8% of total)') +
    scale_fill_viridis_c(labels = comma) +
    customTheme +
    theme_map()


  pplot2 <- nga_pop %>%
    #pivot_longer(cols = -c('State', 'ADM1_NAME')) %>%
    ggplot() +
    geom_col(aes(x = reorder(State, -total_pop), y = total_pop / 1000, col = 'total', fill = total_pop), position = position_dodge()) + #+ facet_wrap(~name, scales='free_y')
    geom_col(aes(x = reorder(State, -total_pop), y = total_pop * prop_popU2_totalpop / 1000, col = 'U2', fill = total_pop), position = position_dodge(), fill = NA) +
    scale_y_continuous(labels = comma) +
    scale_color_manual(values = c('darkgrey', 'darkred')) +
    labs(color = 'Population', x = 'State (ordered by population)', y = expr('Population 10'^3)) +
    scale_fill_viridis_c(labels = comma) +
    customTheme_nogrid +
    guides(fill = "none") +
    theme(legend.position = c(0.8, 0.9)) +
    coord_flip()


  pplot <- plot_grid(pplot0, pplot1, ncol = 1, labels = c('A', 'B'))
  pplot <- plot_grid(pplot, pplot2, rel_widths = c(1, 0.8), labels = c('', 'C'))
  f_save_plot(pplot, plot_name = paste0('SI_demography'), width = 10, height = 6, plot_dir = plot_dir)


}


### Treatment coverage

dat <- fread(file.path(pmc_path, 'simulation_inputs/nigeria/projection_csvs/CM', 'HS_by_LGA_v5.csv')) %>%
  mutate(NAME_1 = State) %>%
  filter(year == 2018) %>%
  group_by(NAME_1) %>%
  summarize(U5_coverage = mean(U5_coverage)) %>%
  add_nga_admin() %>%
  filter(southern_nga == 'Southern States') %>%
  dplyr::select(NAME_1, U5_coverage) %>%
  mutate(State = NAME_1, NAME_1 = ifelse(NAME_1 == 'Akwa Ibom', 'Akwa lbom', NAME_1)) %>%
  arrange(State) %>%
  mutate(ADM1_NAME = State)


summary(dat$U5_coverage)

pplot0 <- admin1_sp.f %>%
  left_join(dat, by = "ADM1_NAME") %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group, fill = U5_coverage), color = "white", size = 0.35) +  #microscopy_2018
  geom_polygon(aes(x = long, y = lat, group = group), color = "white", fill = NA, size = 0.35) +
  geom_polygon(data = admin1_sp_disolv.f, aes(x = long, y = lat, group = group), color = "black", fill = NA, size = 0.35) +
  theme(legend.position = "right") +
  labs(fill = 'Clinical treatment\ncoverage (U5)') +
  scale_fill_viridis_c(labels = comma, lim = c(0.2, 0.5), breaks = seq(0.2, 0.5, 0.1)) +
  customTheme +
  theme_map()


pplot2 <- ggplot(data = dat) +
  geom_col(aes(x = reorder(State, -U5_coverage), y = U5_coverage, fill = U5_coverage)) +
  labs(color = 'Population', x = 'State (ordered by treatment coverage)', y = 'Effective clinical\ntreatment coverage U5') +
  geom_hline(yintercept = 0.3193) +
  customTheme_nogrid +
  scale_fill_viridis_c(labels = comma, lim = c(0.2, 0.5), breaks = seq(0.2, 0.5, 0.1)) +
  theme(legend.position = 'None') +
  coord_flip()

pplot <- plot_grid(NULL, pplot0, NULL, ncol = 1, rel_heights = c(0.4, 1, 0.4))
pplot <- plot_grid(pplot, pplot2, ncol = 2, rel_widths = c(1, 0.8), labels = c('A', 'B'))
f_save_plot(pplot, plot_name = paste0('SI_CM'), width = 10, height = 6, plot_dir = plot_dir)
