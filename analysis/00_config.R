pckg <- c("dplyr","data.table", "wrapr", "dplyr", "tidyr", "lubridate", "zoo", "ggthemes", "cowplot", "scales", "ggplot2", "RColorBrewer", "see",'pracma')
a <- lapply(pckg, require, character.only = TRUE)
rm(a)

source(file.path('analysis','00_plotter_helpers.R'))
source(file.path('analysis','00_load_generic_simdat.R'))
drive <- file.path(Sys.getenv('USERPROFILE'), 'NU Malaria Modeling Dropbox')

ipti_path <- file.path(drive,'projects', 'ipti_pmc')
pmc_path <- ipti_path
plot_dir <- file.path(ipti_path, 'project_notes/publication/1_pmc_rtss_generic/raw_figures/')
data_path <- file.path(ipti_path, 'data')
if (!dir.exists(plot_dir))dir.create(plot_dir)
simout_dir <- file.path(ipti_path, 'simulation_output', '_pmc_rtss_generic')

###----------------------------------------------------
SAVE <- TRUE
eir_max <- 200

###-------------------  Define custom objects
max_years <- c(1, 2, 5, 10)
age_groups <- paste0('U', max_years)
age_label_values <- seq(0, 10, 2) # seq(0, 8, length.out = 5)
age_labels <- age_label_values # paste0(age_label_values, '-', (age_label_values + 1))
default_eir <- 32
default_seasons <- c('season1')  # c('comstant','season1','season2','season3')
default_Uage <- 'U1'
ipti_rtss_levels <- c('None', 'IPTi', 'RTS,S', 'IPTi + RTS,S')


###-------------------  Define custom themes and colors
theme_set(theme_bw())
customTheme <- f_getCustomTheme() + guides(panel.spacing = unit(1.3, "lines"))
customTheme_nogrid <- customTheme + theme(panel.grid = element_blank())

device_format <- c('pdf', 'png')

#https://coolors.co/palette/03045e-023e8a-0077b6-0096c7-00b4d8-48cae4-90e0ef-ade8f4-bdedf6
rtss_col <- '#FBB040'
pmc_cols <- c('#4EA3D1', '#007BBD', '#CBDB47', '#99CA3C', '#208B3A')
rtss_pmc_cols <- c('#fa8e9e', '#d97373')
getPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))(9)[3:9]
custom_cols <- c("#4EA3D1", "#007BBD", "#51571C", "#99CA3C", "#208B3A", rtss_col, '#d97373')
CustomPalette_sub <- c('#55b748', '#00BFFF', '#FFA500') # if only 1 IPTi mode
Dark2palette <- colorRampPalette(brewer.pal(8, "Dark2"))(12)
RdBupalette <- colorRampPalette(brewer.pal(8, "RdBu"))(12)
RdYlBupalette <- colorRampPalette(brewer.pal(8, "RdYlBu"))(12)
BrBgpalette <- colorRampPalette(brewer.pal(8, "BrBG"))(12)
YlGnpalette <- colorRampPalette(brewer.pal(9, "YlGnBu"))(9)
ORpalette <- colorRampPalette(brewer.pal(9, "Oranges"))(9)
REpalette <- colorRampPalette(brewer.pal(9, "Reds"))(9)




###-------------------  Generic
format_num = function(x, d = 0) {
  x = format(round(x, d), format = "d", big.mark = ",")
  return(x)
}

dat_table <- function(dat, grp_vars) {
  dat %>%
    group_by_at(.vars = grp_vars) %>%
    unique() %>%
    tally() %>%
    pivot_wider(names_from = grp_vars[length(grp_vars)], values_from = n) %>%
    kbl() %>%
    kable_styling(bootstrap_options = tbl_opts, full_width = F, position = "left", fixed_thead = T)
}


gen_pmc_mode_fct <- function(dat) {
  dat$pmc_mode_fct <- factor(dat$pmc_mode,
                              levels = c('3tp', '4tp', '4tp2ndyr', '5tp', '5tp2ndyr', '6tp2ndyr', '7tp2ndyr', 'RTS,S'),
                              labels = c('PMC-3', 'PMC-4', 'PMC-4', 'PMC-5', 'PMC-5', 'PMC-6', 'PMC-7', 'RTS,S'))

  dat$pmc_mode_fct2 <- factor(dat$pmc_mode,
                               levels = c('3tp', '4tp', '4tp2ndyr', '5tp', '5tp2ndyr', '6tp2ndyr', '7tp2ndyr', 'RTS,S'),
                               labels = c('PMC-3', 'PMC-4a', 'PMC-4b', 'PMC-5a', 'PMC-5b', 'PMC-6', 'PMC-7', 'RTS,S'))

  dat$pmc_mode_ab <- factor(dat$pmc_mode,
                             levels = c('3tp', '4tp', '4tp2ndyr', '5tp', '5tp2ndyr', '6tp2ndyr', '7tp2ndyr', 'RTS,S'),
                             labels = c('infants', 'infants', 'extended', 'infants', 'extended', 'extended', 'extended', 'RTS,S'))

  return(dat)
}