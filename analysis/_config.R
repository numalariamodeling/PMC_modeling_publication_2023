##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## _config.R
##---------------------

tbl_opts <- c("striped", "hover", "condensed", "responsive") # kableExtrapackrat::init("~/projects/babynames")
pckg <- c("dplyr", "data.table", "wrapr", "tidyr", "lubridate", "zoo", "ggthemes", "cowplot",
          "scales", "ggplot2", "RColorBrewer", "see", 'pracma') # kableExtra
a <- lapply(pckg, require, character.only = TRUE)
rm(a)

plot_dir <- file.path('figures')
data_path <- file.path('data_files')
simout_dir <- file.path('postprocessed_output')

###-------------------  Define custom objects
eir_max = 200
max_years <- c(1, 2, 5, 10)
age_groups <- paste0('U', max_years)
age_label_values <- seq(0, 10, 2)
age_labels <- age_label_values
default_eir <- 32
default_seasons <- c('season1')  # c('comstant','season1','season2','season3')
default_Uage <- 'U1'
ipti_rtss_levels <- c('None', 'IPTi', 'RTS,S', 'IPTi + RTS,S')
scenario_labels <- c('None', 'PMC-3', 'PMC-4', 'PMC-5', 'PMC-6', 'PMC-7', 'RTS,S', 'PMC-3 + RTS,S')

###-------------------

touchpoints <- list()
touchpoints[['3tp']] <- round(c(10 / 4, 14 / 4, 9) * (365 / 12), 0)
touchpoints[['5tp2ndyr']] <- sort(round(c(touchpoints[['3tp']], c(12, 15) * (365 / 12)), 0))
touchpoints[['7tp']] <- sort(round(c(touchpoints[['3tp']], c(6, 12, 15, 18) * (365 / 12)), 0))


###-------------------  Define custom themes and colors
device_format <- c('pdf', 'png')

rtss_col <- '#FBB040'
pmc_cols <- c('#4EA3D1', '#007BBD', '#CBDB47', '#99CA3C', '#208B3A')
rtss_pmc_cols <- c('#fa8e9e', '#d97373')
getPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))(9)[3:9]
custom_cols <- c("#4EA3D1", "#007BBD", "#51571C", "#99CA3C", "#208B3A", rtss_col, '#d97373')
custom_cols2 <- c("#4EA3D1", "#51571C", "#208B3A", rtss_col, '#d97373')
RdYlBupalette <- colorRampPalette(brewer.pal(8, "RdYlBu"))(12)
YlGnpalette <- colorRampPalette(brewer.pal(9, "YlGnBu"))(9)

###------------------- Custom functions
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

len <- function(x) { return(length(x)) }

f_save_plot <- function(pplot, plot_name, plot_dir, width = 14, height = 8, units = 'in', device_format = c('pdf', 'png')) {
  if ('png' %in% device_format) {
    ggsave(paste0(plot_name, ".png"), plot = pplot, path = plot_dir,
           width = width, height = height, units = units, device = "png")
  }
  if ('pdf' %in% device_format) {
    if (!dir.exists(file.path(plot_dir, "pdf"))) { dir.create(file.path(plot_dir, "pdf")) }
    ggsave(paste0(plot_name, ".pdf"), plot = pplot, path = file.path(plot_dir, "pdf"),
           width = width, height = height, units = units, device = "pdf", useDingbats = FALSE)
  }
}

plot_combine <- function(plist, ncol = 1, legend_position = 'right', rel_dims = 1, leg_dim = 0.25, labels = c('', '')) {
  plegend <- get_legend(plist[[1]])

  for (i in c(1:length(plist))) {
    plist[[i]] <- plist[[i]] + theme(legend.position = 'None')
  }

  if (sum(lengths(rel_dims)) == 1)rel_dims <- rep(rel_dims, length(plist))
  pplot <- cowplot::plot_grid(plotlist = plist, ncol = ncol, labels = labels, rel_widths = rel_dims)
  if (tolower(legend_position) == 'right')pplot <- plot_grid(pplot, plegend, ncol = 2, rel_widths = c(1, leg_dim))
  if (tolower(legend_position) == 'left')pplot <- plot_grid(plegend, pplot, ncol = 1, rel_heights = c(leg_dim, 1))
  if (tolower(legend_position) == 'bottom')pplot <- plot_grid(pplot, plegend, ncol = 1, rel_heights = c(1, leg_dim))
  if (tolower(legend_position) == 'top')pplot <- plot_grid(plegend, pplot, ncol = 1, rel_heights = c(leg_dim, 1))
  if (tolower(legend_position) == 'none')pplot <- pplot
  return(pplot)
}


f_getCustomTheme <- function(fontscl = 1) {

  customTheme <- theme(
    strip.text.x = element_text(size = 12 * fontscl, face = "bold"),
    strip.text.y = element_text(size = 12 * fontscl, face = "bold"),
    strip.background = element_blank(),
    plot.title = element_text(size = 14 * fontscl, vjust = -1, hjust = 0, color = 'black'),
    plot.subtitle = element_text(size = 12 * fontscl, color = 'black'),
    plot.caption = element_text(size = 9 * fontscl, color = 'black'),
    legend.title = element_text(size = 12 * fontscl, color = 'black'),
    legend.text = element_text(size = 12 * fontscl, color = 'black'),
    axis.title.x = element_text(size = 12 * fontscl, color = 'black'),
    axis.text.x = element_text(size = 12 * fontscl, color = 'black'),
    axis.title.y = element_text(size = 12 * fontscl, color = 'black'),
    axis.text.y = element_text(size = 12 * fontscl, color = 'black'),
    axis.ticks = element_line(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )

  return(customTheme)
}


theme_set(theme_bw())
customTheme <- f_getCustomTheme() + guides(panel.spacing = unit(1.3, "lines"))
customTheme_nogrid <- customTheme + theme(panel.grid = element_blank())


gen_pmc_mode_fct <- function(dat) {
  pmc_levels <- c('3tp', '4tp', '4tp2ndyr', '5tp', '5tp2ndyr', '6tp2ndyr', '7tp2ndyr', 'RTS,S')
  pmc_labels <- c('PMC-3', 'PMC-4', 'PMC-4', 'PMC-5', 'PMC-5', 'PMC-6', 'PMC-7', 'RTS,S')
  dat$pmc_mode_fct <- factor(dat$pmc_mode, levels = pmc_levels, labels = pmc_labels)
  return(dat)
}