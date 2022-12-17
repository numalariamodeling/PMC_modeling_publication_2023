##---------------------
## __run_all.R
##---------------------

fig01 <- T
fig02 <- T
fig03 <- T
fig04 <- T
fig05 <- T
figSI <- T

source(file.path('analysis', '_config.R'))

## Figure 1-4
if (fig01)source(file.path('analysis', 'Fig_01.R'))
if (fig02)source(file.path('analysis', 'Fig_02.R'))
if (fig03)source(file.path('analysis', 'Fig_03.R'))
if (fig04)source(file.path('analysis', 'Fig_04.R'))

## Figure 5
if (fig05) {
  source(file.path('analysis', '_fig05_helper_functions.R'))
  source(file.path('analysis', 'Fig_05a_extract_DHS.R'))
  source(file.path('analysis', 'Fig_05b_EPI_coverage.R'))
  source(file.path('analysis', 'Fig_05c_pfpr_to_eir.R'))
  source(file.path('analysis', 'Fig_05d_scenario_csvs.R'))
  source(file.path('analysis', 'Fig_05e_maps.R'))
  source(file.path('analysis', 'Fig_05f_figures.R'))
}


## Supplementary figures
if (figSI) {

}