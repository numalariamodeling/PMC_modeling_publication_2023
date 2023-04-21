##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## __run_all.R
##---------------------

fig01 <- T
fig02 <- T
fig03 <- T
fig04 <- T
fig05 <- T
figSI <- T

source(file.path('analysis', '_config.R'))

## Figures 1-4
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
  source(file.path('analysis', 'SI_Navrongo.R'))
  source(file.path('analysis', 'SI_fig_01_supp_rtss_by_eir.R'))
  source(file.path('analysis', 'SI_fig_01_supp_pmc_single.R'))
  source(file.path('analysis', 'SI_fig_02_scen_EIR_U1.R'))
  source(file.path('analysis', 'SI_fig_02_supp_age_incidence.R'))
  source(file.path('analysis', 'SI_fig_02_supp_agegroup.R'))
  source(file.path('analysis', 'SI_fig_04_supp.R'))
  source(file.path('analysis', 'SI_fig_incremental_impact_EIR.R'))
  source(file.path('analysis', 'SI_pmc_cm_agevariation.R'))
  source(file.path('analysis', 'SI_pmc_matAntibodies.R'))
  source(file.path('analysis', 'SI_population_demography.R'))
}