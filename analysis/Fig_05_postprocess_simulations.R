tbl_opts <- c("striped", "hover", "condensed", "responsive")
setwd(file.path(getwd(), 'manuscript_01'))
source(file.path('00_config.R'))
source(file.path('Fig_05_helper_functions.R'))

simout_dir <- gsub('_generic', '_NGA', simout_dir)
plot_dir <- file.path(plot_dir, 'fig_5_additional')
outcome_channels = c('PE_clinical_incidence', 'PE_severe_incidence')

exp_name_counterfactual <- "NGA_counterfactual_vaccSP_IIV"
exp_names <- unique(c(exp_name_counterfactual, list.files(simout_dir)))
print(exp_names)

df_c <- f_counterfactual_pfpr(exp_name_counterfactual)

df_list = list()
for (exp_name in c(exp_names)) {
  if (file.exists(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv'))) {
    df_list[[length(df_list) + 1]] <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
      filter(age_group %in% c('U5', 'U2', 'U1')) %>%
      filter(name %in% c('PfHRP2_Prevalence', 'clinical_cases', 'severe_cases')) %>%
      pivot_longer(cols = c(mean_val, median_val, low_val, up_val, min_val, max_val), names_to = 'statistic') %>%
      pivot_wider(names_from = name,
                  values_from = c(value),
                  names_glue = "{name}_{.value}") %>%
      rename_with(~gsub("_value", "", .x)) %>%
      mutate(pmc_coverage = as.character(pmc_coverage),
             rtss_coverage = as.character(rtss_coverage),
             scen = gsub('_vaccSP_IIV', '', gsub('NGA_', '', exp_name))) %>%
      mutate(Annual_EIR = round(Annual_EIR, 2),
             cm_coverage = round(cm_coverage, 2))
  }
}

#Joining, by = c("age_group", "seasonality", "Annual_EIR", "cm_coverage", "statistic")

df <- df_list %>%
  bind_rows() %>%
  left_join(df_c) %>%
  mutate(PE_clinical_cases = 1 - (clinical_cases / clinical_cases_c),
         PE_severe_cases = 1 - (severe_cases / severe_cases_c),
         clinical_cases_averted = clinical_cases_c - clinical_cases,
         severe_cases_averted = severe_cases_c - severe_cases)

tapply(df$clinical_cases, df$scen, summary)
tapply(df$severe_cases, df$scen, summary)


target_scen <- unique(df$scen[grep('targetcov', df$scen)])
df <- df %>%
  mutate(coveragemode = ifelse(scen %in% target_scen, 'target', 'operational'),
         scen = gsub('_targetcov', '', scen))

df$scen <- factor(df$scen,
                  levels = c('counterfactual', 'pmc_3tp', 'pmc_5tp2ndyr', 'pmc_7tp2ndyr', 'rtss', 'rtss_pmc_3tp'),
                  labels = c('counterfactual', 'PMC-3', 'PMC-5', 'PMC-7', 'RTS,S', 'PMC-3 + RTS,S'))
fwrite(df, file.path(simout_dir, 'NGA_simdat_aggr_agegroup.csv'))


###------------------
prop_popU2_totalpop <- 0.068
require(readxl)
nga_pop_all <- read_xlsx(file.path(pmc_path, 'NGA', 'population_size', 'GEOPODE_NGA_population_ages0to100_version2022.xlsx'), skip = 6, trim_ws = T,
                         .name_repair = "universal") %>%
  dplyr::select_if(~!(all(is.na(.)) | all(. == ""))) %>%
  filter(!(is.na(LGA))) %>%
  rename_with(~gsub('[...]', '', .x)) %>%
  group_by(State) %>%
  summarise(total_pop = sum(Total7)) %>%
  mutate(pop_U2 = round(total_pop * prop_popU2_totalpop, 0),
         State = gsub('Akwa lbom', 'Akwa Ibom', State))
unique(nga_pop$State)

nga_pop <- nga_pop_all %>% filter(State %in% df$seasonality)
unique(nga_pop$State)

sum(nga_pop$total_pop)
sum(nga_pop$pop_U2)

unique(nga_pop$State)
unique(df$seasonality)

df <- fread(file.path(simout_dir, 'NGA_simdat_aggr_agegroup.csv'))
df$scen <- factor(df$scen,
                  levels = c('counterfactual', 'PMC-3', 'PMC-5', 'PMC-7', 'RTS,S', 'PMC-3 + RTS,S'),
                  labels = c('counterfactual', 'PMC-3', 'PMC-5', 'PMC-7', 'RTS,S', 'PMC-3 + RTS,S'))


df <- df %>%
  mutate(State = seasonality) %>%
  filter(statistic == 'mean_val' &
           age_group == 'U2') %>%
  left_join(nga_pop) %>%
  mutate(clinicalinc_ppa = clinical_cases / 1000,
         clinicalcases_pop = clinical_cases / 1000 * pop_U2,
         severeinc_ppa = severe_cases / 1000,
         severecases_pop = severe_cases / 1000 * pop_U2,
         clinical_cases_averted_pop = clinical_cases_averted / 1000 * pop_U2, ,
         severe_cases_averted_pop = severe_cases_averted / 1000 * pop_U2) %>%
  mutate(ADM1_NAME = gsub('Akwa Ibom', 'Akwa lbom', State))

tapply(df$clinical_cases_averted_pop, df$scen, summary)


table(df$scen)
tapply(df$PE_clinical_cases, df$scen, summary)
tapply(df$PE_severe_cases, df$scen, summary)

tapply(df$PE_clinical_cases, df$age_group, summary)
tapply(df$PE_severe_cases, df$age_group, summary)


ggplot(data = subset(df, age_group == 'U2')) +
  geom_boxplot(aes(x = scen, y = PE_clinical_cases, group = interaction(coveragemode, scen), fill = coveragemode))

ggplot(data = subset(df, age_group == 'U2')) +
  geom_boxplot(aes(x = scen, y = clinical_cases_averted, group = interaction(coveragemode, scen), fill = coveragemode))

grp_vars <- qc(age_group, scen, pmc_mode, coveragemode)
outcome_vars <- qc(clinical_cases, severe_cases, clinical_cases_averted, severe_cases_averted,
                   PE_clinical_cases, PE_severe_cases,
                   clinical_cases_averted_pop, severe_cases_averted_pop)

plotdat <- df %>%
  filter(age_group == 'U2' &
           scen != 'counterfactual' &
           statistic == 'mean_val') %>%
  dplyr::select_at(.vars = c(grp_vars, outcome_vars)) %>%
  tidyr::pivot_longer(col = outcome_vars) %>%
  dplyr::group_by_at(.vars = c(grp_vars, 'name')) %>%
  dplyr::summarise(n.val = n(),
                   sd.val = sd(value, na.rm = TRUE),
                   mean_val = mean(value),
                   median_val = median(value),
                   low_val = quantile(value, probs = 0.05, na.rm = TRUE),
                   up_val = quantile(value, probs = 0.95, na.rm = TRUE),
                   min_val = min(value),
                   max_val = max(value)) %>%
  dplyr::mutate(
    se.val = sd.val / sqrt(n.val),
    lower.ci.val = mean_val - qt(1 - (0.05 / 2), n.val - 1) * se.val,
    upper.ci.val = mean_val + qt(1 - (0.05 / 2), n.val - 1) * se.val,
    weighted = 0)

line_cols <- rep(custom_cols2, 2)
fill_cols <- line_cols
fill_cols[c(1:5)] <- 'white' #fill_cols[c(6:10)] <- 'white'

figleg1 <- get_legend(ggplot(data = subset(plotdat, name == 'clinical_cases_averted')) +
                        geom_col(aes(x = scen, y = mean_val, fill = scen)) +
                        scale_y_continuous(lim = c(0, 600), expand = c(0, 0)) +
                        scale_fill_manual(values = custom_cols2) +
                        labs(color = 'Scenario', fill = 'Scenario'))

figleg2 <- get_legend(ggplot(data = subset(plotdat, name == 'clinical_cases_averted')) +
                        geom_col(aes(x = scen, y = mean_val, fill = coveragemode, col = coveragemode)) +
                        scale_y_continuous(lim = c(0, 600), expand = c(0, 0)) +
                        scale_color_manual(values = c('black', 'black')) +
                        scale_fill_manual(values = c('white', 'black')) +
                        labs(color = 'Coverage', fill = 'Coverage'))

pp <- ggplot(data = subset(plotdat, name == 'clinical_cases_averted_pop')) +
  geom_col(aes(x = scen, y = mean_val,
               group = interaction(coveragemode, scen),
               col = interaction(scen, coveragemode),
               fill = interaction(scen, coveragemode)),
           position = position_dodge(width = 0.6), width = 0.7, show.legend = F) +
  geom_errorbar(aes(x = scen, ymin = lower.ci.val, ymax = upper.ci.val,
                    group = interaction(coveragemode, scen)),
                position = position_dodge(width = 0.6), width = 0) +
  #scale_y_continuous(lim = c(0, 600), expand = c(0, 0),labels=comma) +
  scale_y_continuous(lim = c(0, 180000), expand = c(0, 0), labels = comma) +
  scale_fill_manual(values = fill_cols) +
  scale_color_manual(values = line_cols) +
  labs(x = '', y = 'Clinical cases per total population U2') +
  customTheme_nogrid

figleg <- plot_grid(NULL, figleg1, figleg2, NULL, ncol = 1, rel_heights = c(1, 0.6, 0.6, 1), align = 'v')

pplot <- plot_grid(pp, figleg2, rel_widths = c(1, 0.3))
pplot
f_save_plot(pplot, paste0('Fig5_clinical_cases'), file.path(plot_dir), width = 8, height = 3, units = 'in', device_format = device_format)


pp <- ggplot(data = subset(plotdat, name == 'severe_cases_averted_pop')) +
  geom_col(aes(x = scen, y = mean_val,
               group = interaction(coveragemode, scen),
               col = interaction(scen, coveragemode),
               fill = interaction(scen, coveragemode)),
           position = position_dodge(width = 0.6), width = 0.7, show.legend = F) +
  geom_errorbar(aes(x = scen, ymin = lower.ci.val, ymax = upper.ci.val,
                    group = interaction(coveragemode, scen)),
                position = position_dodge(width = 0.6), width = 0) +
  scale_y_continuous(lim = c(0, 8), expand = c(0, 0)) +
  scale_y_continuous(lim = c(0, 2500), expand = c(0, 0), labels = comma) +
  scale_fill_manual(values = fill_cols) +
  scale_color_manual(values = line_cols) +
  labs(x = '', y = 'Severe cases per total population') +
  customTheme_nogrid

figleg <- plot_grid(NULL, figleg1, figleg2, NULL, ncol = 1, rel_heights = c(1, 0.6, 0.6, 1), align = 'v')
pplot <- plot_grid(pp, figleg2, rel_widths = c(1, 0.3))
f_save_plot(pplot, paste0('Fig5_severe_cases'), file.path(plot_dir), width = 8, height = 3, units = 'in', device_format = device_format)





