##---------------------
## Perennial malaria chemoprevention with and without malaria vaccination to reduce malaria burden in young children: a modeling analysis
## _fig05_helper_functions.R
##---------------------

pckg <- c("rdhs", "raster", "malariaAtlas", "exactextractr", "wpgpDownloadR", "RColorBrewer")
a <- lapply(pckg, require, character.only = TRUE)
rm(a)

##---------------
## Helper functions
##---------------
add_nga_admin <- function(df) {

  pmc_states <- c("Abia", "Akwa lbom", "Akwa Ibom", "Anambra", "Bayelsa", "Benue", "Enugu", "Cross River", "Delta", "Oyo", "Ebonyi", "Edo", "Ekiti", "Imo",
                  "Kogi", "Lagos", "Ogun", "Ondo", "Osun", "Rivers", "Taraba")

  df <- df %>% mutate(southern_nga = ifelse(NAME_1 %in% pmc_states, 'Southern States', 'Northern States'),
                      geo_zone = case_when(NAME_1 %in% qc(Oyo, Ogun, Lagos, Osun, Ondo, Ekiti) ~ 'South-West',
                                           NAME_1 %in% qc(Edo, Delta, Bayelsa, Rivers, `Akwa Ibom`, `Akwa lbom`, `Cross River`) ~ 'South-South',
                                           NAME_1 %in% qc(Anambra, Enugu, Imo, Abia, Ebonyi) ~ 'South-East',
                                           NAME_1 %in% qc(Niger, Kwara, Kogi, `FCT Abuja`, Nasarawa, Benue, Plateau) ~ 'North-Central',
                                           NAME_1 %in% qc(Kebbi, Sokoto, Zamfara, Katsina, Kano, Jigawa, Kaduna) ~ 'North-West',
                                           NAME_1 %in% qc(Bauchi, Gombe, Yobe, Borno, Adamawa, Taraba) ~ 'North-East',

                      ))
  return(df)
}

dhs_get_name1 <- function(dhs_df) {
  levels_to_keep <- grep("[..]", dhs_df$CharacteristicLabel)
  dhs_df <- dhs_df[levels_to_keep,]
  dhs_df$NAME_1 <- gsub("[..]", "", dhs_df$CharacteristicLabel)
  return(dhs_df)
}

EPI_to_pmc_cov <- function(x = NULL, get_scaling_factors = T) {
  #https://malariajournal.biomedcentral.com/articles/10.1186/s12936-021-03615-3
  epicov_mean <- c(0.804, 0.654, 0.522)
  ipticov_mean <- c(0.674, 0.625, 0.364)
  downscaling <- ipticov_mean / epicov_mean

  if (get_scaling_factors == FALSE & !is.null(x)) {
    if (length(x) != length(downscaling)) {
      print("Not same length!")
    }else {
      print("Returning downscaled coverage values")
      x <- x * downscaling
    }
  }else {
    print("Returning downscaling factors")
    x <- downscaling
  }

  return(x)

}


##---------------
## Interpolations
##---------------
f_spline <- function(dat, x, y) {
  dat <- as.data.frame(dat)
  sp_model <- smooth.spline(dat[, x], dat[, y])
  return(sp_model)
}

f_get_pred_dat <- function(pdat) {
  dfGLM <- pdat %>%
    group_by(cm_coverage) %>%
    do(fitglm = f_spline(dat = ., x = 'pfpr_U5', y = 'mean_val'))

  ### Make predictions
  pred_list <- list()
  for (i in c(1:nrow(dfGLM))) {
    pred_dat <- data.frame('x' = seq(min(pdat[, 'pfpr_U5']), max(pdat[, 'pfpr_U5']), 0.025), 'y_pred' = NA)
    sp_model <- dfGLM[i, 'fitglm'][[1]]
    pred_dat$y_pred <- predict(sp_model[[1]], x = pred_dat$x)$y
    pred_dat$cm_coverage <- dfGLM[i, 'cm_coverage'][[1]]
    pred_list[[length(pred_list) + 1]] <- pred_dat
  }

  ### Combine data-list
  pred_dat <- pred_list %>%
    bind_rows() %>%
    dplyr::rename(pfpr_U5_mean = x, mean_val = y_pred) %>%
    as.data.table()

  return(pred_dat)
}


get_interpolatedgrid <- function(dat, cov1, cov2, outcome, return_model = F) {
  ### Define regression model per group
  dfLM <- dat %>%
    as.data.table() %>%
    rename(cov1 = as.name(cov1),
           cov2 = as.name(cov2)) %>%
    dplyr::group_by(cm_coverage, pmc_mode) %>%
    do(fitglm = glm(get(outcome) ~ cov1 + cov2, family = 'gaussian', data = .))

  ### Make predictions
  pred_list <- list()
  for (i in c(1:nrow(dfLM))) {
    pred_dat <- data.frame(expand.grid('cov1' = seq(0, 1, 0.05), 'cov2' = seq(0, 1, 0.05)))
    pred_dat$y_pred <- NA

    model <- dfLM[i, 'fitglm'][[1]]
    pred_dat$y_pred <- predict(model[[1]], newdata = pred_dat, type = "response")
    pred_dat$cm_coverage <- dfLM[i, 'cm_coverage'][[1]]
    pred_dat$pmc_mode <- dfLM[i, 'pmc_mode'][[1]]
    pred_dat[, outcome] <- pred_dat[, 'y_pred']
    pred_dat[, cov1] <- pred_dat[, 'cov1']
    pred_dat[, cov2] <- pred_dat[, 'cov2']
    pred_dat <- pred_dat %>% dplyr::select(-cov1, -cov2, -y_pred)
    pred_list[[length(pred_list) + 1]] <- pred_dat
  }

  ### Combine data-list
  pred_dat <- pred_list %>% bind_rows()

  if (return_model)outdat <- dfLM
  if (return_model == FALSE)outdat <- pred_dat
  return(outdat)
}


##---------------
## Data postprocessing, effect size
##---------------
f_counterfactual_pfpr <- function(exp_name, SAVE = T) {
  counterfactual_pfpr <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
    filter(age_group %in% c('U5', 'U2', 'U1') &
             pmc_coverage == 0 &
             rtss_coverage == 0) %>%
    filter(name %in% c('PfHRP2_Prevalence', 'clinical_cases', 'severe_cases')) %>%
    mutate(name = paste0(name, '_c')) %>%
    pivot_longer(cols = c(mean_val, median_val, low_val, up_val, min_val, max_val), names_to = 'statistic') %>%
    pivot_wider(names_from = name,
                values_from = value,
                names_glue = "{name}_{.value}") %>%
    dplyr::select(-pmc_coverage, -rtss_coverage, -pmc_mode, -pmc_rtss_cov) %>%
    rename_with(~gsub("_value", "", .x)) %>%
    mutate(Annual_EIR = round(Annual_EIR, 2),
           cm_coverage = round(cm_coverage, 2))

  if (SAVE)fwrite(counterfactual_pfpr, file.path(simout_dir, exp_name, 'counterfactual_pfpr.csv'))
  return(counterfactual_pfpr)
}

f_effect_size <- function() {
  outcome_channels = c('clinical_cases', 'severe_cases')
  pmc_effect_dat <- fread(file.path(simout_dir, exp_name, 'simdat_aggr_agegroup.csv')) %>%
    filter(age_group == 'U2' & (pmc_coverage != 0 | rtss_coverage != 0)) %>%
    filter(name %in% outcome_channels) %>%
    pivot_longer(cols = c(mean_val, median_val, low_val, up_val, min_val, max_val), names_to = 'statistic') %>%
    pivot_wider(names_from = age_group, values_from = value,
                names_glue = "{age_group}_{.value}") %>%
    pivot_wider(names_from = name,
                values_from = c(U2_value),
                names_glue = "{name}_{.value}") %>%
    rename_with(~gsub("_value", "", .x)) %>%
    left_join(counterfactual_pfpr)

  return(pmc_effect_dat)

}