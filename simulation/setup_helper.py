import os
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import scipy.stats as stats
import copy
from dtk.interventions.input_EIR import add_InputEIR, monthly_to_daily_EIR
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from malaria.site.input_EIR_by_site import mAb_vs_EIR
from malaria.reports.MalariaReport import add_filtered_report
from malaria.reports.MalariaReport import add_event_counter_report
from hbhi.utils import add_monthly_parasitemia_rep_by_year, add_nmf_trt, generate_multiyr_df, read_main_dfs, tryread_df
from intervention_suite_cohort import InterventionSuite, add_all_interventions
from malaria.interventions.malaria_drugs import get_drug_param


def tryread_df2(input_path, scen_df, id, fname):
    try:
        df = tryread_df(os.path.join(input_path, scen_df.at[id, fname]))
    except:
        # print(f"WARNING: {fname} not in scen_df.")
        df = pd.DataFrame()
    return df


def setup_setting(cb,  scen_df, id, eir_monthly_multipliers, EIR_scale='monthly', cohort_month_shift=0, Maternal_Antibody_Protection=0.1327):

    if 'pmc_coverage' not in scen_df.columns:
        scen_df['pmc_coverage'] = 0
    if 'rtss_coverage' not in scen_df.columns:
        scen_df['rtss_coverage'] = 0
    if 'rtss_mode' not in scen_df.columns:
        scen_df['rtss_mode'] = 'constant'
    if 'rtss_booster1_min_age' not in scen_df.columns:
        scen_df['rtss_booster1_min_age'] = 730
    if 'pmc_mode' not in scen_df.columns:
        scen_df['pmc_mode'] = ''
    if 'rtss_target_group' not in scen_df.columns:
        scen_df['rtss_target_group'] = 'random'
    if 'cm_target_group' not in scen_df.columns:
        scen_df['cm_target_group'] = 'random'
    if 'pmc_target_group' not in scen_df.columns:
        scen_df['pmc_target_group'] = 'random'

    scen_row = scen_df[scen_df['setting_id'] == id]

    # DEMOGRAPHICS
    cb.update_params({'Demographics_Filenames': [os.path.join('generic_cohort_ds',scen_row['demographics_filename'][0])],
                      'Age_Initialization_Distribution_Type': 'DISTRIBUTION_SIMPLE'})

    annual_eir = float(scen_row['annual_EIR'][0])
    monthly_eir_scalers0 = eir_monthly_multipliers[scen_row['seasonality'][0]].tolist()
    # shift monthly EIRs to account for cohort born cohort_month_shift months into the year
    if cohort_month_shift > 0:
        monthly_eir_scalers = monthly_eir_scalers0[cohort_month_shift:12] + monthly_eir_scalers0[0:cohort_month_shift]
    else:
        monthly_eir_scalers = monthly_eir_scalers0
    eir_sum = sum([x for x in monthly_eir_scalers])

    if annual_eir is None:
        annual_eir = eir_sum
        monthly_eir = monthly_eir_scalers
    else:
        annual_eir = annual_eir
        monthly_eir = [(x / eir_sum) * annual_eir for x in monthly_eir_scalers]

    if EIR_scale == 'monthly':
        add_InputEIR(cb, start_day=0, monthlyEIRs=monthly_eir)  # , EIR_type='MONTHLY', , age_dependence='OFF'
    if EIR_scale == 'daily':
        daily_eir = monthly_to_daily_EIR(monthly_eir)
        daily_eir = [(x / sum(daily_eir)) * annual_eir for x in daily_eir]
        add_InputEIR(cb, start_day=0, EIR_type='DAILY', dailyEIRs=daily_eir)  # , age_dependence='OFF'

    """Adjust Maternal_Antibody_Protection for transmission intensity"""
    #Maternal_Antibody_Protection = 0.1327
    mAb = Maternal_Antibody_Protection * mAb_vs_EIR(annual_eir)
    cb.update_params({'Maternal_Antibody_Protection': mAb})

    return {'seasonality': scen_row['seasonality'][0],
            'Annual EIR': annual_eir,
            'mAb': mAb,
            'cm_coverage': scen_row['cm_coverage'][0],
            'rtss_coverage': scen_row['rtss_coverage'][0],
            'pmc_coverage': scen_row['pmc_coverage'][0],
            'rtss_target_group': scen_row['rtss_target_group'][0],
            'cm_target_group': scen_row['cm_target_group'][0],
            'pmc_target_group': scen_row['pmc_target_group'][0],
            'frac_high_access': scen_row['frac_high_access'][0],
            'rtss_mode': scen_row['rtss_mode'][0],
            'minBoostAge': scen_row['rtss_booster1_min_age'][0],
            'pmc_mode': scen_row['pmc_mode'][0]
            }


def setup_simulation(years):
    # BASIC SETUP
    cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')

    # Run time
    cb.update_params({'Simulation_Duration': years * 365 + 1})

    # Logging
    cb.update_params({
        'logLevel_JsonConfigurable': 'ERROR',
        'logLevel_VectorHabitat': 'ERROR',
        'logLevel_StandardEventCoordinator': 'ERROR',
        'logLevel_SusceptibilityMalaria': 'ERROR'
    })

    # Demographics
    cb.update_params({
        'Enable_Birth': 0,
        'Enable_Births': 0,
        'Enable_Demographics_Birth': 0,
        'Enable_Demographics_Risk': 0,
        'Birth_Rate_Dependence': 'FIXED_BIRTH_RATE',
        'Enable_Initial_Prevalence': 0,
        'Age_Dependent_Biting_Risk_Type': 'OFF',
        'x_Birth': 1,
        'x_Base_Population': 1,
        'Maternal_Antibodies_Type': 'CONSTANT_INITIAL_IMMUNITY',
        "Age_Initialization_Distribution_Type": 'DISTRIBUTION_COMPLEX',
        'Disable_IP_Whitelist': 1
    })

    # Serialization (none)
    cb.update_params({
        'Serialization_Type': 'NONE',
        'Serialized_Population_Writing_Type': 'NONE',
        'Serialized_Population_Reading_Type': 'NONE',
    })

    # Vector and climate
    cb.update_params({
        "Vector_Species_Names": [],
        'x_temporary_Larval_Habitat': 0,
        'Climate_Model': 'CLIMATE_CONSTANT'  # "CLIMATE_BY_DATA"
    })

    # Reporting
    cb.update_params({
        'Enable_Default_Reporting': 0,
        'Enable_Property_Output': 0,
        'Enable_Vector_Species_Report': 0,
        'Report_Detection_Threshold_Blood_Smear_Parasites': 10,  # 50
        "Parasite_Smear_Sensitivity": 0.02,  # 50/uL
        'RDT_Sensitivity': 0.1
    })
    cb.update_params({
        "Report_Event_Recorder": 1,
        "Report_Event_Recorder_Individual_Properties": ['VaccineStatus', 'AccessToInterventions'],
        "Report_Event_Recorder_Ignore_Events_In_List": 0,
        "Report_Event_Recorder_Events": ['Births', 'PropertyChange'],
        'Custom_Individual_Events': ['Received_Treatment', 'Received_Severe_Treatment',
                                     'Received_NMF_Treatment', 'Received_Vaccine', 'Received_Campaign_Drugs']
    })

    # Filtered report (all years):
    num_year = 1
    start = 1  # 1 + (years - num_year) * 365
    end = 1 + years * 365
    add_filtered_report(cb, start=start, end=end)

    # CUSTOM REPORTS
    add_event_counter_report(cb, event_trigger_list=['Received_Treatment', 'Received_Severe_Treatment',
                                                     'Received_Vaccine', 'Received_Campaign_Drugs'],
                             duration=years * 365 + 1)
    #
    # for year in range(years):
    #     add_monthly_parasitemia_rep_by_year(cb, num_year=years, tot_year=years,
    #                                         sim_start_year=2020,
    #                                         yr_plusone=True, prefix='Monthly')
    # add_monthly_parasitemia_rep_by_year(cb, num_year=years, tot_year=years,
    #                                     sim_start_year=2020,
    #                                     yr_plusone=True,
    #                                     age_bins=[1, 5, 120],
    #                                     prefix='FineMonthly')
    return cb


def add_generic_interventions(cb, input_path, scen_df, id, ds_name='run_col', cohort_month_shift=0):
    # INTERVENTIONS
    int_suite = InterventionSuite()

    int_suite.hs_ds_col = ds_name
    int_suite.pmc_ds_col = ds_name
    int_suite.rtss_ds_col = ds_name

    hs_df = tryread_df2(input_path, scen_df, id, 'CM_filename')
    pmc_df = tryread_df2(input_path, scen_df, id, 'PMC_filename')
    rtss_df = tryread_df2(input_path, scen_df, id, 'RTSS_filename')

    if not pmc_df.empty:
        pmc_df['drug_code'] = scen_df['drug_code'].unique()[0]
        if 'num_IIV_groups' in scen_df.columns:
            pmc_df['num_IIV_groups'] = scen_df['num_IIV_groups'].unique()[0]

    add_all_interventions(cb,
                          int_suite,
                          my_ds='run',  # to not subset dataframe as each intervention csv is generic
                          high_access_ip_frac=scen_df.at[id, 'frac_high_access'],
                          rtss_target_group=scen_df.at[id, 'rtss_target_group'],
                          cm_target_group=scen_df.at[id, 'cm_target_group'],
                          pmc_target_group=scen_df.at[id, 'pmc_target_group'],
                          rtss_booster1_min_age=scen_df.at[id, 'rtss_booster1_min_age'],
                          hs_df=hs_df,
                          pmc_df=pmc_df,
                          rtss_df=rtss_df,
                          cohort_month_shift=cohort_month_shift)

    return {'Scenario_id': id}


def add_generic_interventions_nga(cb, input_path, scen_df, id, ds_name='run_col', cohort_month_shift=0):
    # INTERVENTIONS
    int_suite = InterventionSuite()

    int_suite.hs_ds_col = ds_name
    int_suite.pmc_ds_col = ds_name
    int_suite.rtss_ds_col = ds_name

    hs_df = tryread_df2(input_path, scen_df, id, 'CM_filename')
    pmc_df = tryread_df2(input_path, scen_df, id, 'PMC_filename')
    rtss_df = tryread_df2(input_path, scen_df, id, 'RTSS_filename')

    try:
        hs_df = hs_df.loc[hs_df['State'] == id]
    except:
        pass
    try:
        pmc_df = pmc_df.loc[pmc_df['State'] == id]
    except:
        pass
    try:
        rtss_df = rtss_df.loc[rtss_df['State'] == id]
    except:
        pass

    if not pmc_df.empty:
        pmc_df['drug_code'] = scen_df['drug_code'].unique()[0]
        if 'num_IIV_groups' in scen_df.columns:
            pmc_df['num_IIV_groups'] = scen_df['num_IIV_groups'].unique()[0]

    add_all_interventions(cb,
                          int_suite,
                          my_ds='run',  # to not subset dataframe as each intervention csv is generic
                          high_access_ip_frac=scen_df.at[id, 'frac_high_access'],
                          rtss_target_group=scen_df.at[id, 'rtss_target_group'],
                          cm_target_group=scen_df.at[id, 'cm_target_group'],
                          pmc_target_group=scen_df.at[id, 'pmc_target_group'],
                          rtss_booster1_min_age=scen_df.at[id, 'rtss_booster1_min_age'],
                          hs_df=hs_df,
                          pmc_df=pmc_df,
                          rtss_df=rtss_df,
                          cohort_month_shift=cohort_month_shift)

    return {'Scenario_id': id}
