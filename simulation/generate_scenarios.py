import os
import sys
import numpy as np
import pandas as pd
import itertools

sys.path.append('../')
from load_paths import load_box_paths


def get_parameter_space():
    """all parameter until deployment_parameters are full factorial combinations """
    parameter_space = {
        'pop_size': 30000,  # 20000
        'annual_EIR': [5, 10, 30],
        'cm_coverage': [0.6],
        'seasonality': ['high_unimodal'],  # ['constant', 'moderate_unimodal', 'high_unimodal'],
        'rtss_target_group': ['random'],  # specify whether RTS,S is targeted toward an access group or given at random
        'cm_target_group': ['random'],  # specify whether CM is targeted toward an access group or given at random
        'pmc_target_group': ['random'],  # specify whether PMC is targeted toward an access group or given at random
        'rtss_coverage': [0, 0.8],
        'rtss_mode': ['constant'],  # options are: 'constant' (EPI main dose and booster dose),
        #              'campboost' (EPI main dose and 1 campaign booster dose),
        #              'campboost2' (EPI main dose and 2 campaign booster doses)
        'pmc_coverage': [0, 0.8],
        'rtss_age_days': [274, 730],  # per default age 9 and 24 months for 2st EPI round, and 1st booster
        'rtss_booster1_min_age': [730],  # per default minimum age 24 months to receive 1st booster
        'booster_coverage_rtss': 0.8,
        'max_age_rtss': 5,
        'decay_time_constant_rtss': 592.4066512,
        'decay_class_rtss': 'WaningEffectExponential',
        'pmc_mode': ['3tp'],  # options are: '3tp', '5tp' or custom
        'pmc_touchpoints': [int(round(x * (365 / 12), 0)) for x in [10 / 4, 14 / 4, 9]],
        'season_start_month': 7,
    }
    return parameter_space


def f_cov(x):
    return int(float(x) * 100)


def create_scenarios_mastercsv(param_dic, fname_out, high_access_frac=0, return_csv=False):
    annual_EIR = param_dic['annual_EIR']
    cm_coverage = param_dic['cm_coverage']
    seasonality = param_dic['seasonality']
    rtss_coverage = param_dic['rtss_coverage']
    rtss_mode = param_dic['rtss_mode']
    rtss_booster1_min_age = param_dic['rtss_booster1_min_age']
    pmc_coverage = param_dic['pmc_coverage']
    pmc_mode = param_dic['pmc_mode']
    rtss_target_group = param_dic['rtss_target_group']
    cm_target_group = param_dic['cm_target_group']
    pmc_target_group = param_dic['pmc_target_group']

    df_array = np.array(
        list(itertools.product(annual_EIR, cm_coverage, seasonality,
                               rtss_target_group, cm_target_group, pmc_target_group,
                               rtss_mode, rtss_coverage, pmc_coverage, pmc_mode,
                               rtss_booster1_min_age)))
    df = pd.DataFrame(df_array)
    # Optionally drop certain combinations
    # [...]

    # Prepare and save df
    df.columns = ['annual_EIR', 'cm_coverage', 'seasonality',
                  'rtss_target_group', 'cm_target_group', 'pmc_target_group',
                  'rtss_mode', 'rtss_coverage', 'pmc_coverage', 'pmc_mode', 'rtss_booster1_min_age']

    df.reset_index(inplace=True)
    df = df.rename(columns={'index': 'setting_id'})
    df['setting_id'] = 'HX' + df['setting_id'].astype(str)

    # for the correlated demographics file, use the minimum of the coverages as the fraction of people in the high-access group unless otherwise specified
    if type(high_access_frac) == int or type(high_access_frac) == float:
        df['frac_high_access'] = high_access_frac
    elif high_access_frac == 'cm_coverage':
        df['frac_high_access'] = df['cm_coverage']
        df['frac_high_access'] = df['frac_high_access'].astype(str).astype(float)
        df.loc[(df['frac_high_access'] == 1.0), 'frac_high_access'] = 0
    elif high_access_frac == 'rtss_coverage':
        df['frac_high_access'] = df['rtss_coverage']
        df['frac_high_access'] = df['frac_high_access'].astype(str).astype(float)
        df.loc[(df['frac_high_access'] == 1.0), 'frac_high_access'] = 0
    elif high_access_frac == 'pmc_coverage':
        df['frac_high_access'] = df['pmc_coverage']
        df['frac_high_access'] = df['frac_high_access'].astype(str).astype(float)
        df.loc[(df['frac_high_access'] == 1.0), 'frac_high_access'] = 0
    else:
        df['frac_high_access'] = [np.min([float(df['cm_coverage'][yy]),
                                          (float(df['rtss_coverage'][yy]) if float(df['rtss_coverage'][yy]) > 0 else 1),
                                          (float(df['pmc_coverage'][yy]) if float(
                                              df['pmc_coverage'][yy]) > 0 else 1)])
                                  for yy in range(len(df))]

    # add demographics, seasonality, and intervention csv filepaths
    df['demographics_filename'] = [os.path.join('Demographics',
                                                f'generic_demographics_cohort_correlated{f_cov(df["frac_high_access"][yy])}_IIV.json')
                                   for yy in range(len(df))]
    df.loc[df.frac_high_access == 0, 'demographics_filename'] = os.path.join('Demographics',
                                                                             'generic_demographics_cohort_uncorrelated_IIV.json')

    df['CM_filename'] = [os.path.join('CM', f'CM_constant_{f_cov(yy)}coverage.csv') for yy in
                         df['cm_coverage']]

    try:
        df['RTSS_filename'] = [os.path.join('RTSS', f'RTSS_{xx}_{f_cov(yy)}coverage.csv') for xx, yy in
                               zip(df['rtss_mode'], df['rtss_coverage'])]
    except:
        df['RTSS_filename'] = [os.path.join('RTSS', f'RTSS_{xx}_{yy}coverage.csv') for xx, yy in
                               zip(df['rtss_mode'], df['rtss_coverage'])]
    try:
        df['PMC_filename'] = [os.path.join('PMC', f'PMC_{xx}_{f_cov(yy)}coverage.csv') for xx, yy in
                              zip(df['pmc_mode'], df['pmc_coverage'])]
    except:
        df['PMC_filename'] = [os.path.join('PMC', f'PMC_{xx}_{(yy)}coverage.csv') for xx, yy in
                              zip(df['pmc_mode'], df['pmc_coverage'])]

    print(f'Saving  {fname_out}.csv')
    ## check whether all intervention csvs exists, print out those not exciting
    print(check_intervention_csvs_exist(df, base_scenario_filepath))
    df.to_csv(os.path.join(base_scenario_filepath, f'{fname_out}.csv'),
              index=False)

    if return_csv:
        print(f'Returning dataframe for {fname_out}.csv')
        return df


def check_intervention_csvs_exist(df, fpath):
    filenames = [col for col in df.columns if 'filename' in col]
    filenames_not_exits = []
    for fname in filenames:
        unique_files = list(df[fname].unique())
        filenames_not_exits = [x for x in unique_files if not os.path.exists(os.path.join(fpath, x))]

    return f'filenames_not_exits: {filenames_not_exits}'


def intervention_inputs(param_dic, CM=False, RTSS=False, PMC=False):
    if CM == True:
        cm_coverage = param_dic['cm_coverage']
        for cm in cm_coverage:
            if cm == 0:
                df = pd.DataFrame()
            else:
                df = pd.DataFrame({'U5_coverage': [cm], 'adult_coverage': [cm], 'severe_cases': np.min([0.8, cm * 2]),
                                   'simday': [0], 'duration': [-1], 'start_day_override': [-1], 'run_col': 'run'})
            print(f'Writing CM_constant_{int(100 * float(cm))}coverage.csv')
            df.to_csv(os.path.join(base_scenario_filepath, 'CM', 'CM_constant_%icoverage.csv' % (100 * float(cm))))

    if RTSS == True:
        # RTS,S parameters
        for rtss_mode in param_dic['rtss_mode']:
            rtss_age_days = param_dic['rtss_age_days']
            if 'campboost' in rtss_mode:
                if rtss_mode == 'campboost' or rtss_mode == 'campboost_7m':
                    num_rtss_rounds = 2  # main EPI round and 1 campaign booster
                    rtss_types = ['simple', 'booster1']
                else:
                    if 'campboost_flexible' in rtss_mode:
                        num_rtss_rounds = 4  # main EPI round and 1 booster within 3 years with defined coverages=
                        rtss_types = ['simple', 'booster1', 'booster2',
                                      'booster3']  # 1,2,3 for years to receive booster within 2-5 years
                    elif 'campboost2' in rtss_mode:
                        if 'campboost2_flexible' in rtss_mode:
                            num_rtss_rounds = 5  # main EPI round and 1 booster within 3 years with defined coverages, +1 booster in yr 4-5
                            rtss_types = ['simple', 'booster1', 'booster2', 'booster3',
                                          'booster3']  # last deployment will get an property_restrictions 'GotBooster1'
                        else:
                            num_rtss_rounds = 3  # main EPI round and 3 boosters
                            rtss_types = ['simple', 'booster1', 'booster2']
                    elif 'campboost3' in rtss_mode:
                        num_rtss_rounds = 4  # main EPI round and 3 boosters
                        rtss_types = ['simple', 'booster1', 'booster2', 'booster3']
                    else:
                        num_rtss_rounds = 3  # main EPI round and 2 campaign boosters

                season_start_month = param_dic['season_start_month']
                rtss_campaign_booster = (round((
                                                       season_start_month - 1) * 30.4)) - 7 + 365  # 1 week before season start (adjusted for min age 24/36 months in set_up_interventions)

                ## seasonal campaign to be repeated every year ('campboost')
                rtss_booster_days = [rtss_campaign_booster + 365 * x for x in range(num_rtss_rounds - 1)]
                rtss_age_days = [rtss_age_days[
                                     0]] + rtss_booster_days  # initial EPI dose as well as dates of campaign booster doses

                try:
                    rtss_types
                except:
                    rtss_types = ['simple', 'booster1'] + ['booster2'] * (len(rtss_booster_days) - 1)

                effect_rtss = [0.8, 0.4] + [0.4] * (len(rtss_booster_days) - 1)
                deploy_types = ['EPI_cohort', 'campboost'] + ['campboost'] * (len(rtss_booster_days) - 1)
                tsteps_btwn_repetitions = [-1] * num_rtss_rounds

            elif 'noboost' in rtss_mode:
                num_rtss_rounds = 1
                rtss_age_days = rtss_age_days[0]
                rtss_types = ['simple']
                effect_rtss = [0.8]
                deploy_types = ['EPI_cohort']
                tsteps_btwn_repetitions = [-1]
            else:
                num_rtss_rounds = len(rtss_age_days)

                ## Single deployment at defined ages (simulation days)
                rtss_types = ['simple'] + ['booster'] * (len(rtss_age_days) - 1)
                effect_rtss = [0.8] + [0.4] * (len(rtss_age_days) - 1)
                deploy_types = ['EPI_cohort'] * (len(rtss_age_days))
                tsteps_btwn_repetitions = [-1] * (len(rtss_age_days))

            booster_coverage_rtss = param_dic['booster_coverage_rtss']
            max_age_rtss = param_dic['max_age_rtss']
            decay_time_constant_rtss = param_dic['decay_time_constant_rtss']
            decay_class_rtss = param_dic['decay_class_rtss']

            # iterate through RTS,S scenarios, creating csvs for each (one row for initial dose, one row for booster)
            rtss_coverage = param_dic['rtss_coverage']
            for rtss in rtss_coverage:
                if rtss == 0:
                    df = pd.DataFrame()
                else:
                    if rtss_mode == 'campboost_flexible':
                        coverage_levels = [rtss] + [0.5] * 3  # corresponding to [0.5, 0.25, 0.125]
                    elif rtss_mode == 'campboost2_flexible':
                        coverage_levels = [rtss] + [0.5] * 3 + [booster_coverage_rtss]
                    else:
                        if rtss == 'custom':
                            coverage_levels = [np.mean([0.673, 0.697, 0.493])] + [0.493]
                        else:
                            coverage_levels = [rtss] + [booster_coverage_rtss] * (num_rtss_rounds - 1)
                    df = pd.DataFrame(
                        {'coverage_levels': coverage_levels,
                         'round': list(range(1, num_rtss_rounds + 1)),
                         'rtss_types': rtss_types,
                         'initial_killing': effect_rtss,
                         'decay_time_constant': [decay_time_constant_rtss] * num_rtss_rounds,
                         'decay_class': [decay_class_rtss] * num_rtss_rounds,
                         'rtss_touchpoints': ['NA'] * num_rtss_rounds,
                         'RTSS_day': rtss_age_days,
                         'repetitions': [1] * num_rtss_rounds,
                         'tsteps_btwn_repetitions': tsteps_btwn_repetitions,
                         'deploy_type': deploy_types,
                         'agemin': [0] * num_rtss_rounds,
                         'agemax': [max_age_rtss] * num_rtss_rounds,
                         'distribution_name': ['CONSTANT_DISTRIBUTION'] * num_rtss_rounds,
                         'distribution_std': [1] * num_rtss_rounds,
                         'run_col': ['run'] * num_rtss_rounds})
                if rtss == 'custom':
                    rtss_csv = f'RTSS_{rtss_mode}_{rtss}coverage.csv'
                else:
                    rtss_csv = f'RTSS_{rtss_mode}_{int(rtss * 100)}coverage.csv'
                    """Add individual property restrictions (optional)"""
                    if rtss > 0 and 'flexible' in rtss_mode:
                        df['rtss_property_restrictions'] = ''
                        df.at[len(df) - 1, 'rtss_property_restrictions'] = 'GotBooster1'
                print(f'Writing {rtss_csv}')
                df.to_csv(os.path.join(base_scenario_filepath, 'RTSS', rtss_csv))

    if PMC == True:
        pmc_3tp = [int(round(x * (365 / 12), 0)) for x in [10 / 4, 14 / 4, 9]]
        pmc_touchpoints_dic = {
            '1tp': pmc_3tp[0],
            '3tp': pmc_3tp,
            '4tp': sorted(pmc_3tp + [int(round(6 * (365 / 12), 0))]),
            '4tp2ndyr': sorted(pmc_3tp + [int(round(12 * (365 / 12), 0))]),
            '5tp': sorted(pmc_3tp + [int(round(x * (365 / 12), 0)) for x in [6, 12]]),
            '5tp2ndyr': sorted(pmc_3tp + [int(round(x * (365 / 12), 0)) for x in [12, 15]]),
            '6tp2ndyr': sorted(pmc_3tp + [int(round(x * (365 / 12), 0)) for x in [6, 12, 15]]),
            '7tp2ndyr': sorted(pmc_3tp + [int(round(x * (365 / 12), 0)) for x in [6, 12, 15, 18]])
        }
        pmc_customcoverage_dic = {
            '3tp': [0.673, 0.697, 0.493],
            '4tp': [0.673, 0.697, 0.584, 0.493],
            '4tp2ndyr': [0.673, 0.697, 0.584, 0.493],
            '5tp': [0.673, 0.697, 0.584, 0.493, 0.493],
            '5tp2ndyr': [0.673, 0.697, 0.584, 0.493, 0.493],
            '6tp2ndyr': [0.673, 0.697, 0.584, 0.493, 0.493, 0.493],
            '7tp2ndyr': [0.673, 0.697, 0.584, 0.493, 0.493, 0.493, 0.493]
        }
        for i, pmc_mode in enumerate(param_dic['pmc_mode']):
            pmc_coverage = param_dic['pmc_coverage']
            pmc_3tp = [int(round(x * (365 / 12), 0)) for x in [10 / 4, 14 / 4, 9]]
            if pmc_mode in pmc_touchpoints_dic.keys():
                pmc_touchpoints = pmc_touchpoints_dic[pmc_mode]
            else:
                pmc_touchpoints = param_dic['pmc_touchpoints'][i]
            ntp = len(pmc_touchpoints)
            # EPI_cohort using campaign style deployment, distribution_name and distribution_std not used
            for pmc in pmc_coverage:
                if pmc == 0:
                    df = pd.DataFrame()
                else:
                    if pmc == 'custom':
                        coverage_levels = pmc_customcoverage_dic[pmc_mode]
                    else:
                        coverage_levels = [pmc] * ntp
                    df = pd.DataFrame({'round': range(1, ntp + 1),
                                       'coverage_levels': coverage_levels,
                                       'pmc_touchpoints': pmc_touchpoints,
                                       'PMC_day': pmc_touchpoints,
                                       'IPTi_day': pmc_touchpoints,
                                       'repetitions': [1] * ntp,
                                       'tsteps_btwn_repetitions': [-1] * ntp,
                                       'agemin': [0] * ntp,
                                       'agemax': [2] * ntp,
                                       'deploy_type': ['EPI_cohort'] * ntp,
                                       'drug_code': ['SP'] * ntp,
                                       'run_col': ['run'] * ntp})

                # fname = f'PMC_{len(pmc_touchpoints)}tp_constant_{int(100 * pmc)}coverage.csv'
                try:
                    fname = f'PMC_{pmc_mode}_{int(100 * pmc)}coverage.csv'
                except:
                    fname = f'PMC_{pmc_mode}_{pmc}coverage.csv'
                print(f'Writing {fname}')
                df.to_csv(os.path.join(base_scenario_filepath, 'PMC', fname),
                          index=False)


def combine_input_files(scen_csv1, scen_csv2, combined_scen_csv):
    scen_df1 = pd.read_csv(os.path.join(base_scenario_filepath, scen_csv1),
                           encoding='latin')
    scen_df2 = pd.read_csv(os.path.join(base_scenario_filepath, scen_csv2),
                           encoding='latin')
    scen_df = pd.concat([scen_df1, scen_df2]).reset_index(drop=True)

    # re-label settings
    scen_df['setting_id'] = ['HX' + str(ii) for ii in list(range(scen_df.shape[0]))]
    scen_df.set_index(keys='setting_id', inplace=True)

    scen_df.to_csv(os.path.join(base_scenario_filepath, combined_scen_csv),
                   encoding='latin')


def remove_duplicate_zerocoverage_scenarios(scen_csv, overwrite=False, return_csv=False):
    # remove rows that are effectively duplicates: these are created from taking the full factorial of parameter sets -
    #    when an intervention coverage is zero, it is not necessary to have separate simulations for different access
    #    groups or schedule types of that intervention
    scen_df = pd.read_csv(os.path.join(project_path, 'simulation_inputs/generic_cohort_ds', scen_csv),
                          encoding='latin')
    # remove duplicate RTS,S rows when RTS,S coverage is zero
    rtss0_subset_cols = ['annual_EIR', 'seasonality', 'cm_coverage', 'pmc_coverage', 'frac_high_access',
                         'cm_target_group', 'pmc_target_group', 'pmc_mode']  # columns to identify duplicates
    rtss0_subset_cols = list(set(rtss0_subset_cols).intersection(scen_df.columns))
    scen_df_rtss0 = scen_df[scen_df.rtss_coverage == 0]
    scen_df_rtss_remainder = scen_df[scen_df.rtss_coverage != 0]
    scen_df_rtss0_fixed = scen_df_rtss0[~scen_df_rtss0.duplicated(subset=rtss0_subset_cols, keep='first')]
    scen_df = pd.concat([scen_df_rtss0_fixed, scen_df_rtss_remainder])

    # remove duplicate PMC rows when PMC coverage is zero
    pmc0_subset_cols = ['annual_EIR', 'seasonality', 'cm_coverage',
                        'rtss_coverage', 'frac_high_access', 'rtss_mode', 'minBoostAge', 'rtss_booster1_min_age',
                        'cm_target_group', 'rtss_target_group']  # columns to identify duplicates
    pmc0_subset_cols = list(set(pmc0_subset_cols).intersection(scen_df.columns))
    scen_df_pmc0 = scen_df[scen_df.pmc_coverage == 0]
    scen_df_pmc_remainder = scen_df[scen_df.pmc_coverage != 0]
    scen_df_pmc0_fixed = scen_df_pmc0[~scen_df_pmc0.duplicated(subset=pmc0_subset_cols, keep='first')]
    scen_df = pd.concat([scen_df_pmc0_fixed, scen_df_pmc_remainder])

    # re-sort rows
    scen_df = scen_df.sort_values(
        by=['annual_EIR', 'seasonality', 'cm_coverage', 'pmc_coverage', 'rtss_coverage', 'rtss_target_group'])

    # re-label settings
    scen_df['setting_id'] = ['HX' + str(ii) for ii in list(range(scen_df.shape[0]))]
    scen_df.set_index(keys='setting_id', inplace=True)

    if not overwrite:
        scen_csv = scen_csv.replace('.csv', '_cleaned.csv')
    scen_df.to_csv(os.path.join(project_path, 'simulation_inputs/generic_cohort_ds', scen_csv), encoding='latin')

    if return_csv:
        print(f'Returning dataframe for {scen_csv}')
        return scen_df


if __name__ == "__main__":
    home_path, data_path, git_dir, project_path = load_box_paths()
    base_scenario_filepath = os.path.join(project_path, 'simulation_inputs/generic_cohort_ds')

    """
    View parameter defaults and options
    """
    p = get_parameter_space()
    print(p)
    n_scen = []

    # 1tp: 10 weeks
    # 2tp:     14 weeks, 9 moths
    # 3tp: 10, 14 weeks, 9 moths
    # 4tp: + 6
    # 4tp2ndyr: + 12
    # 5tp: + 6,12
    # 5tp2ndyr: + 12, 15
    # 6tp2ndyr: + 6, 12, 15
    # 7tp2ndyr: + 6, 12, 15, 1

    param_dic = get_parameter_space()
    param_dic.update({'cm_coverage': np.linspace(0, 1, 11)})
    param_dic.update({'pmc_coverage': np.linspace(0, 1, 11)})
    param_dic.update({'rtss_coverage': np.linspace(0, 1, 11)})
    param_dic.update({'pmc_mode': ['1tp', '3tp', '4tp', '4tp2ndyr', '5tp', '5tp2ndyr', '6tp2ndyr', '7tp2ndyr']})
    intervention_inputs(param_dic, PMC=True, RTSS=False, CM=False)

    param_dic = get_parameter_space()
    param_dic.update({'cm_coverage': np.linspace(0, 1, 9)})
    param_dic.update({'pmc_coverage': np.linspace(0, 1, 9)})
    param_dic.update({'rtss_coverage': np.linspace(0, 1, 9)})
    param_dic.update({'pmc_mode': ['1tp', '3tp', '4tp', '4tp2ndyr', '5tp', '5tp2ndyr', '6tp2ndyr', '7tp2ndyr']})
    intervention_inputs(param_dic, PMC=True, RTSS=True, CM=True)

    """
    0 - single PMC and/or RTS,S 
    """
    param_dic = get_parameter_space()
    param_dic.update({'pmc_coverage': [0, 1]})
    param_dic.update({'cm_coverage': [0.4, 0.6]})
    param_dic.update({'rtss_coverage': [0]})
    param_dic.update({'pmc_mode': ['1tp']})
    param_dic.update({'pmc_touchpoints': [76]})
    param_dic.update({'annual_EIR': [1, 4, 8, 16, 32, 64, 128]})
    param_dic.update({'seasonality': ['constant']})  # ,'season2','season3'
    create_scenarios_mastercsv(param_dic, fname_out='generic_single_PMC')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_single_PMC.csv', overwrite=True, return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    param_dic = get_parameter_space()
    param_dic.update({'pmc_coverage': [0]})
    param_dic.update({'rtss_coverage': [0, 0.2, 0.4, 0.6, 0.8, 1]})
    param_dic.update({'annual_EIR': [32]})
    param_dic.update({'seasonality': ['constant']})  # ,'season2','season3'
    create_scenarios_mastercsv(param_dic, fname_out='generic_single_RTSS')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_single_RTSS.csv', overwrite=True, return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    """
    1 - PMC modes and/or RTS,S (main, Fig 2, 3)
    """
    param_dic = get_parameter_space()
    param_dic.update({'pmc_coverage': [0, 0.2, 0.4, 0.6, 0.8, 1]})
    param_dic.update({'rtss_coverage': [0, 0.8]})
    param_dic.update({'pmc_mode': ['3tp', '4tp', '4tp2ndyr', '5tp', '5tp2ndyr', '6tp2ndyr', '7tp2ndyr']})
    param_dic.update({'annual_EIR': [32]})
    param_dic.update({'seasonality': ['season1']})  # ,'season2','season3'
    create_scenarios_mastercsv(param_dic, fname_out='generic_PMCmode_RTSS')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_PMCmode_RTSS.csv', overwrite=True, return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    param_dic = get_parameter_space()
    param_dic.update({'pmc_coverage': [0, 0.25, 0.5, 0.75, 1]})
    param_dic.update({'rtss_coverage': [0, 0.25, 0.5, 0.75, 1]})
    param_dic.update({'pmc_mode': ['3tp', '5tp', '7tp2ndyr']})
    param_dic.update({'annual_EIR': [32]})
    param_dic.update({'seasonality': ['season1']})  # ,'season2','season3'
    create_scenarios_mastercsv(param_dic, fname_out='generic_PMCmode_RTSS_cov')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_PMCmode_RTSS_cov.csv', overwrite=True,
                                                 return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    """
    2 - PMC and RTS,S by EIR sweep   (main, Fig 4)
    """
    param_dic = get_parameter_space()
    param_dic.update({'pmc_coverage': [0, 0.8]})
    param_dic.update({'rtss_coverage': [0, 0.8]})
    param_dic.update({'pmc_mode': ['3tp', '5tp', '7tp2ndyr']})
    param_dic.update({'annual_EIR': [1, 4, 8, 16, 32, 64, 128]})
    param_dic.update({'seasonality': ['season1']})  # 'constant','season2','season3'
    create_scenarios_mastercsv(param_dic, fname_out='generic_PMC_RTSS_EIR')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_PMC_RTSS_EIR.csv', overwrite=True, return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    ## PMC only
    param_dic = get_parameter_space()
    param_dic.update({'pmc_coverage': [0, 0.8]})
    param_dic.update({'rtss_coverage': [0]})
    param_dic.update({'pmc_mode': ['3tp', '4tp', '4tp2ndyr', '5tp', '5tp2ndyr', '6tp2ndyr', '7tp2ndyr']})
    param_dic.update({'annual_EIR': [1, 4, 16, 32, 64, 128]})
    param_dic.update({'seasonality': ['season1']})  # 'constant','season2','season3'
    create_scenarios_mastercsv(param_dic, fname_out='generic_PMCmode_EIR')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_PMCmode_EIR.csv', overwrite=True, return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    """
    2b - coverage - EIR _constant
    """
    param_dic = get_parameter_space()
    param_dic.update({'cm_coverage': [0.2, 0.4, 0.6]})
    param_dic.update({'pmc_coverage': [0, 0.2, 0.4, 0.6, 0.8, 1]})
    param_dic.update({'rtss_coverage': [0]})
    param_dic.update({'pmc_mode': ['3tp']})
    param_dic.update({'annual_EIR': [1, 4, 16, 32, 64, 128]})
    param_dic.update({'seasonality': ['constant']})  # 'constant','season2','season3'
    create_scenarios_mastercsv(param_dic, fname_out='generic_PMCcov_EIR_PMC3_constant')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_PMCcov_EIR_PMC3_constant.csv', overwrite=True,
                                                 return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    param_dic = get_parameter_space()
    param_dic.update({'cm_coverage': [0.2, 0.4, 0.6]})
    param_dic.update({'pmc_coverage': [0, 0.2, 0.4, 0.6, 0.8, 1]})
    param_dic.update({'rtss_coverage': [0]})
    param_dic.update({'pmc_mode': ['5tp2ndyr']})
    param_dic.update({'annual_EIR': [1, 4, 16, 32, 64, 128]})
    param_dic.update({'seasonality': ['constant']})  # 'constant','season2','season3'
    create_scenarios_mastercsv(param_dic, fname_out='generic_PMCcov_EIR_PMC5tp2ndyr_constant')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_PMCcov_EIR_PMC5tp2ndyr_constant.csv', overwrite=True,
                                                 return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    param_dic = get_parameter_space()
    param_dic.update({'cm_coverage': [0.2, 0.4, 0.6]})
    param_dic.update({'pmc_coverage': [0, 0.2, 0.4, 0.6, 0.8, 1]})
    param_dic.update({'rtss_coverage': [0]})
    param_dic.update({'pmc_mode': ['7tp2ndyr']})
    param_dic.update({'annual_EIR': [1, 4, 16, 32, 64, 128]})
    param_dic.update({'seasonality': ['constant']})  # 'constant','season2','season3'
    create_scenarios_mastercsv(param_dic, fname_out='generic_PMCcov_EIR_PMC7_constant')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_PMCcov_EIR_PMC7_constant.csv', overwrite=True,
                                                 return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    param_dic = get_parameter_space()
    param_dic.update({'cm_coverage': [0.2, 0.4, 0.6]})
    param_dic.update({'pmc_coverage': [0]})
    param_dic.update({'rtss_coverage': [0, 0.2, 0.4, 0.6, 0.8, 1]})
    param_dic.update({'pmc_mode': ['7tp2ndyr']})
    param_dic.update({'annual_EIR': [1, 4, 16, 32, 64, 128]})
    param_dic.update({'seasonality': ['constant']})  # 'constant','season2','season3'
    create_scenarios_mastercsv(param_dic, fname_out='generic_RTSScov_EIR_constant')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_RTSScov_EIR_constant.csv', overwrite=True,
                                                 return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    param_dic = get_parameter_space()
    param_dic.update({'cm_coverage': [0.2, 0.4, 0.6]})
    param_dic.update({'pmc_coverage': [0]})
    param_dic.update({'rtss_coverage': [0, 0.2, 0.4, 0.6, 0.8, 1]})
    param_dic.update({'pmc_mode': ['3tp']})
    param_dic.update({'annual_EIR': [1, 4, 16, 32, 64, 128]})
    param_dic.update({'seasonality': ['constant']})  # 'constant','season2','season3'
    create_scenarios_mastercsv(param_dic, fname_out='generic_PMC3_RTSScov_EIR_constant')

    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_PMC3_RTSScov_EIR_constant.csv', overwrite=True,
                                                 return_csv=True)
    df['pmc_coverage'] = df['rtss_coverage']
    df['PMC_filename'] = [x.replace('RTSS_constant_', 'PMC_3tp_') for x in df['RTSS_filename']]
    df['PMC_filename'] = [x.replace('RTSS', 'PMC') for x in df['PMC_filename']]
    df.to_csv(
        os.path.join(project_path, 'simulation_inputs/generic_cohort_ds', 'generic_PMC3_RTSScov_EIR_constant.csv'))
    n_scen = n_scen + [len(df)]
    print([len(df)])

    """
    2c - operational coverage - EIR _constant 
    """
    param_dic = get_parameter_space()
    param_dic.update({'pmc_coverage': ['custom']})
    param_dic.update({'rtss_coverage': ['custom']})
    param_dic.update({'pmc_mode': ['3tp', '4tp', '4tp2ndyr', '5tp', '5tp2ndyr', '6tp2ndyr', '7tp2ndyr']})
    intervention_inputs(param_dic, PMC=True, RTSS=True)

    param_dic = get_parameter_space()
    param_dic.update({'cm_coverage': [0.2, 0.4, 0.6]})
    param_dic.update({'pmc_coverage': [0, 'custom']})
    param_dic.update({'rtss_coverage': [0]})
    param_dic.update({'pmc_mode': ['3tp']})
    param_dic.update({'annual_EIR': [1, 4, 16, 32, 64, 128]})
    param_dic.update({'seasonality': ['constant']})  # 'constant','season2','season3'
    create_scenarios_mastercsv(param_dic, fname_out='generic_PMCcov_EIR_PMC3_constant_operational')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_PMCcov_EIR_PMC3_constant_operational.csv',
                                                 overwrite=True,
                                                 return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    param_dic = get_parameter_space()
    param_dic.update({'cm_coverage': [0.2, 0.4, 0.6]})
    param_dic.update({'pmc_coverage': [0, 'custom']})
    param_dic.update({'rtss_coverage': [0]})
    param_dic.update({'pmc_mode': ['5tp2ndyr']})
    param_dic.update({'annual_EIR': [1, 4, 16, 32, 64, 128]})
    param_dic.update({'seasonality': ['constant']})  # 'constant','season2','season3'
    create_scenarios_mastercsv(param_dic, fname_out='generic_PMCcov_EIR_PMC5tp2ndyr_constant_operational')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_PMCcov_EIR_PMC5tp2ndyr_constant_operational.csv',
                                                 overwrite=True,
                                                 return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    param_dic = get_parameter_space()
    param_dic.update({'cm_coverage': [0.2, 0.4, 0.6]})
    param_dic.update({'pmc_coverage': [0, 'custom']})
    param_dic.update({'rtss_coverage': [0]})
    param_dic.update({'pmc_mode': ['7tp2ndyr']})
    param_dic.update({'annual_EIR': [1, 4, 16, 32, 64, 128]})
    param_dic.update({'seasonality': ['constant']})  # 'constant','season2','season3'
    create_scenarios_mastercsv(param_dic, fname_out='generic_PMCcov_EIR_PMC7_constant_operational')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_PMCcov_EIR_PMC7_constant_operational.csv',
                                                 overwrite=True,
                                                 return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    param_dic = get_parameter_space()
    param_dic.update({'cm_coverage': [0.2, 0.4, 0.6]})
    param_dic.update({'pmc_coverage': [0]})
    param_dic.update({'rtss_coverage': [0, 'custom']})
    param_dic.update({'pmc_mode': ['7tp2ndyr']})
    param_dic.update({'annual_EIR': [1, 4, 16, 32, 64, 128]})
    param_dic.update({'seasonality': ['constant']})  # 'constant','season2','season3'
    create_scenarios_mastercsv(param_dic, fname_out='generic_RTSScov_EIR_constant_operational')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_RTSScov_EIR_constant_operational.csv',
                                                 overwrite=True,
                                                 return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    param_dic = get_parameter_space()
    param_dic.update({'cm_coverage': [0.2, 0.4, 0.6]})
    param_dic.update({'pmc_coverage': [0]})
    param_dic.update({'rtss_coverage': [0, 'custom']})
    param_dic.update({'pmc_mode': ['3tp']})
    param_dic.update({'annual_EIR': [1, 4, 16, 32, 64, 128]})
    param_dic.update({'seasonality': ['constant']})  # 'constant','season2','season3'
    create_scenarios_mastercsv(param_dic, fname_out='generic_PMC3_RTSScov_EIR_constant_operational')

    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_PMC3_RTSScov_EIR_constant_operational.csv',
                                                 overwrite=True,
                                                 return_csv=True)
    df['pmc_coverage'] = df['rtss_coverage']
    df['PMC_filename'] = [x.replace('RTSS_constant_', 'PMC_3tp_') for x in df['RTSS_filename']]
    df['PMC_filename'] = [x.replace('RTSS', 'PMC') for x in df['PMC_filename']]
    df.to_csv(
        os.path.join(project_path, 'simulation_inputs/generic_cohort_ds',
                     'generic_PMC3_RTSScov_EIR_constant_operational.csv'))
    n_scen = n_scen + [len(df)]
    print([len(df)])

    """
    3 - Standard PMC and RTS,S by CM sweep  (supplement)
    """
    param_dic = get_parameter_space()
    param_dic.update({'pmc_coverage': [0, 0.2, 0.4, 0.6, 0.8]})
    param_dic.update({'rtss_coverage': [0]})
    param_dic.update({'pmc_mode': ['3tp']})
    param_dic.update({'annual_EIR': [4, 32, 64, 128]})
    param_dic.update({'cm_coverage': [0.2, 0.4, 0.6, 0.8]})
    param_dic.update({'seasonality': ['season1']})
    create_scenarios_mastercsv(param_dic, fname_out='generic_PMC_CM')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_PMC_CM.csv', overwrite=True, return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    """
    4 - High low access group
    """
    param_dic = get_parameter_space()
    param_dic.update({'pmc_coverage': [0.2, 0.4, 0.6, 0.8]})
    param_dic.update({'rtss_coverage': [0]})
    param_dic.update({'pmc_mode': ['3tp', '7tp2ndyr']})
    param_dic.update({'annual_EIR': [32]})
    param_dic.update({'cm_coverage': [0.6]})
    param_dic.update({'seasonality': ['season1']})
    param_dic.update({'cm_target_group': ['high', 'random', 'low']})
    param_dic.update({'pmc_target_group': ['high', 'random', 'low']})
    create_scenarios_mastercsv(param_dic, high_access_frac='pmc_coverage', fname_out='generic_PMC_CM_accesscor')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_PMC_CM_accesscor.csv', overwrite=True,
                                                 return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    print(f'Number of scenarios to run: {n_scen}')
    nseeds = 5
    print(f'Number of scenarios to run including {nseeds} seeds and 12 seasonal_birth_months :'
          f' {[x * nseeds * 12 for x in n_scen]}')

    """
    5 - Additional custom scenarios (supplement)
    """

    """
    5a - Addition exploration role of 10 weeks
    """
    param_dic = get_parameter_space()
    param_dic.update({'pmc_coverage': [0, 0.8]})
    param_dic.update({'rtss_coverage': [0]})
    param_dic.update({'pmc_mode': ['3tp', '2tp', '3atp', '3btp', '3ctp']})
    param_dic.update(
        {'pmc_touchpoints': [[76, 106, 274], [106, 274], [106, 182, 274], [106, 274, 365], [106, 274, 456]]})
    param_dic.update({'annual_EIR': [1, 4, 8, 16, 32, 64, 128]})
    param_dic.update({'seasonality': ['season1']})  # ,'season2','season3'
    intervention_inputs(param_dic, PMC=True)
    create_scenarios_mastercsv(param_dic, fname_out='generic_PMC3mode_redistribution')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_PMC3mode_redistribution.csv', overwrite=True,
                                                 return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    """
    5b - Addition exploration start at 4 vs 6 months
    """
    param_dic = get_parameter_space()
    param_dic.update({'pmc_coverage': [0, 0.8]})
    param_dic.update({'rtss_coverage': [0]})
    param_dic.update({'pmc_mode': ['3tp', '3tp4mth', '3tp6mth']})
    param_dic.update(
        {'pmc_touchpoints': [[76, 106, 274], [122, 274, 365], [182, 274, 365]]})
    param_dic.update({'annual_EIR': [1, 4, 8, 16, 32, 64, 128]})
    param_dic.update({'seasonality': ['season1']})  # ,'season2','season3'
    intervention_inputs(param_dic, PMC=True)
    create_scenarios_mastercsv(param_dic, fname_out='generic_PMC3mode_shifted')
    df = remove_duplicate_zerocoverage_scenarios(scen_csv='generic_PMC3mode_shifted.csv', overwrite=True,
                                                 return_csv=True)
    n_scen = n_scen + [len(df)]
    print([len(df)])

    """
    6 - Operational coverage
    """
