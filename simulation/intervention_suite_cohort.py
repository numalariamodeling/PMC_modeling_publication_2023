import pandas as pd
import numpy as np
import math
from dtk.interventions.property_change import change_individual_property
from malaria.interventions.adherent_drug import configure_adherent_drug
from malaria.interventions.health_seeking import add_health_seeking
from malaria.interventions.malaria_drug_campaigns import add_drug_campaign
from malaria.interventions.malaria_vaccine import add_vaccine
from malaria.reports.MalariaReport import add_event_counter_report


def calc_high_low_access_coverages(coverage_all, high_access_frac):
    if (high_access_frac < 1) & (coverage_all >= high_access_frac):
        coverage_high = 1
        coverage_low = (coverage_all - high_access_frac) / (1 - high_access_frac)
    else:
        coverage_high = coverage_all / high_access_frac
        coverage_low = 0
    return [coverage_high, coverage_low]


def set_IIV_vacc(eff_lower, eff_upper, eff_SD, params):
    decay_params_IIV = copy.deepcopy(params)
    cur_eff = params['Waning_Config']['Initial_Effect']
    new_eff = stats.truncnorm.rvs((eff_lower - cur_eff) / eff_SD, (eff_upper - cur_eff) / eff_SD, loc=cur_eff,
                                  scale=eff_SD)
    decay_params_IIV['Waning_Config']['Initial_Effect'] = new_eff
    return decay_params_IIV


class InterventionSuite:
    # hs
    hs_ds_col = 'repDS'
    hs_duration = None

    hs_start_col = 'simday'
    hs_coverage_age = {  # column: [agemin, agemax]
        'U5_coverage': [0, 5],
        'adult_coverage': [5, 100]
    }
    hs_severe_coverage_age = {
        'severe_cases': [0, 100]
    }
    hs_rates = 0.3
    hs_severe_rates = 0.5

    # rtss
    rtss_ds_col = 'DS_Name'
    rtss_type_col = 'rtss_types'
    rtss_start_col = 'RTSS_day'
    rtss_coverage_col = 'coverage_levels'
    rtss_touchpoint_col = 'rtss_touchpoints'  # days since births!
    rtss_deploy_type_col = 'deploy_type'
    rtss_distribution_col = 'distribution_name'
    rtss_std_col = 'distribution_std'
    rtss_min_age_col = 'agemin'
    rtss_max_age_col = 'agemax'
    rtss_repetitions = 'repetitions'
    rtss_tsteps_btwn_repetitions = 'tsteps_btwn_repetitions'
    rtss_init_eff_col = 'initial_killing'
    rtss_decay_t_col = 'decay_time_constant'
    rtss_decay_class_col = 'decay_class'
    rtss_property_restrictions = 'rtss_property_restrictions'

    # pmc
    pmc_ds_col = 'DS_Name'
    pmc_start_col = 'PMC_day'
    pmc_coverage_col = 'coverage_levels'
    pmc_touchpoint_col = 'pmc_touchpoints'  # days since birth
    pmc_distribution_col = 'distribution_name'
    pmc_std_col = 'distribution_std'
    pmc_agemin = 'agemin'
    pmc_agemax = 'agemax'
    pmc_repetitions = 'repetitions'
    pmc_tsteps_btwn_repetitions = 'tsteps_btwn_repetitions'

    def add_ds_hs(self, cb, hs_df, my_ds, high_access_ip_frac=0, cm_target_group='random'):
        ds_col = self.hs_ds_col
        duration = self.hs_duration
        df = hs_df[hs_df[ds_col] == my_ds]
        for r, row in df.iterrows():
            self.add_hs_from_file(cb, row, duration=duration, high_access_ip_frac=high_access_ip_frac,
                                  cm_target_group=cm_target_group)

        return len(df)

    def add_hs_from_file(self, cb, row, duration, high_access_ip_frac=0, cm_target_group='random'):
        rates = self.hs_rates
        severe_rates = self.hs_severe_rates

        start_day = row[self.hs_start_col]  # if start_day_override < 0 else start_day_override
        if duration is None:
            duration = row['duration']

        # Uncomplicated - access depends on IP
        if 'custom_age_coverage' in row.index:
            self.hs_coverage_age = {'custom_age_coverage': [int(row['agemin']), int(row['agemax'])]}
            self.hs_severe_coverage_age = {'severe_cases': [int(row['agemin']), int(row['agemax'])]}

        if (high_access_ip_frac > 0) and (cm_target_group != 'random'):
            for key, value in self.hs_coverage_age.items():
                if cm_target_group == 'low':  # CM is preferentially given to the 'low-access' group
                    high_low_coverages = calc_high_low_access_coverages(coverage_all=row[key],
                                                                        high_access_frac=1 - high_access_ip_frac)
                    high_low_coverages = [high_low_coverages[1], high_low_coverages[0]]
                elif cm_target_group == 'high':
                    high_low_coverages = calc_high_low_access_coverages(coverage_all=row[key],
                                                                        high_access_frac=high_access_ip_frac)
                else:
                    print('WARNING: name for CM access-group targeting not recognized.')

                # # high-access coverage
                targets = [{'trigger': 'NewClinicalCase',
                            'coverage': high_low_coverages[0],
                            'agemin': value[0],
                            'agemax': value[1],
                            'seek': 1,
                            'rate': rates
                            }]
                add_health_seeking(cb, start_day=start_day,
                                   targets=targets,
                                   drug=['Artemether', 'Lumefantrine'], duration=duration,
                                   ind_property_restrictions=[{'AccessToInterventions': 'higher'}])
                # low-access coverage
                targets = [{'trigger': 'NewClinicalCase',
                            'coverage': high_low_coverages[1],
                            'agemin': value[0],
                            'agemax': value[1],
                            'seek': 1,
                            'rate': rates
                            }]
                add_health_seeking(cb, start_day=start_day,
                                   targets=targets,
                                   drug=['Artemether', 'Lumefantrine'], duration=duration,
                                   ind_property_restrictions=[{'AccessToInterventions': 'lower'}])
        else:
            targets = []
            for key, value in self.hs_coverage_age.items():
                targets.append({
                    'trigger': 'NewClinicalCase',
                    'coverage': row[key],
                    'agemin': value[0],
                    'agemax': value[1],
                    'seek': 1,
                    'rate': rates
                })
            add_health_seeking(cb, start_day=start_day,
                               targets=targets,
                               drug=['Artemether', 'Lumefantrine'], duration=duration)

        # Severe - same coverage for all access-IPS
        targets = []
        for key, value in self.hs_severe_coverage_age.items():
            targets.append({
                'trigger': 'NewSevereCase',
                'coverage': row[key],
                'agemin': value[0],
                'agemax': value[1],
                'seek': 1,
                'rate': severe_rates
            })
        add_health_seeking(cb, start_day=start_day,
                           targets=targets,
                           drug=['Artemether', 'Lumefantrine'], duration=duration,
                           broadcast_event_name='Received_Severe_Treatment')

    def change_rtss_ips(self, cb, change_booster_IPs=False):
        change_individual_property(cb,
                                   target_property_name='VaccineStatus',
                                   target_property_value='GotVaccine',
                                   ind_property_restrictions=[{'VaccineStatus': 'None'}],
                                   trigger_condition_list=['Received_Vaccine'],
                                   blackout_flag=False)
        # UPDATE 2021-09-15: Currently, the generic model is set up to give Booster1 in the first year and Booster2 in the second year (this differs from the Nigeria setup).
        #     Booster2 is given to anyone who received the first vaccine, regardless of whether they also received Booster1.
        if change_booster_IPs:
            change_individual_property(cb,
                                       target_property_name='VaccineStatus',
                                       target_property_value='GotBooster1',
                                       ind_property_restrictions=[{'VaccineStatus': 'GotVaccine'}],
                                       trigger_condition_list=['Received_Vaccine'],
                                       blackout_flag=False)
            change_individual_property(cb,
                                       target_property_name='VaccineStatus',
                                       target_property_value='GotBooster2',
                                       ind_property_restrictions=[{'VaccineStatus': 'GotBooster1'}],
                                       trigger_condition_list=['Received_Vaccine'],
                                       blackout_flag=False)

    def add_ds_rtss(self, cb, rtss_df, my_ds, high_access_ip_frac=0, rtss_target_group='random',
                    rtss_booster1_min_age=730, cohort_month_shift=0):
        rtss_df = rtss_df[rtss_df[self.rtss_ds_col].str.upper() == my_ds.upper()]
        change_booster_IPs = False  # updated in booster statements if applicable

        rtss_booster2_min_age = rtss_booster1_min_age + 365
        rtss_booster3_min_age = rtss_booster2_min_age + 365

        if 'rtss_property_restrictions' in rtss_df.columns:
            change_booster_IPs = True  # change booster IPs when specifying coverage per booster dose
        """Use campaign-style deployment targeted to specific ages (since births disabled in simulation)"""
        for r, row in rtss_df.iterrows():
            try:
                Waning_Class = rtss_df['rtss_decay_class_col'].unique()[0]
            except:
                Waning_Class = "WaningEffectExponential"

            vaccine_params = {"Waning_Config": {"Initial_Effect": row[self.rtss_init_eff_col],
                                                "Decay_Time_Constant": row[self.rtss_decay_t_col],
                                                "class": Waning_Class}}

            vtype = row[self.rtss_type_col]
            if vtype == 'booster' or vtype == 'booster1':
                """everyone who received the original vaccine has the same probability of getting a booster, regardless of IP 
                (otherwise, there will be different booster coverages depending on the fraction of people in the high IP)
                if campboost, receive booster between 24-35 month
                """
                if row[self.rtss_deploy_type_col] == 'EPI_cohort':
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[row[self.rtss_start_col]],
                                coverage=row[self.rtss_coverage_col],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                ind_property_restrictions=[{'VaccineStatus': 'GotVaccine'}],
                                target_group={'agemin': row[self.rtss_min_age_col],
                                              'agemax': row[self.rtss_max_age_col]})
                elif row[self.rtss_deploy_type_col] == 'campboost':
                    # calculate RTS,S booster day, given cohort month shift
                    start_day0 = row[self.rtss_start_col]
                    start_day = start_day0 - round(30.4 * cohort_month_shift)
                    # if booster would occur before the eligible age, wait until the next year
                    while start_day < rtss_booster1_min_age:
                        start_day = start_day + 365
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[start_day],
                                coverage=row[self.rtss_coverage_col],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                ind_property_restrictions=[{'VaccineStatus': 'GotVaccine'}],
                                target_group={'agemin': row[self.rtss_min_age_col],
                                              'agemax': row[self.rtss_max_age_col]})

            elif vtype == 'booster2':
                """everyone who received the original vaccine has the same probability of getting a booster, regardless of IP 
                (otherwise, there will be different booster coverages depending on the fraction of people in the high IP)
                if campboost, receive booster between 36-47 month
                """
                try:
                    booster_restr = row[self.rtss_property_restrictions]
                    if booster_restr == 'GotBooster1':
                        change_booster_IPs = True
                    if len(booster_restr) < 1:
                        booster_restr = 'GotVaccine'
                except:
                    booster_restr = 'GotVaccine'
                if row[self.rtss_deploy_type_col] == 'EPI_cohort':
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[row[self.rtss_start_col]],
                                coverage=row[self.rtss_coverage_col],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                ind_property_restrictions=[{'VaccineStatus': booster_restr}],
                                # even if someone didn't get booster1 during the 24-month EPI visit, they can still get a booster during an EPI visit the following year
                                target_group={'agemin': row[self.rtss_min_age_col],
                                              'agemax': row[self.rtss_max_age_col]})
                elif row[self.rtss_deploy_type_col] == 'campboost':
                    # calculate RTS,S booster day, given cohort month shift
                    start_day0 = row[self.rtss_start_col]
                    start_day = start_day0 - round(30.4 * cohort_month_shift)
                    # if booster would occur before the eligible age, wait until the next year
                    while start_day < rtss_booster2_min_age:
                        start_day = start_day + 365
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[start_day],
                                coverage=row[self.rtss_coverage_col],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                ind_property_restrictions=[{'VaccineStatus': booster_restr}],
                                # even if someone didn't get booster1 during the first campaign, they can still get a booster during the following campaign
                                target_group={'agemin': row[self.rtss_min_age_col],
                                              'agemax': row[self.rtss_max_age_col]})
            elif vtype == 'booster3':
                """ if campboost, receive booster between 48-59 month
                """
                try:
                    booster_restr = row[self.rtss_property_restrictions]
                    if len(booster_restr) < 1:
                        booster_restr = 'GotVaccine'
                except:
                    booster_restr = 'GotVaccine'
                if row[self.rtss_deploy_type_col] == 'EPI_cohort':
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[row[self.rtss_start_col]],
                                coverage=row[self.rtss_coverage_col],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                ind_property_restrictions=[{'VaccineStatus': booster_restr}],
                                target_group={'agemin': row[self.rtss_min_age_col],
                                              'agemax': row[self.rtss_max_age_col]})
                elif row[self.rtss_deploy_type_col] == 'campboost':
                    # calculate RTS,S booster day, given cohort month shift
                    start_day0 = row[self.rtss_start_col]
                    start_day = start_day0 - round(30.4 * cohort_month_shift)
                    # if booster would occur before the eligible age, wait until the next year
                    while start_day < rtss_booster3_min_age:
                        start_day = start_day + 365
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[start_day],
                                coverage=row[self.rtss_coverage_col],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                ind_property_restrictions=[{'VaccineStatus': booster_restr}],
                                # even if someone didn't get booster1 during the previous campaigns, they can still get a booster during the following campaign
                                target_group={'agemin': row[self.rtss_min_age_col],
                                              'agemax': row[self.rtss_max_age_col]})
            else:
                if high_access_ip_frac > 0 and rtss_target_group != 'random':  # different coverages in different access groups
                    if rtss_target_group == 'low':
                        # RTSS is preferentially given to the 'low-access' group
                        high_low_coverages = calc_high_low_access_coverages(coverage_all=row[self.rtss_coverage_col],
                                                                            high_access_frac=1 - high_access_ip_frac)
                        high_low_coverages = [high_low_coverages[1], high_low_coverages[0]]
                    elif rtss_target_group == 'high':
                        high_low_coverages = calc_high_low_access_coverages(coverage_all=row[self.rtss_coverage_col],
                                                                            high_access_frac=high_access_ip_frac)
                    else:
                        print('WARNING: name for RTS,S access-group targeting not recognized.')
                    # high-access coverage
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[row[self.rtss_start_col]],
                                coverage=high_low_coverages[0],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                target_group={'agemin': row[self.rtss_min_age_col],
                                              'agemax': row[self.rtss_max_age_col]},
                                ind_property_restrictions=[{'AccessToInterventions': 'higher'}])

                    # low-access coverage
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[row[self.rtss_start_col]],
                                coverage=high_low_coverages[1],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                target_group={'agemin': row[self.rtss_min_age_col],
                                              'agemax': row[self.rtss_max_age_col]},
                                ind_property_restrictions=[{'AccessToInterventions': 'lower'}])

                else:  # uniform probability of getting the vaccine
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[row[self.rtss_start_col]],
                                coverage=row[self.rtss_coverage_col],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                target_group={'agemin': row[self.rtss_min_age_col],
                                              'agemax': row[self.rtss_max_age_col]})

            self.change_rtss_ips(cb, change_booster_IPs=change_booster_IPs)
            cb.update_params({
                "Report_Event_Recorder_Events": ['Births', 'PropertyChange', 'Received_Vaccine', 'Received_Treatment']
            })

        return len(rtss_df)

    def add_ds_pmc(self, cb, pmc_df, my_ds, high_access_ip_frac=0, pmc_target_group='random', vaccSP_offset=10):
        pmc_df = pmc_df[pmc_df[self.pmc_ds_col].str.upper() == my_ds.upper()]

        try:
            num_IIV_groups = pmc_df['num_IIV_groups'].unique()[0]  ## add to pmc_df in simulation script for using IIV
            pmc_IIV = True
        except:
            num_IIV_groups = None
            pmc_IIV = False

        #  used if  drug_code == 'vaccSP'
        decay_dur = 10
        decay_params = {"Waning_Config": {"Initial_Effect": 0.8,
                                          "Box_Duration": 32,
                                          "Decay_Time_Constant": (decay_dur) / math.log(2),
                                          "class": 'WaningEffectBoxExponential'}}

        """Use campaign-style deployment"""
        if not pmc_IIV:
            for r, row in pmc_df.iterrows():

                if high_access_ip_frac > 0 and pmc_target_group != 'random':  # different coverages in different access groups
                    if pmc_target_group == 'low':
                        # PMC is preferentially given to the 'low-access' group
                        high_low_coverages = calc_high_low_access_coverages(
                            coverage_all=row[self.pmc_coverage_col],
                            high_access_frac=1 - high_access_ip_frac)
                        high_low_coverages = [high_low_coverages[1], high_low_coverages[0]]
                    elif pmc_target_group == 'high':
                        high_low_coverages = calc_high_low_access_coverages(
                            coverage_all=row[self.pmc_coverage_col],
                            high_access_frac=high_access_ip_frac)
                    else:
                        print('WARNING: name for PMC access-group targeting not recognized.')

                    # higher access
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=decay_params,
                                start_days=[row[self.pmc_start_col] - vaccSP_offset],
                                coverage=high_low_coverages[0],
                                repetitions=1,
                                tsteps_btwn_repetitions=-1,
                                receiving_vaccine_event_name='Received_PMC',
                                ind_property_restrictions=[{'AccessToInterventions': 'higher'}],
                                target_group={'agemin': row[self.pmc_agemin], 'agemax': row[self.pmc_agemax]})

                    # lower access
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=decay_params,
                                start_days=[row[self.pmc_start_col] - vaccSP_offset],
                                coverage=high_low_coverages[1],
                                repetitions=1,
                                tsteps_btwn_repetitions=-1,
                                receiving_vaccine_event_name='Received_PMC',
                                ind_property_restrictions=[{'AccessToInterventions': 'lower'}],
                                target_group={'agemin': row[self.pmc_agemin], 'agemax': row[self.pmc_agemax]})

                else:
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=decay_params,
                                start_days=[row[self.pmc_start_col] - vaccSP_offset],
                                coverage=row[self.pmc_coverage_col],
                                repetitions=1,
                                tsteps_btwn_repetitions=-1,
                                receiving_vaccine_event_name='Received_PMC',
                                target_group={'agemin': row[self.pmc_agemin], 'agemax': row[self.pmc_agemax]})

        else:
            """When running with IIV birth triggered event needs to happen separate from drug campaigns"""
            IIV_groups = ["Group%d" % x for x in range(num_IIV_groups)]
            for r, row in pmc_df.iterrows():

                for index, val in enumerate(IIV_groups):
                    decay_params_IIV = set_IIV_vacc(eff_lower=0.75, eff_upper=0.9, eff_SD=0.025,
                                                    params=decay_params)
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=decay_params_IIV,
                                start_days=[row[self.pmc_start_col] - vaccSP_offset],
                                coverage=row[self.pmc_coverage_col],
                                repetitions=1,
                                tsteps_btwn_repetitions=-1,
                                receiving_vaccine_event_name='Received_PMC',
                                ind_property_restrictions=[{'DrugResponseGroup': val}],
                                target_group={'agemin': row[self.pmc_agemin], 'agemax': row[self.pmc_agemax]})

        return len(pmc_df)


def add_all_interventions(cb, int_suite, my_ds, high_access_ip_frac=0,
                          rtss_target_group='random',
                          cm_target_group='random',
                          pmc_target_group='random',
                          rtss_booster1_min_age=730,
                          hs_df=pd.DataFrame(),
                          rtss_df=pd.DataFrame(),
                          pmc_df=pd.DataFrame(),
                          cohort_month_shift=0):
    event_list = ['Received_NMF_Treatment']

    if not hs_df.empty:
        """per default CM is included in intervention access correlation, otherwise set high_access_ip_frac_cm to 0"""
        has_cm = int_suite.add_ds_hs(cb, hs_df, my_ds, high_access_ip_frac, cm_target_group)
        if has_cm:
            event_list.append('Received_Treatment')
            event_list.append('Received_Severe_Treatment')

    if not rtss_df.empty:
        has_rtss = int_suite.add_ds_rtss(cb, rtss_df, my_ds, high_access_ip_frac, rtss_target_group,
                                         rtss_booster1_min_age, cohort_month_shift)
        if has_rtss > 0:
            event_list.append('Received_Vaccine')

    if not pmc_df.empty:
        has_pmc = int_suite.add_ds_pmc(cb, pmc_df, my_ds, high_access_ip_frac, pmc_target_group)
        if has_pmc > 0:
            event_list.append('Received_PMC')

    event_list = list(np.unique(event_list))
    add_event_counter_report(cb, event_trigger_list=event_list)
    return {}
