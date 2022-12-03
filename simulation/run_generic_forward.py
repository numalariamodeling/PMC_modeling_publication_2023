import os
import sys
import logging
import argparse
import pandas as pd
import time
import random
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.SetupParser import SetupParser
from hbhi.set_up_general import set_spaq_params
from malaria.interventions.malaria_drugs import set_drug_param, get_drug_param



sys.path.append('../')
from pmc_rtss_generic.setup_helper import setup_setting, setup_simulation, add_generic_interventions
from simulations.intervention_helpers import set_SDX_PYR_with_IIV, set_SP_with_IIV, set_SDX_PYR_params
from simulations.nucluster_helpers import create_submisson_scripts_2
from load_paths import load_box_paths

logger = logging.getLogger(__name__)
if os.name == "posix":
    SetupParser.default_block = 'NUCLUSTER'
    rtss_dir = 'rtss-scenario_IO'
    x_pop_scale = 1
else:
    SetupParser.default_block = 'HPC'
    rtss_dir = 'rtss_scenarios'
    x_pop_scale = 1


def parse_args():
    description = "Simulation specifications"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-si",
        "--scen_index",
        type=int,
        help="Index of simulation experiment (see order in scenario_csv_list)",
        default=None

    )
    parser.add_argument(
        "-sc",
        "--scen_csv",
        type=str,
        help="Name of simulation experiment",
        default=None

    )
    parser.add_argument(
        "--test",
        action='store_true',
        help="If true, runs as test with smaller population, fewer seeds and single birth cohort",
    )
    return parser.parse_args()


args = parse_args()
testrun = args.test
scen_csv = args.scen_csv

_, _, _, projectpath = load_box_paths()
input_path = os.path.join(projectpath, 'simulation_inputs/generic_cohort_ds')

num_seeds = 5
years = 6
use_12_cohorts_flag = True  # simulate one cohort born in each month (seasonality dates updated accordingly)
pmc_mode = 'withIIV'  # 'noIIV' # 'withIIV'
num_IIV_groups = 100
pmc_drug = 'vaccSP'  # 'SP' or 'SDX_PYR', or 'vaccSP'
PDsweep = False
num_params_to_run = 15  # if PDsweep True

scen_csv = 'generic_PMC3mode_shifted.csv'
if scen_csv is None:
    si = args.scen_index
    scenario_csv_list = ['generic_single_RTSS.csv',
                         'generic_single_PMC.csv',
                         'generic_PMC_RTSS_EIR.csv',
                         'generic_PMCmode_RTSS.csv',
                         'generic_PMCmode_EIR.csv',
                         'generic_PMCmode_RTSS_cov.csv',
                         'generic_PMC_CM_accesscor.csv',
                         'generic_PMC_CM_accesscor_fixedFracH.csv',
                         'generic_PMC3mode_shifted']
    scen_csv = scenario_csv_list[si]

if scen_csv =='generic_single_RTSS.csv':
    years = 6
  
if scen_csv =='generic_single_PMC.csv':
   years = 3
      
exp_name = scen_csv.replace('.csv', '')
exp_name = exp_name + f'_{pmc_drug}'
if pmc_mode == 'withIIV':
    exp_name = exp_name + '_IIV'
if PDsweep:
    exp_name = exp_name + '_PDsweep'

# setup basic config for simulations
cb = setup_simulation(years)
scen_df = pd.read_csv(os.path.join(input_path, scen_csv), encoding='latin')

### CUSTOM DRUG SAMPLES FOR IPTI
scen_df['drug_code'] = pmc_drug

if pmc_mode == 'withIIV':
  scen_df['num_IIV_groups'] = num_IIV_groups

if PDsweep:

    if pmc_drug == 'SDX_PYR':
        set_SDX_PYR_with_IIV(cb, num_IIV_groups)
        top_param_df = pd.read_csv(os.path.join(projectpath, 'simulation_inputs', 'SP_PDsweep', 'distscoreo.csv'))
    elif pmc_drug == 'SP':
        set_SP_with_IIV(cb, num_IIV_groups)
        top_param_df = pd.read_csv(
            os.path.join(projectpath, 'simulation_inputs', 'SP_PDsweep', 'distscoreo.csv'))  # FIXME
    else:
        raise ValueError('SPdrug needs to be "SP" or "SDX_PYR"')

scen_df['DS_Name'] = scen_df['setting_id'].str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode(
    'utf-8')
scen_df = scen_df.set_index('DS_Name')
eir_monthly_multipliers = pd.read_csv( os.path.join(input_path, 'Seasonality','seasonality_eir_NGApmc_multipliers.csv'))

if use_12_cohorts_flag:
    cohort_month_shift_values = range(12)
else:
    cohort_month_shift_values = [0]

# REPORTS
summaryreport = False
if summaryreport:
    from malaria.reports.MalariaReport import add_summary_report

    add_summary_report(cb, start=0, interval=365, duration_days=365, age_bins=[0.167, 1, 125], description='U1')
    add_summary_report(cb, start=0, interval=365 * 2, duration_days=365 * 2, age_bins=[0.167, 2, 125], description='U2')
    add_summary_report(cb, start=0, interval=365 * 5, duration_days=365 * 5, age_bins=[0.167, 5, 125], description='U5')
    add_summary_report(cb, start=0, interval=365 * 10, duration_days=365 * 10, age_bins=[0.167, 10, 125],
                       description='U10')
    add_summary_report(cb, start=0, interval=30, duration_days=365 * 5, age_bins=[0.167, 5, 125],
                       description='U5_monthly')


def set_sp_params(cb, c50_scl=1, kill_scl=1, drug_name='SP'):
    if drug_name == 'SP':
        set_drug_param(cb, 'Sulfadoxine', "Drug_PKPD_C50", 0.2 * c50_scl)
        set_drug_param(cb, 'Sulfadoxine', "Max_Drug_IRBC_Kill", 0.506 * kill_scl)  # 0.675
        set_drug_param(cb, 'Pyrimethamine', "Drug_PKPD_C50", 2)  # *6
        set_drug_param(cb, 'Pyrimethamine', "Max_Drug_IRBC_Kill", 0.6)  # *0.6417
    if drug_name == 'SDX_PYR':
        set_drug_param(cb, 'SulfadoxinePyrimethamine', "Drug_PKPD_C50", 0.2 * c50_scl)  # 54.1
        set_drug_param(cb, 'SulfadoxinePyrimethamine', "Max_Drug_IRBC_Kill", 0.506 * kill_scl)  # 3.120

    return {'drug_name': drug_name,
            'c50_scl': c50_scl,
            'kill_scl': kill_scl}


if testrun:
  print('Running simulation as TEST')
  num_seeds = 3
  years = 5
  use_12_cohorts_flag = False 
  cb.update_params({'x_Base_Population': 0.5})
  exp_name = exp_name + '_testrun'

matAbs = [0.1327]
if scen_csv =='generic_PMC_matAb.csv':
  matAbs = [0.0, 0.1327, 0.2, 0.4, 0.6, 0.8, 1.0]
  cohort_month_shift_values = [0]

# BUILDER
builder = ModBuilder.from_list([[ModFn(setup_setting,
                                       scen_df=scen_df,
                                       id=id,
                                       eir_monthly_multipliers=eir_monthly_multipliers,
                                       EIR_scale='daily',
                                       cohort_month_shift=cohort_month_shift,
                                       Maternal_Antibody_Protection=mAb_Protection),
                                 ModFn(add_generic_interventions,
                                       input_path=input_path,
                                       scen_df=scen_df,
                                       id=id,
                                       cohort_month_shift=cohort_month_shift),
                                 # ModFn(set_spaq_params),  ###<-- currently the dtk-tools-malaria parameters are the most updated
                                 # ModFn(set_sp_params, c50_scl=c50_scl,kill_scl=kill_scl, drug_name=pmc_drug),
                                 # ModFn(set_SDX_PYR_params, c50=row['C50_SDX_PYR'], kmax=row['kmax_SDX_PYR'], num_groups=num_IIV_groups),
                                 # ModFn(set_SP_params, row=row, num_groups=num_IIV_groups, pmc_drug= pmc_drug),
                                 # ModFn(DTKConfigBuilder.set_param, 'Sample_Number', r),
                                 ModFn(DTKConfigBuilder.set_param, 'Setting_id', id),
                                 ModFn(DTKConfigBuilder.set_param, 'matAbs', mAb_Protection),
                                 ModFn(DTKConfigBuilder.set_param, 'Run_Number', x),
                                 ModFn(DTKConfigBuilder.set_param, 'Cohort_birth_month', cohort_month_shift),
                                 ]
                                # for id in ['HX6']  # Arbitrary setting for testing
                                for id in scen_df.index
                                for mAb_Protection in matAbs
                                for cohort_month_shift in cohort_month_shift_values
                                # for c50_scl in [1, 6, 54.1]
                                # for kill_scl in [1, 0.675, 3.120]
                                # for r, row in top_param_df[:num_params_to_run].iterrows()
                                for x in range(num_seeds)
                                ])

if __name__ == "__main__":

    run_sim_args = {
        'exp_name': f'mrm9534_{exp_name}',
        'config_builder': cb,
        'exp_builder': builder
    }

    SetupParser.init()
    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)

    if os.name == "posix":
        create_submisson_scripts_2(exp_name=exp_name, exp_id=exp_manager.experiment.exp_id,
                                   analyzer_script="analyzer_cohort_generic.py")
        time.sleep(20)
        exp_manager.wait_for_finished(verbose=True)
        assert (exp_manager.succeeded())
