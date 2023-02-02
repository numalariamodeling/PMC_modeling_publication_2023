import os
import sys
import logging
import argparse
import pandas as pd
import time
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.SetupParser import SetupParser

sys.path.append('../')
from pmc_rtss_generic.setup_helper import setup_setting, setup_simulation, add_generic_interventions
from simulations.nucluster_helpers import create_submisson_scripts_2
from load_paths import load_box_paths

logger = logging.getLogger(__name__)
SetupParser.default_block = 'NUCLUSTER'


def parse_args():
    description = "Simulation specifications"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-sc",
        "--scen_csv",
        type=str,
        help="Name of simulation experiment",
        default=None)
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

exp_name = scen_csv.replace('.csv', '') + '_vaccSP'

if pmc_mode == 'withIIV':
    exp_name = exp_name + '_IIV'

# setup basic config for simulations
cb = setup_simulation(years)

scen_df = pd.read_csv(os.path.join(input_path, scen_csv), encoding='latin')

if pmc_mode == 'withIIV':
    scen_df['num_IIV_groups'] = num_IIV_groups

scen_df['DS_Name'] = scen_df['setting_id'].str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode(
    'utf-8')
scen_df = scen_df.set_index('DS_Name')
eir_monthly_multipliers = pd.read_csv(os.path.join(input_path, 'Seasonality', 'seasonality_eir_NGApmc_multipliers.csv'))

if use_12_cohorts_flag:
    cohort_month_shift_values = range(12)
else:
    cohort_month_shift_values = [0]

matAbs = [0.1327]
if scen_csv == 'generic_PMC_matAb.csv':
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
                                 ModFn(DTKConfigBuilder.set_param, 'Setting_id', id),
                                 ModFn(DTKConfigBuilder.set_param, 'matAbs', mAb_Protection),
                                 ModFn(DTKConfigBuilder.set_param, 'Run_Number', x),
                                 ModFn(DTKConfigBuilder.set_param, 'Cohort_birth_month', cohort_month_shift),
                                 ]
                                for id in scen_df.index
                                for mAb_Protection in matAbs
                                for cohort_month_shift in cohort_month_shift_values
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
