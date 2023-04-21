import argparse
import os
import sys
import pandas as pd
import numpy as np
import datetime
import numpy as np
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser

sys.path.append('../../')
from load_paths import load_box_paths
from simulations.analyzer.analyzer_collection import parse_args, make_simout_dir, dailyTreatedCasesAnalyzer, \
    yearlyTreatedCasesAnalyzer

if os.name == "posix":
    SetupParser.default_block = 'NUCLUSTER'
else:
    SetupParser.default_block = 'HPC'
SetupParser.init()


class weeklyTreatedCasesAnalyzer(BaseAnalyzer):

    def __init__(self, expt_name, channels=None, sweep_variables=None, working_dir=".", start_year=2020, end_year=None, filter_exists=False):
        super(weeklyTreatedCasesAnalyzer, self).__init__(working_dir=working_dir,
                                                         filenames=["output/ReportEventCounter.json",
                                                                    "output/ReportMalariaFiltered.json"]
                                                         )
        self.sweep_variables = sweep_variables or ["admin_name", "Run_Number"]
        if channels is None:
            self.channels = ['Received_Treatment', 'Received_Severe_Treatment', 'Received_NMF_Treatment']
        else:
            self.channels = channels
        self.inset_channels = ['Statistical Population', 'New Clinical Cases', 'New Severe Cases', 'PfHRP2 Prevalence']
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year
        self.filter_exists = filter_exists

    def filter(self, simulation):
        if self.filter_exists:
            file = os.path.join(simulation.get_path(), self.filenames[0])
            return os.path.exists(file)
        else:
            return True

    def select_simulation_data(self, data, simulation):

        simdata = pd.DataFrame({x: data[self.filenames[0]]['Channels'][x]['Data'][:1825] for x in self.channels})
        simdata['Time'] = simdata.index

        d = pd.DataFrame({x: data[self.filenames[1]]['Channels'][x]['Data'][:1825] for x in self.inset_channels})
        d['Time'] = d.index

        if len(self.channels) > 0:
            simdata = pd.merge(left=simdata, right=d, on='Time')
        else:
            simdata = d
        simdata['Day'] = simdata['Time'] % 365
        simdata['year'] = simdata['Time'].apply(lambda x: int(x / 365) + self.start_year)
        simdata['date'] = simdata.apply(
            lambda x: datetime.date(int(x['year']), 1, 1) + datetime.timedelta(x['Day']), axis=1)
        # simdata = simdata.loc[simdata['year'] <= self.end_year]

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        if not os.path.exists(os.path.join(self.working_dir, self.expt_name)):
            os.mkdir(os.path.join(self.working_dir, self.expt_name))

        adf = pd.concat(selected).reset_index(drop=True)
        adf['date'] = pd.to_datetime(adf['date'])
        adf['date'] = adf['date'] - pd.to_timedelta(adf['date'].dt.dayofweek, unit='d') + 7
        #adf['year'] = adf['date'].dt.year
        #adf['date'] = adf['date'].dt.week

        sum_channels = self.channels + ['New Clinical Cases', 'New Severe Cases']
        mean_channels = ['Statistical Population', 'PfHRP2 Prevalence']

        df = adf.groupby(self.sweep_variables + ['date'])[sum_channels].agg(np.sum).reset_index()
        pdf = adf.groupby(self.sweep_variables + ['date'])[mean_channels].agg(np.mean).reset_index()

        adf = pd.merge(left=pdf, right=df, on=(self.sweep_variables + ['date']))
        adf.to_csv(os.path.join(self.working_dir, self.expt_name, 'All_Age_weekly_Cases.csv'), index=False)


class monthlyTreatedCasesAnalyzer(BaseAnalyzer):

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, channels=None, sweep_variables=None, working_dir=".", start_year=2020, filter_exists=False):
        super(monthlyTreatedCasesAnalyzer, self).__init__(working_dir=working_dir,
                                                          filenames=["output/ReportEventCounter.json",
                                                                     "output/ReportMalariaFiltered.json"]
                                                          )
        self.sweep_variables = sweep_variables or ["admin_name", "Run_Number"]
        if channels is None:
            self.channels = ['Received_Treatment', 'Received_Severe_Treatment', 'Received_NMF_Treatment']
        else:
            self.channels = channels
        self.inset_channels = ['Statistical Population', 'New Clinical Cases', 'New Severe Cases', 'PfHRP2 Prevalence']
        self.expt_name = expt_name
        self.start_year = start_year
        self.filter_exists = filter_exists

    def filter(self, simulation):
        if self.filter_exists:
            file = os.path.join(simulation.get_path(), self.filenames[0])
            return os.path.exists(file)
        else:
            return True

    def select_simulation_data(self, data, simulation):

        simdata = pd.DataFrame({x: data[self.filenames[0]]['Channels'][x]['Data'] for x in self.channels})
        simdata['Time'] = simdata.index

        d = pd.DataFrame({x: data[self.filenames[1]]['Channels'][x]['Data'] for x in self.inset_channels})
        d['Time'] = d.index

        if len(self.channels) > 0:
            simdata = pd.merge(left=simdata, right=d, on='Time')
        else:
            simdata = d
        simdata['Day'] = simdata['Time'] % 365
        simdata['year'] = simdata['Time'].apply(lambda x: int(x / 365) + self.start_year)
        simdata['date'] = simdata.apply(
            lambda x: datetime.date(int(x['year']), 1, 1) + datetime.timedelta(x['Day'] - 1), axis=1)
        # simdata['date'] =  simdata['date'].apply(lambda x: datetime.round("M"))
        simdata['month'] = simdata['Day'].apply(lambda x: self.monthparser((x + 1) % 365))

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        if not os.path.exists(os.path.join(self.working_dir, self.expt_name)):
            os.mkdir(os.path.join(self.working_dir, self.expt_name))

        adf = pd.concat(selected).reset_index(drop=True)
        adf['date'] = adf.apply(lambda x: datetime.date(x['year'], x['month'], 1), axis=1)

        sum_channels = self.channels + ['New Clinical Cases', 'New Severe Cases']
        mean_channels = ['Statistical Population', 'PfHRP2 Prevalence']

        df = adf.groupby(self.sweep_variables + ['date'])[sum_channels].agg(np.sum).reset_index()
        pdf = adf.groupby(self.sweep_variables + ['date'])[mean_channels].agg(np.mean).reset_index()

        adf = pd.merge(left=pdf, right=df, on=(self.sweep_variables + ['date']))
        adf.to_csv(os.path.join(self.working_dir, self.expt_name, 'All_Age_monthly_Cases.csv'), index=False)


if __name__ == "__main__":

    home_path, data_path, git_dir, project_path = load_box_paths()

    working_dir = os.path.join(project_path, 'simulation_output', '_pmc_rtss_generic')
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)

    """Simulation arguments"""
    args = parse_args()
    exp_name = args.exp_name
    exp_id = args.exp_id
    start_year = args.start_year
    end_year = args.end_year

    """Define experiment sweeps and PMC treatment_channels"""
    if SetupParser.default_block == 'NUCLUSTER':
        from simulations.nucluster_helpers import get_exp_tags, copy_metadata

        filter_exists = True
        tags_dic = get_exp_tags(exp_name, exp_id)
        exp_sweeps = list(tags_dic.keys())
        copy_metadata(exp_name, exp_id, working_dir)

    else:
        filter_exists = False
        exp_sweeps = ['Scenario_id', 'Run_Number', 'Annual EIR', 'seasonality',
                      'cm_coverage', 'rtss_coverage', 'pmc_mode', 'pmc_coverage', 'Cohort_birth_month'
                      ]  # 'cm_target_group', 'pmc_target_group', 'rtss_target_group', 'minBoostAge'  #'Sample_Number',"kmax_SP","C50_SP",

    treatment_channels = ['Received_NMF_Treatment', 'Received_Severe_Treatment', 'Received_Treatment']
    sweep_variables = list(set(["Run_Number"] + exp_sweeps))

    daily_analyzer = [
        dailyTreatedCasesAnalyzer(expt_name=exp_name,
                                  sweep_variables=sweep_variables,
                                  channels=treatment_channels,
                                  working_dir=working_dir,
                                  start_year=start_year,
                                  filter_exists=filter_exists)
    ]
    weekly_analyzer = [
        weeklyTreatedCasesAnalyzer(expt_name=exp_name,
                                   sweep_variables=sweep_variables,
                                   channels=treatment_channels,
                                   working_dir=working_dir,
                                   start_year=start_year,
                                   filter_exists=filter_exists)
    ]
    monthly_analyzer = [
        monthlyTreatedCasesAnalyzer(expt_name=exp_name,
                                    sweep_variables=sweep_variables,
                                    channels=treatment_channels,
                                    working_dir=working_dir,
                                    start_year=start_year,
                                    filter_exists=filter_exists)
    ]
    yearly_analyzer = [
        yearlyTreatedCasesAnalyzer(expt_name=exp_name,
                                   sweep_variables=sweep_variables,
                                   channels=treatment_channels,
                                   working_dir=working_dir,
                                   start_year=start_year,
                                   filter_exists=filter_exists)
    ]
    am = AnalyzeManager(exp_id, analyzers=weekly_analyzer + monthly_analyzer)
    am.analyze()
