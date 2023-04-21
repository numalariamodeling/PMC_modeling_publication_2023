import os
import pandas as pd
import numpy as np
import json
import requests
import time
# from dtk.tools.demographics.DemographicsGeneratorConcern import WorldBankBirthRateConcern, \
#     EquilibriumAgeDistributionConcern, DefaultIndividualAttributesConcern
# from dtk.tools.demographics.DemographicsGenerator import DemographicsGenerator
from input_file_generation.add_properties_to_demographics import generate_demographics_properties
# from dtk.tools.demographics.DemographicsGenerator import DemographicsGenerator  # use existing demographics json
import sys

sys.path.append('../')
from load_paths import load_box_paths

datapath, projectpath = load_box_paths()
inputs_path = os.path.join(projectpath, 'simulation_inputs/generic_cohort_ds/Demographics')


def check_demo(refdemo_fname, overwrite=True):
    with open(os.path.join(inputs_path, refdemo_fname)) as fin:
        demo = json.loads(fin.read())
    # all_nodeids = [x['NodeID'] for x in demo['Nodes']]

    if 'Defaults' not in demo.keys():
        demo['Defaults'] = {'NodeAttributes': {}, 'IndividualAttributes': {}, 'IndividualProperties': {}}
        demo['Defaults']['NodeAttributes'] = demo['Nodes']['NodeID' == 1]['NodeAttributes']
        demo['Defaults']['IndividualAttributes'] = demo['Nodes']['NodeID' == 1]['IndividualAttributes']
        demo['Defaults']['IndividualProperties'] = demo['Nodes']['NodeID' == 1]['IndividualProperties']

        demo_fname = os.path.join(inputs_path, refdemo_fname)
        if overwrite:
            demo_fname2 = os.path.join(inputs_path, refdemo_fname)
        else:
            demo_fname2 = os.path.join(inputs_path, refdemo_fname.replace('.json', '_wDefault.json'))

        def default(o):
            if isinstance(o, np.int64): return int(o)
            raise TypeError

        with open(demo_fname2, 'w') as fout:
            json.dump(demo, fout, sort_keys=True, indent=4, separators=(',', ': '), default=default)
    else:
        pass


def add_Vaccine_IPs(ffile, overwrite=False):
    """Add VaccineStatus IP"""
    IPs = [{'Property': 'VaccineStatus',
            'Values': ['None', 'GotVaccine', 'GotBooster1', 'GotBooster2'],
            'Initial_Distribution': [1.0, 0.0, 0.0, 0.0],
            'Transitions': []}
           ]

    adf = pd.DataFrame({'Property': 'VaccineStatus',
                        'Property_Value': ['None', 'GotVaccine', 'GotBooster1', 'GotBooster2'],
                        'Initial_Distribution': [1.0, 0.0, 0.0, 0.0]}
                       )
    adf['Property_Type'] = 'IP'
    adf['node'] = 1

    demo_fname = os.path.join(inputs_path, ffile)
    if overwrite:
        IP_demo_fname = os.path.join(inputs_path, ffile)
    else:
        IP_demo_fname = os.path.join(inputs_path, ffile.replace('.json', '_wVax.json'))

    generate_demographics_properties(refdemo_fname=demo_fname,
                                     output_filename=IP_demo_fname,
                                     as_overlay=False,
                                     IPs=IPs,
                                     df=adf)


def add_DrugResponseGroup_IPs(ffile, overwrite=False):
    """Add DrugResponseGroups for IIV"""
    IPs = [{'Property': 'DrugResponseGroup',
            'Values': [f'Group{i}' for i in range(100)],
            'Initial_Distribution': [0.01] * 100,
            'Transitions': []}
           ]

    adf = pd.DataFrame({'Property': ['DrugResponseGroup'] * 100,
                        'Property_Value': [f'Group{i}' for i in range(100)],
                        'Initial_Distribution': [0.01] * 100}
                       )
    adf['Property_Type'] = 'IP'
    adf['node'] = 1

    demo_fname = os.path.join(inputs_path, ffile)
    if overwrite:
        IP_demo_fname = os.path.join(inputs_path, ffile)
    else:
        IP_demo_fname = os.path.join(inputs_path, ffile.replace('.json', '_IIV.json'))

    generate_demographics_properties(refdemo_fname=demo_fname,
                                     output_filename=IP_demo_fname,
                                     as_overlay=False,
                                     IPs=IPs,
                                     df=adf)


if __name__ == '__main__':

    VaccineStatus_IP = False
    DrugResponseGroup_IP = True

    """Customize as needed, inputs_path, file pattern and which IPs to add"""
    inputs_path = os.path.join(inputs_path)

    ffiles = [x for x in os.listdir(inputs_path) if '.json' in x and not 'IIV' in x]
    for ffile in ffiles:
        print(ffile)
        check_demo(refdemo_fname=ffile, overwrite=True)  # not critical to overwrite here
        if VaccineStatus_IP:
            add_Vaccine_IPs(ffile=ffile,
                            overwrite=True)  # Note only run once as overwriting to keep same filename

        if DrugResponseGroup_IP:
            add_DrugResponseGroup_IPs(ffile=ffile, overwrite=False)
