import os
import sys
import json
import numpy as np
import pandas as pd


def get_simtools(WDIR=os.getcwd()):
    """Assumes NUCLUSTER is first section in simtools.ini """
    file = open(os.path.join(WDIR, 'simtools.ini'), 'r')
    contents = file.read()
    contents_list = contents.split("\n")

    if os.name == "posix":
        input_root = [item for item in contents_list if "input_root" in item][0].split(" = ")[1]
        sim_root = [item for item in contents_list if "sim_root" in item][0].split(" = ")[1]
    else:
        input_root = [item for item in contents_list if "input_root" in item][1].split(" = ")[1]
        sim_root = [item for item in contents_list if "sim_root" in item][1].split(" = ")[1]

    return input_root, sim_root


def get_most_recent_experiment_id_by_name(exp_name):
    input_root, sim_root = get_simtools()

    files = os.listdir(sim_root)
    exp_ids = [x.split('___')[1] for x in files if exp_name in x]
    latest_exp_id_int = max([int(x) for x in exp_ids])
    latest_exp_id = [exp_id for exp_id in exp_ids if int(exp_id) == latest_exp_id_int][0]

    return latest_exp_id


def get_exp_tags(exp_name, exp_id, wdir=os.getcwd(), exp_dir=None):
    """Get experiment tags (sweeps) from metadata.json"""
    if exp_dir is None:
        input_root, sim_root = get_simtools(WDIR=wdir)
        exp_dir = os.path.join(sim_root, f'{exp_name}___{exp_id}')
    file = open(os.path.join(exp_dir, 'metadata.json'), 'r')
    data = json.load(file)
    tags_dic = data['simulations'][0]['tags']

    return tags_dic


def copy_metadata(exp_name, exp_id, simout_dir):
    input_root, sim_root = get_simtools()
    exp_dir = os.path.join(sim_root, f'{exp_name}___{exp_id}')
    try:
        os.popen(f"cp {os.path.join(exp_dir, 'metadata.json')} {os.path.join(simout_dir,exp_name, 'metadata.json')}")
        print('metadata.json copied')
    except:
        print('metadata.json not copied')


def check_failed_simulations(exp_dir, n_outputfiles=1):
    header_post = shell_header(job_name=f'failedsims', t='01:00:00', memG=3)
    pymodule = '\n\nmodule purge all' \
               '\nmodule load python/anaconda3.6'
    pycommand = f'\ncd /projects/b1139/templates/' \
                f'\npython run_submit_failed_simulations.py --exp_dir {exp_dir} --nfiles {n_outputfiles}'
    file = open(os.path.join(exp_dir, 'run_submit_failed_simulations.sh'), 'w')
    file.write(header_post + pymodule + pycommand)
    file.close()


def create_submisson_scripts_2(exp_name, exp_id, WDIR=os.getcwd(), A='b1139', p='b1139', t='04:00:00',
                               analyzer_script=None):
    if analyzer_script == None:
        analyzer_script = "analyzer_cohort.py"
    userp = '/home/mrm9534/'
    job_name = f'analyze_{exp_name}'
    header = f'#!/bin/bash\n#SBATCH -A {A}\n#SBATCH -p {p}\n#SBATCH -t {t}\n#SBATCH -N 1\n' \
             f'#SBATCH --ntasks-per-node=1\n#SBATCH --mem=80G\n#SBATCH --job-name="{job_name}"\n'
    err = '#SBATCH --error=log/slurm_%A_%a.err\n'
    out = '#SBATCH --output=log/slurm_%A_%a.out\n'
    header_post = header + err + out
    pymodule = f'\n\nmodule purge all\nmodule load python/anaconda3.6\nsource activate {userp}environments/dtk-tools-p36\n'
    pycommand1 = f'\ncd {userp}dtk-tools-p36/helper_tools/'
    pycommand2 = f'\npython local_db_fixer.py  -id {exp_id} --status Succeeded'
    pycommand3 = f'\ncd {userp}gitrepos/ipti_pmc/simulations/analyzer/ \npython {analyzer_script} --exp_name {exp_name} --exp_id {exp_id}'
    file = open(os.path.join(WDIR, f'run_analyzer_{exp_id}.sh'), 'w')
    file.write(header_post + pymodule + pycommand1 + pycommand2 + pycommand3)
    file.close()


def create_submisson_scripts(exp_name, exp_id, WDIR=os.getcwd(), analyzer_script=None):
    """Analyzer"""
    userp = '/home/mrm9534/'
    input_root, sim_root = get_simtools()
    exp_dir = os.path.join(sim_root, f'{exp_name}___{exp_id}')
    header_post = shell_header(job_name=f'analyze_{exp_name}', t='04:00:00', memG=50)
    if analyzer_script is None:
        analyzer_script = "agebins_analyzer_summaryReport.py"
    pymodule = '\n\nmodule purge all' \
               '\nmodule load python/anaconda3.6' \
               f'\nsource activate {userp}environments/dtk-tools-p36\n' \
               f'\ncd {userp}dtk-tools-p36/helper_tools/\n' \
               f'\npython local_db_fixer.py  -id {exp_id} --status Succeeded\n'
    pycommand = f'\ncd {WDIR}/analyzer ' \
                f'\npython {analyzer_script} --exp_name {exp_name} --exp_id {exp_id}'
    file = open(os.path.join(exp_dir, 'run_analyzer.sh'), 'w')
    file.write(header_post + pymodule + pycommand)
    file.close()

    """Check failed simulations"""
    check_failed_simulations(exp_dir=exp_dir, n_outputfiles=2)  # at least 2 (varies depending on experiment)

    """Sample plot"""
    header_post = shell_header(job_name=f'sample_plot_{exp_name}', t='00:20:00', memG=3)
    pymodule = '\n\nmodule purge all' \
               '\nmodule load python/anaconda3.6' \
               '\nsource activate /projects/p30781/anaconda3/envs/team-test-py37\n'
    exp_sweep = 'Run_Number'
    pycommand = f'\ncd ../Py_plotter ' \
                f'\npython sample_plot.py --exp_name {exp_name} --exp_id {exp_id} --exp_sweep {exp_sweep}'
    file = open(os.path.join(exp_dir, 'run_sampleplot.sh'), 'w')
    file.write(header_post + pymodule + pycommand)
    file.close()


def shell_header(A='b1139', p='b1139', t='02:00:00', N=1, ntasks_per_node=1, memG=8, job_name='myjob', arrayJob=None):
    """Requires a 'log' subfolder to write in .err and .out files, alternatively log/ needs to be removed"""
    header = f'#!/bin/bash\n' \
             f'#SBATCH -A {A}\n' \
             f'#SBATCH -p {p}\n' \
             f'#SBATCH -t {t}\n' \
             f'#SBATCH -N {N}\n' \
             f'#SBATCH --ntasks-per-node={ntasks_per_node}\n' \
             f'#SBATCH --mem={memG}G\n' \
             f'#SBATCH --job-name="{job_name}"\n'
    if arrayJob is not None:
        array = arrayJob
        err = '#SBATCH --error=log/slurm_%A_%a.err\n'
        out = '#SBATCH --output=log/slurm_%A_%a.out\n'
        header = header + array + err + out
    else:
        err = f'#SBATCH --error=log/{job_name}.%j.err\n'
        out = f'#SBATCH --output=log/{job_name}.%j.out\n'
        header = header + err + out
    return header
