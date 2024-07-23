"""
Module Name: individual_sigprofiler.py

Description:
This module enables SigProfilerAssignment to be run individually on each of the files within the data_dir
sub-directories that are specified within the list cancer_types.
Once SigProfilerAssignment.Analyzer has run on a file the results are examined and the data in
'Composition After Add-Remove' is added to an all_results.json file, whilst the other results are thrown away.
The all_results.json and all_errors.json files are kept within the analysis_dir, as well as the temporary folders.
The current version works with GrCh38.



Author: Sarah Wooller
"""

import shutil
import os
import pandas as pd
import json
from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerAssignment import Analyzer as Analyze

def clear_directory(directory):
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)

def log_error(e,cancer_type, file, analysis_dir):
    error_path = os.path.join(analysis_dir, 'errors.json')
    print(e)
    if os.path.isfile(error_path):
        with open(error_path) as fp:
            errors = json.load(fp)
    else:
        errors = dict()

    cancer_errors = errors.get(cancer_type,[])
    cancer_errors.append([file.split('.')[0],str(e)])
    errors[cancer_type] = cancer_errors

    with open(error_path, 'w') as fp:
        json.dump(errors, fp)


def fit_file(file, cancer_type, data_dir, analysis_dir):

    all_results = os.path.join(analysis_dir, "all_results")
    temp_input = os.path.join(analysis_dir, "all_input")
    input_file = os.path.join(data_dir, cancer_type, file)
    os.makedirs(all_results, exist_ok=True)
    os.makedirs(temp_input, exist_ok=True)
    clear_directory(temp_input)
    clear_directory(all_results)
    temp_input_file = os.path.join(temp_input, file)
    shutil.copy(input_file, temp_input_file)
    Analyze.cosmic_fit(temp_input,
                       all_results,
                       input_type="vcf",
                       context_type="96",
                       collapse_to_SBS96=True,
                       cosmic_version=3.4,
                       exome=False,
                       genome_build="GRCh38")


def process_results(cancer_type, file, data_dir, analysis_dir):
    sep = '####################################### Composition After Add-Remove #######################################\n'
    all_results = os.path.join(data_dir, "all_results")
    result = os.path.join(all_results,
                          'Assignment_Solution/Solution_Stats/Assignment_Solution_Signature_Assignment_log.txt')
    with open(result) as f:
        text = f.read().split(sep)[1].split('\n')[:2]

    columns = ['SBS1', 'SBS2', 'SBS3', 'SBS4', 'SBS5', 'SBS6', 'SBS7a', 'SBS7b',
               'SBS7c', 'SBS7d', 'SBS8', 'SBS9', 'SBS10a', 'SBS10b', 'SBS10c',
               'SBS10d', 'SBS11', 'SBS12', 'SBS13', 'SBS14', 'SBS15', 'SBS16',
               'SBS17a', 'SBS17b', 'SBS18', 'SBS19', 'SBS20', 'SBS21', 'SBS22a',
               'SBS22b', 'SBS23', 'SBS24', 'SBS25', 'SBS26', 'SBS27', 'SBS28', 'SBS29',
               'SBS30', 'SBS31', 'SBS32', 'SBS33', 'SBS34', 'SBS35', 'SBS36', 'SBS37',
               'SBS38', 'SBS39', 'SBS40a', 'SBS40b', 'SBS40c', 'SBS41', 'SBS42',
               'SBS43', 'SBS44', 'SBS45', 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50',
               'SBS51', 'SBS52', 'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57', 'SBS58',
               'SBS59', 'SBS60', 'SBS84', 'SBS85', 'SBS86', 'SBS87', 'SBS88', 'SBS89',
               'SBS90', 'SBS91', 'SBS92', 'SBS93', 'SBS94', 'SBS95', 'SBS96', 'SBS97',
               'SBS98', 'SBS99']

    remove_space = lambda thelist: [i for i in thelist if i != ""]
    index = remove_space(text[0].split(' '))
    values = remove_space(text[1].split(' ')[1:])
    name = f"{cancer_type}_{file.split('.')[0]}"
    all_results_path = os.path.join(analysis_dir, "all_results.json")
    if not os.path.isfile(all_results_path):
        all_results = dict()
    else:
        with open(all_results_path) as fp:
            all_results = json.load(fp)

    all_results[name] = {'index': index, 'values': values}
    with open(all_results_path, 'w') as fp:
        json.dump(all_results, fp)

def main():
    genInstall.install('GRCh38')
    cancer_types = [
     'lung_WT',
     'lung_G12D',
     'lung_G12R',
     'lung_G12V',
     'lung_G12V_results',
     'large_intestine_WT',
     'large_intestine_G12D',
     'pancreas_G12V',
     'large_intestine_G12C',
     'large_intestine_G12R',
     'pancreas_WT',
     'pancreas_G12C',
     'large_intestine_G12V',
     'lung_G12C']
    analysis_dir = '/home/skw24/krasm'
    data_dir = '/home/skw24/cathy/kras'
    all_results_path = os.path.join(analysis_dir, "all_results.json")

    for cancer_type in cancer_types:
        input_files = os.listdir(os.path.join(data_dir, cancer_type))
        for file in input_files:
            with open(all_results_path) as f:
                all_results = json.load(f)
            name = f"{cancer_type}_{file.split('.')[0]}"
            if name not in all_results:
                print(cancer_type, file)
                try:
                    fit_file(file, cancer_type, data_dir)
                    process_results(cancer_type, file, data_dir, analysis_dir)
                    print(f'Success for {cancer_type} {file}')
                except KeyboardInterrupt:
                    break
                except Exception as e:
                    log_error(e,cancer_type, file, analysis_dir)
            else:
                print(f'{cancer_type} {file} already analysed')



