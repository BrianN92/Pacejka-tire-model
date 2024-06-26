'''
If the left and right tires have different nominal wheel loads, 
you would typically use a weighted average to combine the Pacejka tire parameters 
for a more accurate representation in a bicycle model. 

The weighting factor often corresponds to the wheel load or the force experienced by each tire.
'''

import numpy as np
import matplotlib.pyplot as plt
import yaml
import os
import time

def combineTires(tire1, tire2):
    FNOMIN1 = tire1['FNOMIN']
    FNOMIN2 = tire2['FNOMIN']
    tire_combine = tire1
    for key in tire1.keys():
        if isinstance(tire1[key],float):
            tire_combine[key] = round((tire1[key]*FNOMIN1 + tire2[key]*FNOMIN2) / (FNOMIN1+FNOMIN2),2)
    return tire_combine

def load_params(tire_select, data_path):
    data_ls = [data for data in os.listdir(data_path)]
    idx_ls  = [idx for idx, ltr in enumerate(data_ls) if tire_select in ltr]
    name    = data_ls[idx_ls[0]]
    with open(data_path + name, 'r') as stream:
        try:
            parsed_yaml = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return parsed_yaml
def saveToYaml(params, path, name):
    with open(path + name + '.yaml', 'w') as outfile:
        yaml.dump(params, outfile)
    outfile.close()
    

