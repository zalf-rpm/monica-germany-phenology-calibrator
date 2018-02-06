import os
import spotpy
from collections import defaultdict
import csv
import numpy as np

def representsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def representsfloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def retrieve_data(ids, years, data, obs_sim=""):
    'returns data from selected bkrs/lks (i.e., ids) and years. obs_sim must be either "obs" or "sim"'
    out = []
    for id in ids:
        if id not in data.keys():
            continue
        for yr in years:
            if yr not in data[id].keys():
                continue
            if obs_sim not in data[id][yr].keys():
                continue
            out.append(data[id][yr][obs_sim])
    return out

def find_obs_years(id, data, yr_range):
    'returns a list of years for which observed yields are available in the lk or bkr'
    out = []
    for yr in data[id].keys():
        if yr in yr_range and "obs" in data[id][yr].keys():
            out.append(yr)
    return out
