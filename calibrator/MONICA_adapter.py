import json
import sys
import monica_io
import zmq
import csv
import os
from datetime import date
import collections
import threading
from threading import Thread
from collections import defaultdict
#from util_calibrator import *
import random
import copy
import monica_germany_utils

from run_producer_consumer import prod_cons_calib

class monica_adapter(object):
    def __init__(self, custom_crop, server, preprocessed_data=None):

        self.sim_id = -1
        self.custom_crop = copy.deepcopy(custom_crop)
        self.species_params = self.custom_crop["cropParams"]["species"]
        self.cultivar_params = self.custom_crop["cropParams"]["cultivar"]
        self.server = server
        self.preprocessed_data = preprocessed_data
        
        #read observations
        obs_pheno = monica_germany_utils.read_stations_obspheno("Z:/projects/monica-germany/DWD_stations.csv",
                                                                "Z:/projects/monica-germany/DWD_1995_2012_obs_phases_WW.csv",
                                                                custom_crop["CROP_ID"])
        
        #populate observation list
        self.observations = []#for spotpy

        #order stations and dates within each station (use the same ordering when collecting results)
        stations = sorted(obs_pheno.keys())
        for stat in stations:
            sorted_by_date = sorted(obs_pheno[stat], key=lambda tup: tup[0])
            for pair in sorted_by_date:
                self.observations.append(pair[1])

    def run(self,args):
        return self._run(*args)

    def _run(self, vector, user_params):

        self.sim_id += 1
        evallist = []
        
        def seek_set_param(par, p_value, model_params):
            p_name = par["name"]
            array = par["array"]
            add_index = False
            if isinstance(model_params[p_name], int) or isinstance(model_params[p_name], float):
                add_index = False
            elif len(model_params[p_name]) > 1 and isinstance(model_params[p_name][1], basestring):
                add_index = True #the param contains text (e.g., units)
            if array.upper() == "FALSE":
                if add_index:
                    model_params[p_name][0] = p_value
                else:
                    model_params[p_name] = p_value
            else: #param is in an array (possibly nested)
                array = array.split("_") #nested array
                if add_index:
                    array = [0] + array
                if len(array) == 1:
                    model_params[p_name][int(array[0])] = p_value
                elif len(array) == 2:
                    model_params[p_name][int(array[0])][int(array[1])] = p_value
                elif len(array) == 3:
                    model_params[p_name][int(array[0])][int(array[1])][int(array[2])] = p_value
                else:
                    print "param array too nested, contact developers"
            

        #set params according to spotpy sampling.
        for i in range(len(user_params)):
            if user_params[i]["name"] in self.species_params:
                seek_set_param(user_params[i],
                user_params[i]["derive_function"](vector, self.species_params) if "derive_function" in user_params[i] else vector[i],
                self.species_params)
            if user_params[i]["name"] in self.cultivar_params:
                seek_set_param(user_params[i],
                user_params[i]["derive_function"](vector, self.cultivar_params) if "derive_function" in user_params[i] else vector[i],
                self.cultivar_params)
        
        
        sim_pheno = prod_cons_calib(self.custom_crop, self.server, self.preprocessed_data)
        
        #order sim according to the same criteria applied for self.observations
        stations = sorted(sim_pheno.keys())
        for stat in stations:
            sorted_by_date = sorted(sim_pheno[stat], key=lambda tup: tup[0])
            for pair in sorted_by_date:
                evallist.append(pair[1])

        return evallist
        
        