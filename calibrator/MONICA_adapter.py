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
from util_calibrator import *
import random
import copy

from run_producer_consumer import prod_cons_calib

class monica_adapter(object):
    def __init__(self, custom_crop, server):

        self.all_years = range(1999, 2013) #both obs and sims available for these years!
        self.sim_id = -1
        self.custom_crop = copy.deepcopy(custom_crop)
        self.species_params = self.custom_crop["cropParams"]["species"]
        self.cultivar_params = self.custom_crop["cropParams"]["cultivar"]
        self.server = server
        
        #read observations
        sim_lk = set()
        with open("calculate-indices/landkreise_bkrs.csv") as _:
            reader = csv.reader(_)
            next(reader, None)
            for row in reader:
                sim_lk.add(int(row[0]))

        #data[bkr/lk][year][sim/obs]
        self.yield_data = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
        self.obs_lks = set()
        with open("calculate-indices/official_yields_DE.csv") as _:
            reader = csv.reader(_)
            for i in range(7):
                next(reader, None)

            #3= Winterweizen (WW)
            for row in reader:
                if len(row) > 2 and representsInt(row[1]):
                    lk = int(row[1])
                    if lk in sim_lk:
                        year = int(row[0])
                        if representsfloat(row[3]):
                            if float(row[3]) == 0.0:
                                #consider 0 as nodata
                                continue
                            self.obs_lks.add(lk)
                            obs_yield = float(row[3]) * 100 #kg ha-1
                            self.yield_data[lk][year]["obs"] = obs_yield
        
        #populate observation list
        self.observations = []#for spotpy

        for lk in self.obs_lks:
            yrs = find_obs_years(id=lk, data=self.yield_data, yr_range=self.all_years)
            obs = retrieve_data(ids=[lk], years=yrs, data=self.yield_data, obs_sim="obs")
            self.observations += obs

        #print len(self.observations)


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
        
        
        sim_yield_DE = prod_cons_calib(self.custom_crop, self.server)
        
        #####temporary way to populate sim_yield_DE (for testing)
        #sim_yield_DE = defaultdict(lambda: defaultdict(float))
        #out_dir = "calculate-indices/test_out/"
        #random_multiplier = random.uniform(0, 2)
        #for filename in os.listdir(out_dir):
        #    if "yield" not in filename:
        #        continue
        #    year = int(filename.split("_")[2])
        #    with open(out_dir + filename) as _:
        #        reader = csv.reader(_)
        #        next(reader, None)
        #        for row in reader:
        #            lk = int(row[0])
        #            sim_yield_DE[year][lk] = float(row[1]) * random_multiplier

        for yr in sim_yield_DE:
            for lk in sim_yield_DE[yr]:
                self.yield_data[lk][yr]["sim"] = sim_yield_DE[yr][lk]
        
        #populate evallist
        for lk in self.obs_lks:
            yrs = find_obs_years(id=lk, data=self.yield_data, yr_range=self.all_years)
            sim = retrieve_data(ids=[lk], years=yrs, data=self.yield_data, obs_sim="sim")
            #sim = retrieve_data(ids=[lk], years=yrs, data=self.yield_data, obs_sim="obs") #use the observations as sim
            evallist += sim
        #print len(evallist)

        return evallist
        
        