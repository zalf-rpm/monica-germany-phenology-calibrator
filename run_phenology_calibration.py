#!/usr/bin/python
# -*- coding: UTF-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Michael Berg-Mohnicke <michael.berg@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# This file has been created at the Institute of
# Landscape Systems Analysis at the ZALF.
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

import time
import os
import math
import json
import csv
import copy
from StringIO import StringIO
from datetime import date, datetime, timedelta
from collections import defaultdict
#import types
import sys
sys.path.append("./calibrator")
#print sys.path
import zmq
#print "pyzmq version: ", zmq.pyzmq_version(), " zmq version: ", zmq.zmq_version()

import sqlite3
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from pyproj import Proj, transform

from multiprocessing import Process

from sampler_MONICA import start_calibration
import run_work_producer
import run_work_consumer

def main():

    config = {
        "server": "cluster2",
        "prod-port": "6666",
        "cons-port": "7777",
        "nd-port": "5555",
        "run-ids": "all", #"[9, 10, 11, 12, 13, 14, 15, 16, 25, 26, 27, 28, 29, 30, 31, 32]"
    }
    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            k,v = arg.split("=")
            if k in config:
                config[k] = v

    #load default params
    with open("calibrator/default-params/wheat.json") as _:
        species_params = json.load(_)
    with open("calibrator/default-params/winter-wheat.json") as _:
        cultivar_params = json.load(_)
    with open("calibrator/default-params/wheat-residue.json") as _:
        residue_params = json.load(_)

    #create crop object
    calibrated_custom_crop = {
        "is-winter-crop": True,
        "cropParams": {
            "species": species_params,
            "cultivar": {
                "=": cultivar_params,
                "StageTemperatureSum": [
                    [
                        150, 
                        584, 
                        526, 
                        168, 
                        389, 
                        25
                    ], 
                    "°C d"
                ]
            }
        },
        "residueParams": residue_params
    }

    default_custom_crop = {
        "is-winter-crop": True,
        "cropParams": {
            "species": species_params,
            "cultivar": cultivar_params
        },
        "residueParams": residue_params
    }

    server = {
        "producer": {
            "server": config["server"],
            "port": config["prod-port"],
            "nd-port": config["nd-port"]
        },
        "consumer": {
            "server": config["server"],
            "port": config["cons-port"],
            "nd-port": config["nd-port"]
        }
    }


    custom_crop = default_custom_crop #calibrated_custom_crop

    id_best, vals_params = start_calibration(custom_crop=custom_crop, server=server)
    with open("best_calibrations.csv", "ab") as _:
        writer = csv.writer(_)
        row = []
        row.append(str(run_id))
        row.append(str(id_best))
        row.append(str(vals_params))
        writer.writerow(row)
                

if __name__ == "__main__":
    main()
