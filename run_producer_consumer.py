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

import os
import json
import csv
import copy
from collections import defaultdict
import sys
from glob import glob
#print sys.path

from multiprocessing import Process
import threading
from threading import Thread

import run_work_producer
import run_work_consumer

class FuncThread(threading.Thread):
    def __init__(self, target, *args):
        self._target = target
        self._args = args
        threading.Thread.__init__(self)
 
    def run(self):
        self._target(*self._args)


def prod_cons_calib(custom_crop, server):

    setup = {
        "crop": "WW",
        "groundwater-level": False,
        "impenetrable-layer": False,
        "elevation": True,
        "latitude": True,
        "slope": True,
        "fertilization": False,
        "NitrogenResponseOn": False,
        "irrigation": False,
        "WaterDeficitResponseOn": True,
        "LeafExtensionModifier": False,
        "EmergenceMoistureControlOn": False,
        "EmergenceFloodingControlOn": False
    }
    
    producer = Process(target=run_work_producer.run_producer, args=(setup, custom_crop, server["producer"]))
    consumer = Process(target=run_work_consumer.run_consumer, args=(server["consumer"]))
    #producer = FuncThread(rwp.run_producer, setup, custom_crop)
    #consumer = FuncThread(rwc.run_consumer, path_to_grids_output, True)
    producer.daemon = True
    consumer.daemon = True
    producer.start()
    consumer.start()
    producer.join()
    consumer.join()



    return year_to_lk_to_value
