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
#print sys.path
import zmq
#print "pyzmq version: ", zmq.pyzmq_version(), " zmq version: ", zmq.zmq_version()

import sqlite3
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from pyproj import Proj, transform

import monica_io
import soil_io
import ascii_io
import monica_germany_utils

LOCAL_RUN = False

PATHS = {
    "stella": {
        "include-file-base-path": "C:/Users/stella/Documents/GitHub",
        "path-to-soil-dir": "Z:/data/soil/buek1000/brd/",
        "path-to-climate-csvs-dir": "Z:/data/climate/dwd/csvs/germany/",
        "path-to-data-dir": "Z:/data/",
        "path-to-projects-dir": "Z:/projects/",
        "archive-path-to-climate-csvs-dir": "/archiv-daten/md/data/climate/dwd/csvs/germany/"
    },
    "berg-lc": {
        "include-file-base-path": "C:/Users/berg.ZALF-AD/GitHub",
        "path-to-soil-dir": "N:/soil/buek1000/brd/",
        "path-to-climate-csvs-dir": "N:/climate/dwd/csvs/germany/",
        #"path-to-climate-csvs-dir": "N:/climate/isimip/csvs/germany/",
        "path-to-data-dir": "N:/",
        "path-to-projects-dir": "P:/",
        "archive-path-to-climate-csvs-dir": "/archiv-daten/md/data/climate/dwd/csvs/germany/"
        #"archive-path-to-climate-csvs-dir": "/archiv-daten/md/data/climate/isimip/csvs/germany/"
    },
    "berg-xps15": {
        "include-file-base-path": "C:/Users/berg.ZALF-AD/GitHub",
        "path-to-soil-dir": "D:/soil/buek1000/brd/",
        "path-to-climate-csvs-dir": "D:/climate/dwd/csvs/germany/",
        #"path-to-climate-csvs-dir": "N:/climate/isimip/csvs/germany/",
        "path-to-data-dir": "N:/",
        "path-to-projects-dir": "P:/",
        "archive-path-to-climate-csvs-dir": "/archiv-daten/md/data/climate/dwd/csvs/germany/"
        #"archive-path-to-climate-csvs-dir": "/archiv-daten/md/data/climate/isimip/csvs/germany/"
    }
}

USER = "stella"

paths = PATHS[USER]

def preprocess_data():
    
    def read_header(path_to_ascii_grid_file):
        "read metadata from esri ascii grid file"
        metadata = {}
        header_str = ""
        with open(path_to_ascii_grid_file) as _:
            for i in range(0, 6):
                line = _.readline()
                header_str += line
                sline = [x for x in line.split() if len(x) > 0]
                if len(sline) > 1:
                    metadata[sline[0].strip().lower()] = float(sline[1].strip())
        return metadata, header_str


    wgs84 = Proj(init="epsg:4326")
    #gk3 = Proj(init="epsg:3396")
    gk5 = Proj(init="epsg:31469")

    def create_ascii_grid_interpolator(arr, meta, ignore_nodata=True):
        "read an ascii grid into a map, without the no-data values"

        rows, cols = arr.shape

        cellsize = int(meta["cellsize"])
        xll = int(meta["xllcorner"])
        yll = int(meta["yllcorner"])
        nodata_value = meta["nodata_value"]

        xll_center = xll + cellsize // 2
        yll_center = yll + cellsize // 2
        yul_center = yll_center + (rows - 1)*cellsize

        points = []
        values = []

        for row in range(rows):
            for col in range(cols):
                value = arr[row, col]
                if ignore_nodata and value == nodata_value:
                    continue
                r = xll_center + col * cellsize
                h = yul_center - row * cellsize
                points.append([r, h])
                values.append(value)

        return NearestNDInterpolator(np.array(points), np.array(values))
    
    print("reading DEM info...")
    path_to_dem_grid = paths["path-to-data-dir"] + "/germany/dem_1000_gk5.asc"
    dem_metadata, _ = read_header(path_to_dem_grid)
    dem_grid = np.loadtxt(path_to_dem_grid, dtype=int, skiprows=6)
    dem_gk5_interpolate = create_ascii_grid_interpolator(dem_grid, dem_metadata)
    
    print("reading slope info...")
    path_to_slope_grid = paths["path-to-data-dir"] + "/germany/slope_1000_gk5.asc"
    slope_metadata, _ = read_header(path_to_slope_grid)
    slope_grid = np.loadtxt(path_to_slope_grid, dtype=float, skiprows=6)
    slope_gk5_interpolate = create_ascii_grid_interpolator(slope_grid, slope_metadata)

    print("reading soil info...")
    path_to_soil_grid = paths["path-to-data-dir"] + "/germany/buek1000_1000_gk5.asc"
    soil_meta, _ = read_header(path_to_soil_grid)
    soil_grid = np.loadtxt(path_to_soil_grid, dtype=int, skiprows=6)
    soil_gk5_interpolate = create_ascii_grid_interpolator(soil_grid, soil_meta)
    

    cdict = {}
    def create_climate_gk5_interpolator_from_json_file(path_to_latlon_to_rowcol_file, wgs84, gk5):
        "create interpolator from json list of lat/lon to row/col mappings"
        with open(path_to_latlon_to_rowcol_file) as _:
            points = []
            values = []

            for latlon, rowcol in json.load(_):
                row, col = rowcol
                clat, clon = latlon
                try:
                    cr_gk5, ch_gk5 = transform(wgs84, gk5, clon, clat)
                    cdict[(row, col)] = (round(clat, 4), round(clon, 4))
                    points.append([cr_gk5, ch_gk5])
                    values.append((row, col))
                    #print "row:", row, "col:", col, "clat:", clat, "clon:", clon, "h:", h, "r:", r, "val:", values[i]
                except:
                    continue

            return NearestNDInterpolator(np.array(points), np.array(values))
    print("interpolating climate...")
    climate_gk5_interpolate = create_climate_gk5_interpolator_from_json_file(paths["path-to-climate-csvs-dir"] + "../latlon-to-rowcol.json", wgs84, gk5)

    stat_infos = monica_germany_utils.create_rotation_events_gk5(paths["path-to-projects-dir"] + "monica-germany/DWD_stations.csv",
    paths["path-to-projects-dir"] + "monica-germany/DWD_1995_2012_obs_phases_WW.csv", wgs84, gk5)

    preprocessed_data = {
        "dem_gk5_interpolate": dem_gk5_interpolate,
        "slope_gk5_interpolate": slope_gk5_interpolate,
        "soil_gk5_interpolate": soil_gk5_interpolate,
        "cdict": cdict,
        "climate_gk5_interpolate": climate_gk5_interpolate,
        "stat_infos": stat_infos
    }

    return preprocessed_data


def run_producer(setup, custom_crop, server = {"server": None, "port": None, "nd-port": None}, preprocessed_data=None):
    "main"

    config = {
        "user": USER,
        "port": server["port"] if server["port"] else "66663",
        "no-data-port": server["nd-port"] if server["nd-port"] else "5555",
        "server": server["server"] if server["server"] else "localhost"
    }

    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            k, v = arg.split("=")
            if k in config:
                config[k] = v

    context = zmq.Context()
    socket = context.socket(zmq.PUSH)
    prod_cons_socket = context.socket(zmq.PUSH)

    soil_db_con = sqlite3.connect(paths["path-to-data-dir"] + "germany/buek1000.sqlite")

    #connect producer and consumer directly
    prod_cons_socket.bind("tcp://*:" + str(config["no-data-port"]))

    if LOCAL_RUN:
        socket.connect("tcp://localhost:" + str(config["port"]))
    else:
        socket.connect("tcp://" + config["server"] + ":" + str(config["port"]))

    with open("sim.json") as _:
        sim_json = json.load(_)
        sim_json["include-file-base-path"] = paths["include-file-base-path"]
    
    with open("site.json") as _:
        site_json = json.load(_)

    with open("crop.json") as _:
        crop_json = json.load(_)

    def get_value(list_or_value):
        return list_or_value[0] if isinstance(list_or_value, list) else list_or_value
    
    sent_env_count = 1
    start_time = time.clock()

    if preprocessed_data == None:
        preprocessed_data = preprocess_data()    
    
    dem_gk5_interpolate = preprocessed_data["dem_gk5_interpolate"]
    slope_gk5_interpolate = preprocessed_data["slope_gk5_interpolate"]
    soil_gk5_interpolate = preprocessed_data["soil_gk5_interpolate"]
    cdict = preprocessed_data["cdict"]
    climate_gk5_interpolate = preprocessed_data["climate_gk5_interpolate"]
    stat_infos = preprocessed_data["stat_infos"]

    #put custom crop in place
    for stat_info in stat_infos:
        for cm in stat_info["rotation"]:
            for ws in cm["worksteps"]:
                if "crop" in ws.keys():
                    ws["crop"] = custom_crop

    #send stations list to the consumer
    sim_stats = []
    for stat_info in stat_infos:
        sim_stats.append(int(stat_info["station_id"]))
    prod_cons_socket.send_json(sim_stats)

    for stat_info in stat_infos:

        station_id = stat_info["station_id"]
        r_gk5, h_gk5 = stat_info["gk5_rh"]

        env_template = monica_io.create_env_json_from_json_config({
            "crop": crop_json,
            "site": site_json,
            "sim": sim_json,
            "climate": ""
        })
        env_template["cropRotation"] = stat_info["rotation"]
        env_template["events"] = stat_info["events"]
        
        soil_id = soil_gk5_interpolate(r_gk5, h_gk5)
        if soil_id == -9999:
            continue
        if soil_id < 1 or soil_id > 71:
            #print "row/col:", row, "/", col, "has unknown soil_id:", soil_id
            #unknown_soil_ids.add(soil_id)
            continue
        
        crow, ccol = climate_gk5_interpolate(r_gk5, h_gk5)

        height_nn = dem_gk5_interpolate(r_gk5, h_gk5)
        slope = slope_gk5_interpolate(r_gk5, h_gk5)
        
        clat, clon = cdict[(crow, ccol)]
        #slon, slat = transform(gk5, wgs84, r, h)
        #print "srow:", srow, "scol:", scol, "h:", h, "r:", r, " inter:", inter, "crow:", crow, "ccol:", ccol, "slat:", slat, "slon:", slon, "clat:", clat, "clon:", clon

        env_template["params"]["userCropParameters"]["__enable_T_response_leaf_expansion__"] = setup["LeafExtensionModifier"]

        # set soil-profile
        sp_json = soil_io.soil_parameters(soil_db_con, soil_id)
        soil_profile = monica_io.find_and_replace_references(sp_json, sp_json)["result"]
        #print "soil:", soil_profile
        env_template["params"]["siteParameters"]["SoilProfileParameters"] = soil_profile

        # setting groundwater level
        if setup["groundwater-level"]:
            groundwaterlevel = 20
            layer_depth = 0
            for layer in soil_profile:
                if layer.get("is_in_groundwater", False):
                    groundwaterlevel = layer_depth
                    #print "setting groundwaterlevel of soil_id:", str(soil_id), "to", groundwaterlevel, "m"
                    break
                layer_depth += get_value(layer["Thickness"])
            env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepthMonth"] = 3
            env_template["params"]["userEnvironmentParameters"]["MinGroundwaterDepth"] = [max(0, groundwaterlevel - 0.2) , "m"]
            env_template["params"]["userEnvironmentParameters"]["MaxGroundwaterDepth"] = [groundwaterlevel + 0.2, "m"]
            
        # setting impenetrable layer
        if setup["impenetrable-layer"]:
            impenetrable_layer_depth = get_value(env_template["params"]["userEnvironmentParameters"]["LeachingDepth"])
            layer_depth = 0
            for layer in soil_profile:
                if layer.get("is_impenetrable", False):
                    impenetrable_layer_depth = layer_depth
                    #print "setting leaching depth of soil_id:", str(soil_id), "to", impenetrable_layer_depth, "m"
                    break
                layer_depth += get_value(layer["Thickness"])
            env_template["params"]["userEnvironmentParameters"]["LeachingDepth"] = [impenetrable_layer_depth, "m"]
            env_template["params"]["siteParameters"]["ImpenetrableLayerDepth"] = [impenetrable_layer_depth, "m"]

        if setup["elevation"]:
            env_template["params"]["siteParameters"]["heightNN"] = height_nn

        if setup["slope"]:
            env_template["params"]["siteParameters"]["slope"] = slope / 100.0

        if setup["latitude"]:
            env_template["params"]["siteParameters"]["Latitude"] = clat

        env_template["params"]["simulationParameters"]["UseNMinMineralFertilisingMethod"] = setup["fertilization"]
        env_template["params"]["simulationParameters"]["UseAutomaticIrrigation"] = setup["irrigation"]

        env_template["params"]["simulationParameters"]["NitrogenResponseOn"] = setup["NitrogenResponseOn"]
        env_template["params"]["simulationParameters"]["WaterDeficitResponseOn"] = setup["WaterDeficitResponseOn"]
        env_template["params"]["simulationParameters"]["EmergenceMoistureControlOn"] = setup["EmergenceMoistureControlOn"]
        env_template["params"]["simulationParameters"]["EmergenceFloodingControlOn"] = setup["EmergenceFloodingControlOn"]

        env_template["csvViaHeaderOptions"] = sim_json["climate.csv-options"]

        if LOCAL_RUN:
            env_template["pathToClimateCSV"] = paths["path-to-climate-csvs-dir"] + "row-" + str(crow) + "/col-" + str(ccol) + ".csv"
        else:
            env_template["pathToClimateCSV"] = paths["archive-path-to-climate-csvs-dir"] + "row-" + str(crow) + "/col-" + str(ccol) + ".csv"

        env_template["customId"] = {
            "station_id": int(station_id)
        }

        #with open("envs/env-"+str(sent_env_count)+".json", "w") as _: 
        #    _.write(json.dumps(env))

        socket.send_json(env_template)
        print "sent env ", sent_env_count, " customId: ", env_template["customId"]
        #exit()
        sent_env_count += 1

    stop_time = time.clock()



    print "sending ", (sent_env_count-1), " envs took ", (stop_time - start_time), " seconds"
    #print "ran from ", start, "/", row_cols[start], " to ", end, "/", row_cols[end]
    print "exiting run_producer()"

if __name__ == "__main__":
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

    #load default params
    with open("calibrator/default-params/wheat.json") as _:
        species_params = json.load(_)
    with open("calibrator/default-params/winter-wheat.json") as _:
        cultivar_params = json.load(_)
    with open("calibrator/default-params/wheat-residue.json") as _:
        residue_params = json.load(_)

    #create crop object
    custom_crop = {
        "CROP_ID": "WW", 
        "is-winter-crop": True,
        "cropParams": {
            "species": species_params,
            "cultivar": cultivar_params
        },
        "residueParams": residue_params
    }

    run_producer(setup, custom_crop)