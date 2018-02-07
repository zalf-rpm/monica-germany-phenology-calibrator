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

import sys
#sys.path.insert(0, "C:\\Users\\berg.ZALF-AD\\GitHub\\monica\\project-files\\Win32\\Release")
#sys.path.insert(0, "C:\\Users\\berg.ZALF-AD\\GitHub\\monica\\src\\python")
#sys.path.insert(0, "C:\\Program Files (x86)\\MONICA")
#print sys.path

import gc
import csv
import types
import os
import json
from datetime import datetime
from collections import defaultdict, OrderedDict
import numpy as np

import zmq
#print "pyzmq version: ", zmq.pyzmq_version(), " zmq version: ", zmq.zmq_version()

import monica_io
#print "path to monica_io: ", monica_io.__file__

LOCAL_RUN = False

PATHS = {
    "stella": {
        "local-path-to-output-dir": "out/"
    },
    "berg-lc": {
        "local-path-to-output-dir": "out/"
    },
    "berg-xps15": {
        "local-path-to-output-dir": "out/"
    }
}


def create_output(result):
    "create output structure for single run"

    cm_count_to_vals = defaultdict(dict)
    if len(result.get("data", [])) > 0 and len(result["data"][0].get("results", [])) > 0:

        for data in result.get("data", []):
            results = data.get("results", [])
            oids = data.get("outputIds", [])

            #skip empty results, e.g. when event condition haven't been met
            if len(results) == 0:
                continue

            assert len(oids) == len(results)
            for kkk in range(0, len(results[0])):
                vals = {}

                for iii in range(0, len(oids)):
                    oid = oids[iii]
                    val = results[iii][kkk]

                    name = oid["name"] if len(oid["displayName"]) == 0 else oid["displayName"]

                    if isinstance(val, types.ListType):
                        for val_ in val:
                            vals[name] = val_
                    else:
                        vals[name] = val

                if "CM-count" not in vals:
                    print "Missing CM-count in result section. Skipping results section."
                    continue

                cm_count_to_vals[vals["CM-count"]].update(vals)

    for cmc in sorted(cm_count_to_vals.keys()):
        if "harvest-doy" not in cm_count_to_vals[cmc]:
            del cm_count_to_vals[cmc] 

    return cm_count_to_vals


def write_row_to_grids(row_col_data, row, ncols, header, path_to_output_dir):
    "write grids row by row"

    if not hasattr(write_row_to_grids, "nodata_row_count"):
        write_row_to_grids.nodata_row_count = 0
        write_row_to_grids.list_of_output_files = []

    make_dict_nparr = lambda: defaultdict(lambda: np.full((ncols,), -9999, dtype=np.float))

    output_grids = {
        "yield": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 2}#,
#        "biom-final": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 2},
#        "crop-sum-precip": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
#        "crop-max-LAI": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 2},
#        "crop-avg-transpiration-deficit": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 2},
        #"avg-30cm-sand": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 2},
        #"avg-30cm-clay": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 2},
        #"avg-30cm-silt": {"data" : make_dict_nparr(), "cast-to": "int"},
#        "maturity-doy": {"data" : make_dict_nparr(), "cast-to": "int"},
#        "harvest-doy": {"data" : make_dict_nparr(), "cast-to": "int"},
#        "sowing-doy": {"data" : make_dict_nparr(), "cast-to": "int"},
#        "at-harvest-relative-total-development": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 2},
#        "doy90-to-harvest-sum-precip": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
#        "anthesis-doy": {"data" : make_dict_nparr(), "cast-to": "int"}#,
#        "yearly-avg-tavg": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
#        "yearly-sum-precip": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
#        "crop-avg-tavg": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
#        "crop-sum-precip": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
#        "crop-sum-nfert": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
#        "yearly-sum-nleach": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
    }

    cmc_to_crop = {}

    #is_no_data_row = True
    # skip this part if we write just a nodata line
    if row in row_col_data:
        no_data_cols = 0
        for col in xrange(0, ncols):
            if col in row_col_data[row]:
                rcd_val = row_col_data[row][col]
                if rcd_val == -9999:
                    no_data_cols += 1
                    continue
                else:
                    cmc_and_year_to_vals = defaultdict(lambda: defaultdict(list))
                    for cell_data in rcd_val:
                        for cm_count, data in cell_data.iteritems():
                            for key, val in output_grids.iteritems():
                                if cm_count not in cmc_to_crop:
                                    cmc_to_crop[cm_count] = data["Crop"]

                                if key in data:
                                    cmc_and_year_to_vals[(cm_count, data["Year"])][key].append(data[key])
                                else:
                                    cmc_and_year_to_vals[(cm_count, data["Year"])][key] #just make sure at least an empty list is in there

                    for (cm_count, year), key_to_vals in cmc_and_year_to_vals.iteritems():
                        for key, vals in key_to_vals.iteritems():
                            output_vals = output_grids[key]["data"]
                            if len(vals) > 0:
                                output_vals[(cm_count, year)][col] = sum(vals) / len(vals)
                            else:
                                output_vals[(cm_count, year)][col] = -9999
                                #no_data_cols += 1

        is_no_data_row = no_data_cols == ncols
        if is_no_data_row:
            write_row_to_grids.nodata_row_count += 1

    def write_nodata_rows(file_):
        for _ in range(write_row_to_grids.nodata_row_count):
            rowstr = " ".join(["-9999" for __ in range(ncols)])
            file_.write(rowstr +  "\n")

    for key, y2d_ in output_grids.iteritems():

        y2d = y2d_["data"]
        cast_to = y2d_["cast-to"]
        digits = y2d_.get("digits", 0)
        if cast_to == "int":
            mold = lambda x: str(int(x))
        else:
            mold = lambda x: str(round(x, digits))

        for (cm_count, year), row_arr in y2d.iteritems():

            crop = cmc_to_crop[cm_count]    
            crop = crop.replace("/", "").replace(" ", "")
            path_to_file = path_to_output_dir + crop + "_" + key + "_" + str(year) + "_" + str(cm_count) + ".asc"

            if not os.path.isfile(path_to_file):
                with open(path_to_file, "w") as _:
                    _.write(header)
                    write_row_to_grids.list_of_output_files.append(path_to_file)

            with open(path_to_file, "a") as file_:
                write_nodata_rows(file_)
                rowstr = " ".join(["-9999" if int(x) == -9999 else mold(x) for x in row_arr])
                file_.write(rowstr +  "\n")

    # if we're at the end of the output and just empty lines are left, then they won't be written in the
    # above manner because there won't be any rows with data where they could be written before
    # so add no-data rows simply to all files we've written to before
    if is_no_data_row and write_row_to_grids.list_of_output_files:
        for path_to_file in write_row_to_grids.list_of_output_files:
            with open(path_to_file, "a") as file_:
                write_nodata_rows(file_)
        write_row_to_grids.nodata_row_count = 0
    
    # clear the no-data row count when no-data rows have been written before a data row
    if not is_no_data_row:
        write_row_to_grids.nodata_row_count = 0

    if row in row_col_data:
        del row_col_data[row]


def run_consumer(server = {"server": None, "port": None, "nd-port": None}, path_to_output_dir=None):
    "collect data from workers"

    config = {
        "user": "stella",
        "port": server["port"] if server["port"] else "77773",
        "no-data-port": server["nd-port"] if server["nd-port"] else "5555",
        "server": server["server"] if server["server"] else "localhost"
    }
    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            k,v = arg.split("=")
            if k in config:
                config[k] = v

    paths = PATHS[config["user"]]
    if path_to_output_dir:
        paths["local-path-to-output-dir"] = path_to_output_dir

    try:
        os.makedirs(paths["local-path-to-output-dir"])
    except:
        pass

    data = defaultdict(list)

    print "consumer config:", config

    received_env_count = 1
    context = zmq.Context()
    socket = context.socket(zmq.PULL)
    
    #connect producer and consumer directly
    prod_cons_socket = context.socket(zmq.PULL)
    prod_cons_socket.connect("tcp://localhost:" + config["no-data-port"])
    expected_stats = None
    expected_stats = prod_cons_socket.recv_json(encoding="latin-1")


    if LOCAL_RUN:
        socket.connect("tcp://localhost:" + config["port"])
    else:
        socket.connect("tcp://" + config["server"] + ":" + config["port"])
    socket.RCVTIMEO = 1000
    leave = False
    write_calibration_output_files = True
    

    #for calib out
    to_write = []
    header = ["Station", "Date", "Stage"]
    to_write.append(header)

    while not leave:

        try:
            result = socket.recv_json(encoding="latin-1")
        except:
            continue

        if write_calibration_output_files:

            if result.get("type", "") in ["jobs-per-cell", "no-data", "target-grid-metadata"]:
                print "ignoring", result.get("type", "")
                continue

            print "received work result ", received_env_count, " customId: ", result.get("customId", "")

            custom_id = result["customId"]
            station = custom_id["station_id"]
            
            #save expected station and leave when they are all simulated            
            expected_stats.remove(station)
            print ("missing stations: " + str(expected_stats))
            if len(expected_stats) == 0:
                leave = True
            
            
            for data_ in result.get("data", []):
                results = data_.get("results", [])
                orig_spec = data_.get("origSpec", "")
                output_ids = data_.get("outputIds", [])
                row = []
                row.append(str(station))
                row.append(orig_spec)
                row.append(str(results[0][0]))
                to_write.append(row)
            

            received_env_count = received_env_count + 1
        
        else: #normal files

            if result.get("type", "") in ["jobs-per-cell", "no-data", "target-grid-metadata"]:
                print "ignoring", result.get("type", "")
                continue

            print "received work result ", received_env_count, " customId: ", result.get("customId", "")

            custom_id = result["customId"]
            ci_parts = custom_id.split("|")
            station = int(ci_parts[0])
            
            
            #with open("out/out-" + str(i) + ".csv", 'wb') as _:
            with open("out/out-" + custom_id.replace("|", "_") + ".csv", 'wb') as _:
                writer = csv.writer(_, delimiter=",")

                for data_ in result.get("data", []):
                    results = data_.get("results", [])
                    orig_spec = data_.get("origSpec", "")
                    output_ids = data_.get("outputIds", [])

                    if len(results) > 0:
                        writer.writerow([orig_spec.replace("\"", "")])
                        for row in monica_io.write_output_header_rows(output_ids,
                                                                      include_header_row=True,
                                                                      include_units_row=True,
                                                                      include_time_agg=False):
                            writer.writerow(row)

                        for row in monica_io.write_output(output_ids, results):
                            writer.writerow(row)

                    writer.writerow([])

            received_env_count = received_env_count + 1
    
    if write_calibration_output_files:
        print("writing out file for pheno cal...")
        with open("out/calib_out.csv", "wb") as out_file:
            writer = csv.writer(out_file)
            for row in to_write:
                writer.writerow(row)
                

    print "exiting run_consumer()"
    #debug_file.close()

if __name__ == "__main__":
    run_consumer()
#main()


