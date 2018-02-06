from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import os
import spotpy
import spotpy_setup_MONICA
import csv


def make_lambda(excel):
    return lambda v, p: eval(excel)

def start_calibration(custom_crop, server):

    #read params to be calibrated
    params = []
    with open('calibrator/calibratethese.csv') as paramscsv:
        dialect = csv.Sniffer().sniff(paramscsv.read(), delimiters=';,\t')
        paramscsv.seek(0)
        reader = csv.reader(paramscsv, dialect)
        next(reader, None)  # skip the header
        for row in reader:
            p={}
            p["name"] = row[0]
            p["array"] = row[1]
            p["low"] = float(row[2])
            p["high"] = float(row[3])
            p["stepsize"] = float(row[4])
            p["optguess"] = float(row[5])
            p["minbound"] = float(row[6])
            p["maxbound"] = float(row[7])
            if len(row) == 9 and row[8] != "":
                p["derive_function"] = make_lambda(row[8])
            params.append(p)

    spot_setup = spotpy_setup_MONICA.spot_setup(params, custom_crop, server)
    rep = 10

    sampler = spotpy.algorithms.sceua(spot_setup, dbname='SCEUA', dbformat='ram')
    #sampler.sample(rep, ngs=len(params)+1, kstop=2*len(params)+1, peps=0.0001, pcento=0.00001)
    sampler.sample(rep, ngs=len(params)+1, kstop=2)

    id_best = sampler.status.bestrep
    best_params = sampler.status.params

    print(str(id_best))

    print("finished!")

    return id_best, best_params


#if __name__ == "__main__":
#    start_calibration()
