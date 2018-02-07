import csv
import numpy as np
import pandas as pd
from datetime import date, timedelta
from pyproj import Proj, transform
from scipy.interpolate import NearestNDInterpolator

def template_ws(ref_crop, s_date, h_date):
    return {
        "worksteps": [
            {"date": s_date, "type": "Sowing", "crop": ref_crop},
            {"date": h_date, "type": "Harvest"}
        ]
    }
def template_out(out_date, out_var="Stage"):
    return out_date, [out_var]

def year_doy_to_isodate(year, doy):
    return (date(year, 1, 1) + timedelta(days=doy-1)).isoformat()


def create_rotation_events_gk5(path_to_station_csv_file, path_to_phenology_csv_file, wgs84, gk5, ref_crop="CROP_PLACEHOLDER"):
    #read station info
    DWD_stations = pd.read_csv(path_to_station_csv_file)
    DWD_stations_dict = DWD_stations.set_index('Stations_id').T.to_dict('dict')

    #read mgt and pheno info
    DWD_pheno = pd.read_csv(path_to_phenology_csv_file)

    # iterate over available stations
    DWD_stations = set(DWD_pheno["Stations_id"])
    stat_infos = []
    for station_id in DWD_stations:
        print("processing mgt and pheno of station " + str(station_id))
        station_cp_rot = []
        events = []
        station_pheno = DWD_pheno.loc[(DWD_pheno["Stations_id"] == station_id)]
        #identify available crop cycles
        cm_counts = set(station_pheno["CM-counter"])
        for cm in cm_counts:
            #get sowing and harvest to build crop rotation
            cm_pheno = station_pheno.loc[(station_pheno["CM-counter"] == cm)]
            
            s_year = int(cm_pheno.loc[(cm_pheno[" Phase_id"] == 10), " Referenzjahr"])
            s_doy = int(cm_pheno.loc[(cm_pheno[" Phase_id"] == 10), " Jultag"])
            s_date = year_doy_to_isodate(s_year, s_doy)
            
            h_year = int(cm_pheno.loc[(cm_pheno[" Phase_id"] == 24), " Referenzjahr"])
            h_doy = int(cm_pheno.loc[(cm_pheno[" Phase_id"] == 24), " Jultag"])
            h_date = year_doy_to_isodate(h_year, h_doy)

            my_cm = template_ws(ref_crop, s_date, h_date)
            station_cp_rot.append(my_cm)

            #build event structure
            for index, record in cm_pheno.iterrows():
                if record[" Phase_id"] in [10, 24]:
                    continue #skip sowing and harvest
                r_year = int(record[" Referenzjahr"])
                r_doy = int(record[" Jultag"])
                r_date = year_doy_to_isodate(r_year, r_doy)
                events.append(template_out(r_date, "Stage")[0])
                events.append(template_out(r_date, "Stage")[1])
        
        #store georeferenced info

        lat = float(DWD_stations_dict[station_id]["geograph.Breite"])
        lon = float(DWD_stations_dict[station_id]["geograph.Laenge"])
        r_gk5, h_gk5 = transform(wgs84, gk5, lon, lat)
            
        stat_infos.append({
            "rotation": station_cp_rot,
            "events": events,
            "station_id": station_id,
            "gk5_rh": (r_gk5, h_gk5)
        })

    #print len(points)
    return stat_infos


def DWD_to_MONICA_pheno(DWD_phase, crop_id):
    #keys: crop_id, DWD_phase; vals: MONICA_phase
    crop_pheno = {
        "WW":{
            12: 2,
            15: 3,
            18: 3,
            19: 5,
            21: 5

        }
    }
    return crop_pheno[crop_id][DWD_phase]


def read_stations_obspheno(path_to_station_csv_file, path_to_phenology_csv_file, crop_id):
    #read station info
    DWD_stations = pd.read_csv(path_to_station_csv_file)
    DWD_stations_dict = DWD_stations.set_index('Stations_id').T.to_dict('dict')

    #read mgt and pheno info
    DWD_pheno = pd.read_csv(path_to_phenology_csv_file)

    # iterate over available stations
    DWD_stations = set(DWD_pheno["Stations_id"])
    pheno_info = {}
    for station_id in DWD_stations:
        pheno_obs = []
        print("processing pheno of station " + str(station_id))
        station_pheno = DWD_pheno.loc[(DWD_pheno["Stations_id"] == station_id)]
        #identify available crop cycles
        cm_counts = set(station_pheno["CM-counter"])
        for cm in cm_counts:
            cm_pheno = station_pheno.loc[(station_pheno["CM-counter"] == cm)]            
            for index, record in cm_pheno.iterrows():
                if record[" Phase_id"] in [10, 24]:
                    continue #skip sowing and harvest
                r_year = int(record[" Referenzjahr"])
                r_doy = int(record[" Jultag"])
                r_date = year_doy_to_isodate(r_year, r_doy)
                pheno_obs.append((r_date, DWD_to_MONICA_pheno(record[" Phase_id"], crop_id)))

        pheno_info[station_id] = pheno_obs

    return pheno_info



if __name__ == "__main__":
    wgs84 = Proj(init="epsg:4326")
    gk5 = Proj(init="epsg:31469")

    rotation_events_gk5_interpolate = create_rotation_events_gk5_interpolator("subset_1995-2012/DWD_1995_2012_obs_phases_WW.csv", wgs84, gk5, "WW")

    print("done!")