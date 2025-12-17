#!/home/treynolds/mambaforge/envs/pypeit_dev/bin/python3

import os
import glob
from argparse import ArgumentParser
import subprocess
from astropy.io import fits
import json
import shutil
import logging

pypeit_path = "/home/treynolds/mambaforge/envs/pypeit_dev"
path = "/home/treynolds/data/LCO_spectra_reduction/scripts"

parser = ArgumentParser(
    description="Redoes the flux calibration"
)

parser.add_argument(
    "--input_folder",
    default=None,
    help="Path where the data was downloaded and sorted",
)

parser.add_argument(
    "--sensfunc",
    default=None,
    help="Provide the new sensfunc file",
)

parser.add_argument(
    "--requested_telescope",
    default=None,
    help="Specific LCO north or south",
)

args = parser.parse_args()
input_folder = args.input_folder
sensfunc = args.sensfunc
requested_telescope = args.requested_telescope

sensfunc_file = os.getcwd()+"/"+sensfunc

HERE = os.getcwd()
if os.path.isabs(input_folder):
    input_folder_abs_path = input_folder
else:
    input_folder_abs_path = os.path.join(HERE, input_folder)
all_folders = glob.glob(input_folder_abs_path+"/*")
#print(all_folders)
with open(input_folder_abs_path + "/data_summary_all.json", 'r') as f:
    data_summary_all = json.load(f)
reduction_status_dictionary = {}

for folder in all_folders :
    date = folder.split("/")[-1]
    # Currently, I am excluding any random files in the folder by hand
    # do something smarter than this.
    if date == "data_summary_all.json":
        continue
    if date == "logging_reduction.log":
        continue
    if date == "reduction_diags_all.json":
        continue
    with open(folder + '/reduction_diags_'+ date + '.json', 'r') as f: 
            reduction_status_dictionary[date] = json.load(f)
            
    telescope = data_summary_all[date]["telescope"]
    if telescope == "LCO_north":
        pypeit_spectrograph = "lco_floyds_north"
    elif telescope == "LCO_south" :
        pypeit_spectrograph = "lco_floyds_south"
    if telescope != requested_telescope :
        print(date, "is not the requested telescope, continuing")
        continue
    else :
        print(date, "is requested telescope, redoing flux calibration")

    os.chdir(folder+"/sci")

    # edit the .flux file
    flux_file = glob.glob(folder + "/sci/*.flux")

    new_flux_file_list = []
    with open(flux_file[0],"r") as file:
        flux_file_lines = file.readlines()
        print(flux_file_lines)
        for index,line in enumerate(flux_file_lines) :
            #print(line,index)
            if index == 12 :
                print(sensfunc)
                new_line = "|".join(line.split("|")[:-1] + [" " + sensfunc])
                new_flux_file_list.append(new_line)
            else :
                new_flux_file_list.append(line)
    #print(new_flux_file_list)
    new_flux_file = flux_file[0][:-5]+"_new.flux"
    with open(new_flux_file,"w") as file:
        #print("\n".join(new_flux_file_list))
        file.writelines("\n".join(new_flux_file_list))
        
    shutil.copy(sensfunc_file,
                folder + "/sci/Science/" + sensfunc)
    print("Applying the new sensitivity function...")
    subprocess.call([pypeit_path + "/bin/pypeit_flux_calib",
                 folder + "/sci/" + date + "_new.flux",
                 "-v","0", # turn the verbosity to none.
                     ])
   # break














        

